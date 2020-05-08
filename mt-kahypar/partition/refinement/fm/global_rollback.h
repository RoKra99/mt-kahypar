/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <atomic>
#include <mt-kahypar/partition/metrics.h>
#include <mt-kahypar/utils/timer.h>

#include "fm_commons.h"

#include <boost/dynamic_bitset.hpp>

#include "tbb/parallel_scan.h"

namespace mt_kahypar {

struct BalanceAndBestIndexScan {
  const PartitionedHypergraph& phg;
  const vec<Move>& moves;

  struct GainIndex {
    Gain gain = 0;          /** gain when using all (still valid) moves up to best_index */
    MoveID best_index = 0;  /** local ID of first move to revert */
  };

  std::shared_ptr< tbb::enumerable_thread_specific<GainIndex> > local_best;

  Gain gain_sum = 0;

  vec<HypernodeWeight> part_weights;
  HypernodeWeight max_part_weight;

  BalanceAndBestIndexScan(BalanceAndBestIndexScan& b, tbb::split) :
          phg(b.phg),
          moves(b.moves),
          local_best(b.local_best),
          part_weights(b.part_weights.size(), 0),
          max_part_weight(b.max_part_weight) { }

  BalanceAndBestIndexScan(const PartitionedHypergraph& phg, const vec<Move>& moves,
                          vec<HypernodeWeight>& part_weights, HypernodeWeight max_part_weight) :
          phg(phg),
          moves(moves),
          local_best(std::make_shared< tbb::enumerable_thread_specific<GainIndex> >()),
          part_weights(part_weights),
          max_part_weight(max_part_weight) { }


  void operator()(const tbb::blocked_range<MoveID>& r, tbb::pre_scan_tag ) {
    for (MoveID i = r.begin(); i < r.end(); ++i) {
      const Move& m = moves[i];
      if (m.gain != invalidGain) {  // skip locally reverted moves
        gain_sum += m.gain;
        part_weights[m.from] -= phg.nodeWeight(m.node);
        part_weights[m.to] += phg.nodeWeight(m.node);
      }
    }
  }

  // subranges a | b | c | d . assuming this ran pre_scan on c,
  // then lhs ran pre_scan on b and final_scan of this will be on d
  void reverse_join(BalanceAndBestIndexScan& lhs) {
    for (size_t i = 0; i < part_weights.size(); ++i) {
      part_weights[i] += lhs.part_weights[i];
    }
    gain_sum += lhs.gain_sum;
  }

  void operator()(const tbb::blocked_range<MoveID>& r, tbb::final_scan_tag ) {
    size_t overloaded = 0;
    for (size_t i = 0; i < part_weights.size(); ++i) {
      if (part_weights[i] > max_part_weight) {
        overloaded++;
      }
    }

    Gain best_gain_sum = 0;
    MoveID best_index = 0;
    for (MoveID i = r.begin(); i < r.end(); ++i) {
      const Move& m = moves[i];

      if (m.gain != invalidGain) {  // skip locally reverted moves
        const bool from_overloaded = part_weights[m.from] > max_part_weight;
        part_weights[m.from] -= phg.nodeWeight(m.node);
        if (from_overloaded && part_weights[m.from] <= max_part_weight) {
          overloaded--;
        }

        const bool to_overloaded = part_weights[m.to] > max_part_weight;
        part_weights[m.to] += phg.nodeWeight(m.node);
        if (!to_overloaded && part_weights[m.to] > max_part_weight) {
          overloaded++;
        }

        gain_sum += m.gain;
        if (overloaded == 0 /* in balance */ && gain_sum > best_gain_sum /* better solution */) {
          best_gain_sum = gain_sum;
          best_index = i + 1;
        }
      }
    }

    if (best_index != 0) {
      GainIndex& lb = local_best->local();
      if (best_gain_sum > lb.gain) {
        lb.gain = best_gain_sum;
        lb.best_index = best_index;
      }
    }
  }

  void assign(BalanceAndBestIndexScan& ) {

  }

  GainIndex finalize() {
    GainIndex res = { 0, 0 };
    for (const GainIndex& x : *local_best) {
      if (x.gain > res.gain || (x.gain == res.gain && x.best_index < res.best_index)) {
        res = x;
      }
    }
    return res;
  }
};


class GlobalRollback {
  static constexpr bool enable_heavy_assert = false;
public:
  explicit GlobalRollback(const Hypergraph& hg, PartitionID numParts) :
          numParts(numParts),
          remaining_original_pins(hg.numNonGraphEdges() * numParts),
          first_move_in(hg.numNonGraphEdges() * numParts),
          last_move_out(hg.numNonGraphEdges() * numParts)
  {
    resetStoredMoveIDs();
  }

  HyperedgeWeight revertToBestPrefix(PartitionedHypergraph& phg, FMSharedData& sharedData,
                                     vec<HypernodeWeight>& partWeights, HypernodeWeight maxPartWeight,
                                     bool parallel = true) {
    if (parallel) {
      return revertToBestPrefixParallel(phg, sharedData, partWeights, maxPartWeight);
    } else {
      return revertToBestPrefixSequential(phg, sharedData, partWeights, maxPartWeight);
    }
  }

  HyperedgeWeight revertToBestPrefixSequential(PartitionedHypergraph& phg, FMSharedData& sharedData,
                                     vec<HypernodeWeight>& , HypernodeWeight maxPartWeight) {


    GlobalMoveTracker& tracker = sharedData.moveTracker;
    const MoveID numMoves = tracker.numPerformedMoves();
    const vec<Move>& move_order = tracker.moveOrder;

    // revert all moves
    tbb::parallel_for(0U, numMoves, [&](const MoveID localMoveID) {
      const Move& m = move_order[localMoveID];
      if (tracker.isMoveStillValid(m)) {
        phg.changeNodePartFullUpdate(m.node, m.to, m.from, std::numeric_limits<HypernodeWeight>::max(), [] { });
      }
    });


    size_t num_unbalanced_slots = 0;

    size_t overloaded = 0;
    for (PartitionID i = 0; i < numParts; ++i) {
      if (phg.partWeight(i) > maxPartWeight) {
        overloaded++;
      }
    }

    // roll forward sequentially
    Gain best_gain = 0, gain_sum = 0;
    MoveID best_index = 0;
    for (MoveID localMoveID = 0; localMoveID < numMoves; ++localMoveID) {
      const Move& m = move_order[localMoveID];
      if (!tracker.isMoveStillValid(m)) continue;

      Gain gain = 0;
      for (HyperedgeID e : phg.incidentEdges(m.node)) {
        if (phg.pinCountInPart(e, m.from) == 1) {
          gain += phg.edgeWeight(e);
        }
        if (phg.pinCountInPart(e, m.to) == 0) {
          gain -= phg.edgeWeight(e);
        }
      }
      gain_sum += gain;

      const bool from_overloaded = phg.partWeight(m.from) > maxPartWeight, to_overloaded = phg.partWeight(m.to) > maxPartWeight;
      phg.changeNodePartFullUpdate(m.node, m.from, m.to, std::numeric_limits<HypernodeWeight>::max(), [] { });
      if (from_overloaded && phg.partWeight(m.from) <= maxPartWeight) {
        overloaded--;
      }
      if (!to_overloaded && phg.partWeight(m.to) > maxPartWeight) {
        overloaded++;
      }

      if (overloaded > 0) {
        num_unbalanced_slots++;
      }

      if (overloaded == 0 && gain_sum > best_gain) {
        best_index = localMoveID + 1;
      }
    }

    // revert rejected moves again
    tbb::parallel_for(best_index, numMoves, [&](const MoveID i) {
      const Move& m = move_order[i];
      if (tracker.isMoveStillValid(m)) {
        phg.changeNodePartFullUpdate(m.node, m.to, m.from, std::numeric_limits<HypernodeWeight>::max(), [] { });
      }
    });

    tbb::parallel_for(0U, numMoves, [&](const MoveID i) {
      phg.recomputeMoveFromBenefit(move_order[i].node);
    });

    tracker.reset();

    return best_gain;
  }

  HyperedgeWeight revertToBestPrefixParallel(PartitionedHypergraph& phg, FMSharedData& sharedData,
                                     vec<HypernodeWeight>& partWeights, HypernodeWeight maxPartWeight) {
    const MoveID numMoves = sharedData.moveTracker.numPerformedMoves();
    if (numMoves == 0) return 0;

    const vec<Move>& move_order = sharedData.moveTracker.moveOrder;
    utils::Timer& timer = utils::Timer::instance();

    timer.start_timer("recalculate_gains", "Recalculate Gains");
    recalculateGains(phg, sharedData);
    timer.stop_timer("recalculate_gains");

    timer.start_timer("find_best_prefix_and_balance", "Find Best Balanced Prefix");
    BalanceAndBestIndexScan s(phg, move_order, partWeights, maxPartWeight);
    tbb::parallel_scan(tbb::blocked_range<MoveID>(0, numMoves, 2500), s);
    BalanceAndBestIndexScan::GainIndex b = s.finalize();
    timer.stop_timer("find_best_prefix_and_balance");

    timer.start_timer("revert_and_rem_orig_pin_updates", "Revert Moves and apply updates");

    tbb::parallel_invoke([&] {
      // revert rejected moves
      tbb::parallel_for(b.best_index, numMoves, [&](const MoveID moveID) {
        const Move& m = move_order[moveID];
        if (sharedData.moveTracker.isMoveStillValid(m)) {
          phg.changeNodePartFullUpdate(m.node, m.to, m.from, std::numeric_limits<HypernodeWeight>::max(), []{/* do nothing */});
          for (HyperedgeID e : phg.incidentEdges(m.node)) {
            if (phg.edgeSize(e) > 2) {
              remaining_original_pins[phg.nonGraphEdgeID(e) * numParts + m.from].fetch_add(1, std::memory_order_relaxed);
            }
          }
        }
      });
    }, [&] {
      // apply updates to remaining original pins
      tbb::parallel_for(MoveID(0), b.best_index, [&](const MoveID moveID) {
        const Move& m = move_order[moveID];
        if (sharedData.moveTracker.isMoveStillValid(m)) {
          for (HyperedgeID e : phg.incidentEdges(move_order[moveID].node)) {
            if (phg.edgeSize(e) > 2) {
              remaining_original_pins[phg.nonGraphEdgeID(e) * numParts + m.to].fetch_add(1, std::memory_order_relaxed);
            }
          }
        }
      });
    } );

    timer.stop_timer("revert_and_rem_orig_pin_updates");

    timer.start_timer("recompute_move_from_benefits", "Recompute Move-From Benefits");
    // recompute moveFromBenefit values since they are potentially invalid
    // TODO think of a way to update them
    tbb::parallel_for(MoveID(0), numMoves, [&](MoveID localMoveID) {
      phg.recomputeMoveFromBenefit(move_order[localMoveID].node);
    });
    timer.stop_timer("recompute_move_from_benefits");

    if (sharedData.moveTracker.reset()) {
      resetStoredMoveIDs();
    }

    HEAVY_REFINEMENT_ASSERT([&] {
      for ( const HyperedgeID& he : phg.edges() ) {
        if (phg.edgeSize(he) == 2) {
          continue;
        }
        for ( PartitionID block = 0; block < numParts; ++block ) {
          if ( phg.pinCountInPart(he, block) != remaining_original_pins[phg.nonGraphEdgeID(he) * numParts + block] ) {
            LOG << "Pin Count In Part does not match remaining original pins" << V(he) << V(block)
                << V(phg.pinCountInPart(he, block)) << V(remaining_original_pins[phg.nonGraphEdgeID(he) * numParts + block]);
            return false;
          }
        }
      }
      return true;
    }());
    HEAVY_REFINEMENT_ASSERT(phg.checkTrackedPartitionInformation());
    return b.gain;
  }


  void recalculateGains(PartitionedHypergraph& phg, FMSharedData& sharedData) {
    GlobalMoveTracker& tracker = sharedData.moveTracker;
    vec<Move>& move_order = tracker.moveOrder;
    const MoveID numMoves = tracker.numPerformedMoves();
    const MoveID firstMoveID = tracker.firstMoveID;

    utils::Timer& timer = utils::Timer::instance();
    timer.start_timer("move_id_flagging", "Move Flagging");

    tbb::parallel_for(0U, numMoves, [&](MoveID localMoveID) {
      const Move& m = move_order[localMoveID];
      if (!sharedData.moveTracker.isMoveStillValid(m)) return;

      const MoveID globalMoveID = localMoveID + firstMoveID;
      for (HyperedgeID e_global : phg.incidentEdges(m.node)) {
        if (phg.edgeSize(e_global) > 2) {
          const HyperedgeID e = phg.nonGraphEdgeID(e_global);
          CAtomic<MoveID>& fmi = first_move_in[e * numParts + m.to];
          MoveID expected = fmi.load(std::memory_order_acq_rel);
          // first_move_in = min(first_move_in, this_move)
          while ((tracker.isMoveStale(expected) || expected > globalMoveID)
                 && !fmi.compare_exchange_weak(expected, globalMoveID, std::memory_order_acq_rel)) { }

          CAtomic<MoveID>& lmo = last_move_out[e * numParts + m.from];
          expected = lmo.load(std::memory_order_acq_rel);
          // last_move_out = max(last_move_out, this_move)
          while (expected < globalMoveID && !lmo.compare_exchange_weak(expected, globalMoveID, std::memory_order_acq_rel)) { }

          remaining_original_pins[e * numParts + m.from].fetch_sub(1, std::memory_order_relaxed);
        }
      }
    });

    timer.stop_timer("move_id_flagging");
    timer.start_timer("gain_recalculation", "Recalculate Gains");

    tbb::parallel_for(0U, numMoves, [&](MoveID localMoveID) {
      Move& m = move_order[localMoveID];
      if (!tracker.isMoveStillValid(m)) {
        return;
      }

      MoveID globalMoveID = firstMoveID + localMoveID;
      Gain gain = 0;
      for (HyperedgeID e_global : phg.incidentEdges(m.node)) {
        const HyperedgeWeight edgeWeight = phg.edgeWeight(e_global);

        if (phg.edgeSize(e_global) > 2) {
          const HyperedgeID e = phg.nonGraphEdgeID(e_global);

          const MoveID firstInFrom = firstMoveIn(e, m.from);
          if (remainingPinsFromBeginningOfMovePhase(e, m.from) == 0 && lastMoveOut(e, m.from) == globalMoveID
              && (firstInFrom > globalMoveID || firstInFrom < firstMoveID)) {
            gain += edgeWeight;
          }

          const MoveID lastOutTo = lastMoveOut(e, m.to);
          if (remainingPinsFromBeginningOfMovePhase(e, m.to) == 0
              && firstMoveIn(e, m.to) == globalMoveID && lastOutTo < globalMoveID) {
            gain -= edgeWeight;
          }

        } else {

          const HypernodeID u = m.node;
          const HypernodeID v = phg.graphEdgeHead(e_global, u);
          const MoveID move_id_of_v = tracker.moveOfNode[v];
          const bool v_moved = !tracker.isMoveStale(move_id_of_v) && tracker.isMoveStillValid(move_id_of_v);

          PartitionID pv;
          if (v_moved) {
            const MoveID local_move_id_of_v = move_id_of_v - firstMoveID;
            const Move& move_of_v = move_order[local_move_id_of_v];
            if (local_move_id_of_v < localMoveID) {
              pv = move_of_v.to;
            } else {
              pv = move_of_v.from;
            }
          } else {
            pv = phg.partID(v);
          }

          if (pv == m.to) {
            gain += edgeWeight;
          } else if (pv == m.from) {
            gain -= edgeWeight;
          }

        }

      }

      m.gain = gain;
    });

    timer.stop_timer("gain_recalculation");

    HEAVY_REFINEMENT_ASSERT([&] {
      // verify that all updates were correct (except moveFromBenefit of moved nodes)
      for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
        phg.recomputeMoveFromBenefit(move_order[localMoveID].node);
      }
      phg.checkTrackedPartitionInformation();

      // revert all moves
      for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
        const Move& m = sharedData.moveTracker.moveOrder[localMoveID];
        if (sharedData.moveTracker.isMoveStillValid(m))
          phg.changeNodePartFullUpdate(m.node, m.to, m.from, std::numeric_limits<HypernodeWeight>::max(), []{});
      }

      for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
        if (sharedData.moveTracker.isMoveStillValid(move_order[localMoveID]))
          phg.recomputeMoveFromBenefit(move_order[localMoveID].node);
      }

      // roll forward sequentially and check gains
      for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
        const Move& m = sharedData.moveTracker.moveOrder[localMoveID];
        if (!sharedData.moveTracker.isMoveStillValid(m))
          continue;

        const Gain estimated_gain = phg.km1Gain(m.node, m.from, m.to);
        assert(phg.moveFromBenefit(m.node) == phg.moveFromBenefitRecomputed(m.node));
        assert(phg.moveToPenalty(m.node, m.to) == phg.moveToPenaltyRecomputed(m.node, m.to));
        const HyperedgeWeight km1_before_move = metrics::km1(phg, false);
        phg.changeNodePartFullUpdate(m.node, m.from, m.to, std::numeric_limits<HypernodeWeight>::max(), []{});
        const HyperedgeWeight km1_after_move = metrics::km1(phg, false);
        assert(km1_after_move + estimated_gain == km1_before_move);
        assert(km1_after_move + m.gain == km1_before_move);
        assert(estimated_gain == m.gain);
      }

      for (MoveID localMoveID = 0; localMoveID < sharedData.moveTracker.numPerformedMoves(); ++localMoveID) {
        phg.recomputeMoveFromBenefit(move_order[localMoveID].node);
      }
      return true;
    }());
  }

  MoveID lastMoveOut(HyperedgeID he, PartitionID block) const {
    return last_move_out[he * numParts + block].load(std::memory_order_relaxed);
  }

  MoveID firstMoveIn(HyperedgeID he, PartitionID block) const {
    return first_move_in[he * numParts + block].load(std::memory_order_relaxed);
  }


  HypernodeID remainingPinsFromBeginningOfMovePhase(HyperedgeID he, PartitionID block) const {
    return remaining_original_pins[he * numParts + block].load(std::memory_order_relaxed);
  }

  void resetStoredMoveIDs() {
    tbb::parallel_invoke(
            [&] {
              tbb::parallel_for_each(last_move_out, [](auto& x) { x.store(0, std::memory_order_relaxed); });
            },
            [&] {
              tbb::parallel_for_each(first_move_in, [](auto& x) { x.store(0, std::memory_order_relaxed); });
            }
    );
  }

  void setRemainingOriginalPins(PartitionedHypergraph& phg) {
    if (phg.hypergraph().maxEdgeSize() <= 2) {
      return;
    }

    phg.doParallelForAllEdges([&](const HyperedgeID& he) {
      if (phg.edgeSize(he) > 2) {
        const HyperedgeID he_local = phg.nonGraphEdgeID(he);
        for ( PartitionID block = 0; block < numParts; ++block ) {
          remaining_original_pins[he_local * numParts + block].store(phg.pinCountInPart(he, block), std::memory_order_relaxed);
        }
      }
    });
  }

  PartitionID numParts;

  // ! For each hyperedge and block, the number of pins_in_part
  // ! at the beginning of a move phase minus the number of moved out pins
  // ! A value of zero means at some point a gain improvement might have occured, which then has to be checked
  // ! with the move ID flagging
  vec<CAtomic<HypernodeID>> remaining_original_pins;
  // ! For each hyperedge and each block, the ID of the first move to place a pin in that block
  vec<CAtomic<MoveID>> first_move_in;
  // ! For each hyperedge and each block, the ID of the last move to remove a pin from that block
  vec<CAtomic<MoveID>> last_move_out;
};

}