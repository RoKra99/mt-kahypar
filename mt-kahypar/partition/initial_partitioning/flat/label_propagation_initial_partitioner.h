/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "tbb/task.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/initial_partitioning/flat/initial_partitioning_data_container.h"
#include "mt-kahypar/partition/initial_partitioning/flat/policies/gain_computation_policy.h"

namespace mt_kahypar {

class LabelPropagationInitialPartitioner : public tbb::task {
  using HyperGraph = PartitionedHypergraph<false>;

  using DeltaFunction = std::function<void (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID)>;
  #define NOOP_FUNC [] (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;
  static PartitionID kInvalidPartition;
  static HypernodeID kInvalidHypernode;

  struct MaxGainMove {
    const PartitionID block;
    const Gain gain;
  };

 public:
  LabelPropagationInitialPartitioner(const InitialPartitioningAlgorithm,
                                      InitialPartitioningDataContainer& ip_data,
                                      const Context& context) :
    _ip_data(ip_data),
    _context(context),
    _valid_blocks(context.partition.k),
    _tmp_scores(context.partition.k) { }

  tbb::task* execute() override {
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    HyperGraph& hg = _ip_data.local_partitioned_hypergraph();
    _ip_data.reset_unassigned_hypernodes();

    parallel::scalable_vector<HypernodeID> start_nodes =
      PseudoPeripheralStartNodes::computeStartNodes(_ip_data, _context);
    for ( PartitionID block = 0; block < _context.partition.k; ++block ) {
      if ( hg.partID(start_nodes[block]) == kInvalidPartition ) {
        hg.setNodePart(start_nodes[block], block);
      }
    }
    // Each block is extended with 5 additional vertices which are adjacent
    // to their corresponding seed vertices. This should prevent that block
    // becomes empty after several label propagation rounds.
    for ( PartitionID block = 0; block < _context.partition.k; ++block ) {
      if ( hg.partID(start_nodes[block]) == block ) {
        extendBlockToInitialBlockSize(hg, start_nodes[block], block);
      }
    }

    bool converged = false;
    for ( size_t i = 0; i < _context.initial_partitioning.lp_maximum_iterations && !converged; ++i ) {
      converged = true;

      for ( const HypernodeID& hn : hg.nodes() ) {

        if (hg.nodeDegree(hn) > 0) {
          // Assign vertex to the block where FM gain is maximized
          MaxGainMove max_gain_move = computeMaxGainMove(hg, hn);

          const PartitionID to = max_gain_move.block;
          if ( to != kInvalidPartition ) {
            const PartitionID from = hg.partID(hn);
            if ( from == kInvalidPartition ) {
              ASSERT(fitsIntoBlock(hg, hn, to));

              HEAVY_INITIAL_PARTITIONING_ASSERT([&] {
                Gain expected_gain = CutGainPolicy::calculateGain(hg, hn, to);
                if ( expected_gain != max_gain_move.gain ) {
                  LOG << V(hn);
                  LOG << V(from);
                  LOG << V(to);
                  LOG << V(max_gain_move.gain);
                  LOG << V(expected_gain);
                }
                return true;
              }(), "Gain calculation failed");

              hg.setNodePart(hn, to);
            } else if ( from != to ) {
              ASSERT(fitsIntoBlock(hg, hn, to));

              #ifndef KAHYPAR_ENABLE_HEAVY_INITIAL_PARTITIONING_ASSERTIONS
              hg.changeNodePart(hn, from, to, NOOP_FUNC);
              #else
              Gain expected_gain = 0;
              auto cut_delta = [&](const HyperedgeID he,
                                   const HyperedgeWeight edge_weight,
                                   const HypernodeID,
                                   const HypernodeID pin_count_in_from_part_after,
                                   const HypernodeID pin_count_in_to_part_after) {
                HypernodeID adjusted_edge_size = 0;
                for ( const HypernodeID& pin : hg.pins(he) ) {
                  if ( hg.partID(pin) != kInvalidPartition ) {
                    ++adjusted_edge_size;
                  }
                }
                expected_gain -= HyperGraph::cutDelta(he, edge_weight, adjusted_edge_size,
                  pin_count_in_from_part_after, pin_count_in_to_part_after);
              };
              hg.changeNodePart(hn, from, to, cut_delta);
              ASSERT(expected_gain == max_gain_move.gain, "Gain calculation failed"
                << V(expected_gain) << V(max_gain_move.gain));
              #endif
            }
          }

        } else if ( hg.partID(hn) == kInvalidPartition ) {
          // In case vertex hn is a degree zero vertex we assign it
          // to the block with minimum weight
          assignVertexToBlockWithMinimumWeight(hg, hn);
        }

      }

    }

    // If there are still unassigned vertices left, we assign them to the
    // block with minimum weight.
    while ( _ip_data.get_unassigned_hypernode() != kInvalidHypernode ) {
      const HypernodeID unassigned_hn = _ip_data.get_unassigned_hypernode();
      assignVertexToBlockWithMinimumWeight(hg, unassigned_hn);
    }

    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(end - start).count();
    _ip_data.commit(InitialPartitioningAlgorithm::label_propagation, time);
    return nullptr;
  }

 private:
  bool fitsIntoBlock(HyperGraph& hypergraph,
                     const HypernodeID hn,
                     const PartitionID block) const {
    ASSERT(block != kInvalidPartition && block < _context.partition.k);
    return hypergraph.partWeight(block) + hypergraph.nodeWeight(hn) <=
      _context.partition.max_part_weights[block];
  }

  MaxGainMove computeMaxGainMove(HyperGraph& hypergraph,
                                 const HypernodeID hn) {
    if ( hypergraph.partID(hn) == kInvalidPartition ) {
      return computeMaxGainMoveForUnassignedVertex(hypergraph, hn);
    } else {
      return computeMaxGainMoveForAssignedVertex(hypergraph, hn);
    }
  }

  MaxGainMove computeMaxGainMoveForUnassignedVertex(HyperGraph& hypergraph,
                                                    const HypernodeID hn) {
    ASSERT(hypergraph.partID(hn) == kInvalidPartition);
    ASSERT(std::all_of(_tmp_scores.begin(), _tmp_scores.end(), [](Gain i) { return i == 0; }),
           "Temp gain array not initialized properly");
    _valid_blocks.reset();

    HyperedgeWeight internal_weight = 0;
    for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
      const HyperedgeWeight he_weight = hypergraph.edgeWeight(he);
      if (hypergraph.connectivity(he) == 1) {
        // In case, connectivity is one we would make the hyperedge cut if would
        // assign the vertex to an different block than the one already contained
        // in the hyperedge
        const PartitionID connected_block = *hypergraph.connectivitySet(he).begin();
        _valid_blocks.set(connected_block, true);
        internal_weight += he_weight;
        _tmp_scores[connected_block] += he_weight;
      } else {
        // Otherwise we can assign the vertex to a block already contained
        // in the hyperedge without affecting cut
        for (const PartitionID& target_part : hypergraph.connectivitySet(he)) {
          _valid_blocks.set(target_part, true);
        }
      }
    }

    return findMaxGainMove(hypergraph, hn, internal_weight);
  }

  MaxGainMove computeMaxGainMoveForAssignedVertex(HyperGraph& hypergraph,
                                                  const HypernodeID hn) {
    ASSERT(hypergraph.partID(hn) != kInvalidPartition);
    ASSERT(std::all_of(_tmp_scores.begin(), _tmp_scores.end(), [](Gain i) { return i == 0; }),
           "Temp gain array not initialized properly");
    _valid_blocks.reset();

    const PartitionID from = hypergraph.partID(hn);
    HyperedgeWeight internal_weight = 0;
    for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
      const HyperedgeWeight he_weight = hypergraph.edgeWeight(he);
      const PartitionID connectivity = hypergraph.connectivity(he);
      const HypernodeID pins_in_from_part = hypergraph.pinCountInPart(he, from);

      if ( connectivity == 1 && pins_in_from_part > 1 ) {
        // If connectivity is one and there is more than one vertex in block
        // of hypernode hn, we would make the hyperedge cut, if we would assign
        // hn to an different block.
        internal_weight += he_weight;
      } else if ( connectivity == 2 ) {
        for (const PartitionID& to : hypergraph.connectivitySet(he)) {
          _valid_blocks.set(to, true);
          // In case connectivity is two and hn is the last vertex in hyperedge
          // he of block from, we would make that hyperedge a non-cut hyperedge.
          if ( pins_in_from_part == 1 && hypergraph.pinCountInPart(he, to) > 0 ) {
            _tmp_scores[to] += he_weight;
          }
        }
      } else {
        // Otherwise we can assign the vertex to a block already contained
        // in the hyperedge without affecting cut
        for (const PartitionID& to : hypergraph.connectivitySet(he)) {
          _valid_blocks.set(to, true);
        }
      }
    }

    return findMaxGainMove(hypergraph, hn, internal_weight);
  }

  MaxGainMove findMaxGainMove(HyperGraph& hypergraph,
                              const HypernodeID hn,
                              const HypernodeWeight internal_weight) {
    const PartitionID from = hypergraph.partID(hn);
    PartitionID best_block = from;
    Gain best_score = from == kInvalidPartition ? std::numeric_limits<Gain>::min() : 0;
    for (PartitionID block = 0; block < _context.partition.k; ++block) {
      if (from != block && _valid_blocks[block]) {
        _tmp_scores[block] -= internal_weight;

        // Since we perform size-constraint label propagation, the move to the
        // corresponding block is only valid, if it fullfils the balanced constraint.
        if (fitsIntoBlock(hypergraph, hn, block) && _tmp_scores[block] > best_score) {
          best_score = _tmp_scores[block];
          best_block = block;
        }
      }
      _tmp_scores[block] = 0;
    }
    return MaxGainMove { best_block, best_score };
  }

  void extendBlockToInitialBlockSize(HyperGraph& hypergraph,
                                     const HypernodeID seed_vertex,
                                     const PartitionID block) {
    ASSERT(hypergraph.partID(seed_vertex) == block);
    // We search for _context.initial_partitioning.lp_initial_block_size vertices
    // around the seed vertex to extend the corresponding block
    for ( const HyperedgeID& he : hypergraph.incidentEdges(seed_vertex) ) {
      for ( const HypernodeID& pin : hypergraph.pins(he) ) {
        if ( hypergraph.partID(pin) == kInvalidPartition ) {
          hypergraph.setNodePart(pin, block);
          if ( hypergraph.partSize(block) ==
              _context.initial_partitioning.lp_initial_block_size ) {
            break;
          }
        }
      }
      if ( hypergraph.partSize(block) ==
           _context.initial_partitioning.lp_initial_block_size ) {
        break;
      }
    }

    // If there are last than _context.initial_partitioning.lp_initial_block_size
    // adjacent vertices to the seed vertex, we find a new seed vertex and call
    // this function recursive
    const HypernodeID part_size = hypergraph.partSize(block);
    if ( part_size < _context.initial_partitioning.lp_initial_block_size ) {
      const HypernodeID new_seed_vertex = _ip_data.get_unassigned_hypernode();
      if ( new_seed_vertex != kInvalidHypernode  ) {
        hypergraph.setNodePart(new_seed_vertex, block);
        if ( part_size + 1 < _context.initial_partitioning.lp_initial_block_size ) {
          extendBlockToInitialBlockSize(hypergraph, new_seed_vertex, block);
        }
      }
    }
  }

  void assignVertexToBlockWithMinimumWeight(HyperGraph& hypergraph, const HypernodeID hn) {
    ASSERT(hypergraph.partID(hn) == kInvalidPartition);
    PartitionID minimum_weight_block = kInvalidPartition;
    HypernodeWeight minimum_weight = std::numeric_limits<HypernodeWeight>::max();
    for ( PartitionID block = 0; block < _context.partition.k; ++block ) {
      const HypernodeWeight block_weight = hypergraph.partWeight(block);
      if ( block_weight < minimum_weight ) {
        minimum_weight = block_weight;
        minimum_weight_block = block;
      }
    }
    ASSERT(minimum_weight_block != kInvalidPartition);
    hypergraph.setNodePart(hn, minimum_weight_block);
  }

  InitialPartitioningDataContainer& _ip_data;
  const Context& _context;
  kahypar::ds::FastResetFlagArray<> _valid_blocks;
  parallel::scalable_vector<Gain> _tmp_scores;
};

PartitionID LabelPropagationInitialPartitioner::kInvalidPartition = -1;
HypernodeID LabelPropagationInitialPartitioner::kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

} // namespace mt_kahypar
