/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
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
#include <type_traits>

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

#include "robin_hood.h"

namespace mt_kahypar {
namespace ds {

/*!
 * Special hypergraph data structure for the FM algorithm.
 * In order to perform frequent moves on the hypergraph without affecting other
 * local searches one can wrap this hypergraph around a partitioned hypergraph and
 * perform those moves locally on this hypergraph. The changes are not reflected
 * in the global partitioned hypergraph. Instead, the data structure stores deltas
 * relative to the global partitioned hypergraph for each internal member such that
 * by applying those deltas to the original partitioned hypergraph it would reflect the
 * global state of the partitioned hypergraph after all local moves would be applied
 * to it.
 * The rationale behind this is that the majority of local searches to not yield to
 * an improvement and are immediatly reverted. However, applying them directly to global
 * partitioned hypergraph would affect other local searches running concurrently, which build
 * upon that state. This special partitioned hypergraph allows a local search to hide its
 * current search state from other searches in a space efficient manner.
 */
template <typename PartitionedHypergraph = Mandatory>
class DeltaPartitionedHypergraph {
 private:

  using HypernodeIterator = typename PartitionedHypergraph::HypernodeIterator;
  using HyperedgeIterator = typename PartitionedHypergraph::HyperedgeIterator;
  using IncidenceIterator = typename PartitionedHypergraph::IncidenceIterator;
  using IncidentNetsIterator = typename PartitionedHypergraph::IncidentNetsIterator;

 public:
  static constexpr bool supports_connectivity_set = false;
  static constexpr HyperedgeID HIGH_DEGREE_THRESHOLD = PartitionedHypergraph::HIGH_DEGREE_THRESHOLD;

  DeltaPartitionedHypergraph(const PartitionID k) :
    _k(k),
    _phg(nullptr),
    _part_weights_delta(k, 0),
    _part_ids_delta(),
    _pins_in_part_delta(),
    _move_to_penalty_delta(),
    _move_from_benefit_delta() { }

  DeltaPartitionedHypergraph(const DeltaPartitionedHypergraph&) = delete;
  DeltaPartitionedHypergraph & operator= (const DeltaPartitionedHypergraph &) = delete;

  DeltaPartitionedHypergraph(DeltaPartitionedHypergraph&& other) = default;
  DeltaPartitionedHypergraph & operator= (DeltaPartitionedHypergraph&& other) = default;

  ~DeltaPartitionedHypergraph() = default;

  void setPartitionedHypergraph(PartitionedHypergraph* phg) {
    _phg = phg;
  }

  // ####################### Iterators #######################

  // ! Returns an iterator over the set of active nodes of the hypergraph
  IteratorRange<HypernodeIterator> nodes() const {
    ASSERT(_phg);
    return _phg->nodes();
  }

  // ! Returns an iterator over the set of active edges of the hypergraph
  IteratorRange<HyperedgeIterator> edges() const {
    ASSERT(_phg);
    return _phg->edges();
  }

  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetsIterator> incidentEdges(const HypernodeID u) const {
    ASSERT(_phg);
    return _phg->incidentEdges(u);
  }

  // ! Returns a range to loop over the pins of hyperedge e.
  IteratorRange<IncidenceIterator> pins(const HyperedgeID e) const {
    ASSERT(_phg);
    return _phg->pins(e);
  }

  // ####################### Hypernode Information #######################

  HypernodeWeight nodeWeight(const HypernodeID u) const {
    ASSERT(_phg);
    return _phg->nodeWeight(u);
  }

  HyperedgeID nodeDegree(const HypernodeID u) const {
    ASSERT(_phg);
    return _phg->nodeDegree(u);
  }

  // ####################### Hyperedge Information #######################

  // ! Number of pins of a hyperedge
  HypernodeID edgeSize(const HyperedgeID e) const {
    ASSERT(_phg);
    return _phg->edgeSize(e);
  }

  HyperedgeWeight edgeWeight(const HyperedgeID e) const {
    ASSERT(_phg);
    return _phg->edgeWeight(e);
  }

  // ####################### Partition Information #######################

  // ! Changes the block of hypernode u from 'from' to 'to'.
  // ! Move is successful, if it is not violating the balance
  // ! constraint specified by 'max_weight_to'.
  template<typename DeltaFunc>
  bool changeNodePart(const HypernodeID u,
                      const PartitionID from,
                      const PartitionID to,
                      const HypernodeWeight max_weight_to,
                      DeltaFunc&& delta_func) {
    ASSERT(_phg);
    assert(partID(u) == from);
    assert(from != to);
    const HypernodeWeight wu = _phg->nodeWeight(u);
    if ( partWeight(to) + wu <= max_weight_to ) {
      _part_ids_delta[u] = to;
      _part_weights_delta[to] += wu;
      _part_weights_delta[from] -= wu;
      for ( const HyperedgeID& he : _phg->incidentEdges(u) ) {
        const HypernodeID pin_count_in_from_part_after = decrementPinCountInPart(he, from);
        const HypernodeID pin_count_in_to_part_after = incrementPinCountInPart(he, to);
        delta_func(he, _phg->edgeWeight(he), _phg->edgeSize(he), pin_count_in_from_part_after, pin_count_in_to_part_after);
      }
      return true;
    } else {
      return false;
    }
  }

  // curry
  bool changeNodePart(const HypernodeID u,
                      const PartitionID from,
                      const PartitionID to,
                      const HypernodeWeight max_weight_to) {
    return changeNodePart(u, from, to, max_weight_to, NoOpDeltaFunc());
  }

  bool changeNodePartWithGainCacheUpdate(const HypernodeID u,
                                         const PartitionID from,
                                         const PartitionID to,
                                         const HypernodeWeight max_weight_to) {
    auto delta_gain_func = [&]( HyperedgeID he, HyperedgeWeight edge_weight,
                                HypernodeID ,HypernodeID pcip_from, HypernodeID pcip_to ) {
      gainCacheUpdate(he, edge_weight, from, pcip_from, to, pcip_to);
    };
    return changeNodePart(u, from, to, max_weight_to, delta_gain_func);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void gainCacheUpdate(const HyperedgeID he, const HyperedgeWeight we,
                       const PartitionID from, const HypernodeID pin_count_in_from_part_after,
                       const PartitionID to, const HypernodeID pin_count_in_to_part_after) {

    if (pin_count_in_from_part_after == 1) {
      for (HypernodeID u : pins(he)) {
        if (partID(u) == from) {
          _move_from_benefit_delta[u] += we;
        }
      }
    } else if (pin_count_in_from_part_after == 0) {
      for (HypernodeID u : pins(he)) {
        _move_to_penalty_delta[penalty_index(u, from)] += we;
      }
    }

    if (pin_count_in_to_part_after == 1) {
      for (HypernodeID u : pins(he)) {
        _move_to_penalty_delta[penalty_index(u, to)] -= we;
      }
    } else if (pin_count_in_to_part_after == 2) {
      for (HypernodeID u : pins(he)) {
        if (partID(u) == to) {
          _move_from_benefit_delta[u] -= we;
        }
      }
    }
  }


  // ! Returns the block of hypernode u
  PartitionID partID(const HypernodeID u) const {
    ASSERT(_phg);
    const auto it = _part_ids_delta.find(u);
    return it != _part_ids_delta.cend() ? it->second : _phg->partID(u);
  }

  // ! Returns the total weight of block p
  HypernodeWeight partWeight(const PartitionID p) const {
    ASSERT(_phg);
    ASSERT(p != kInvalidPartition && p < _k);
    return _phg->partWeight(p) + _part_weights_delta[p];
  }

  // ! Returns the number of pins of hyperedge e in block p
  HypernodeID pinCountInPart(const HyperedgeID e, const PartitionID p) const {
    ASSERT(_phg);
    ASSERT(p != kInvalidPartition && p < _k);
    const auto it = _pins_in_part_delta.find(e * _k + p);
    return std::max(static_cast<int32_t>(_phg->pinCountInPart(e, p)) +
      ( it != _pins_in_part_delta.cend() ? it->second : 0 ), 0);
  }

  // ! Returns the sum of all edges incident to u, where u is the last remaining
  // ! pin in its block
  HyperedgeWeight moveFromBenefit(const HypernodeID u) const {
    ASSERT(_phg);
    const auto it = _move_from_benefit_delta.find(u);
    return _phg->moveFromBenefit(u) + ( it != _move_from_benefit_delta.cend() ? it->second : 0 );
  }

  // ! Returns the sum of all edges incident to u, where p is not part of
  // ! their connectivity set.
  HyperedgeWeight moveToPenalty(const HypernodeID u, const PartitionID p) const {
    ASSERT(_phg);
    ASSERT(p != kInvalidPartition && p < _k);
    const auto it = _move_to_penalty_delta.find(penalty_index(u, p));
    return _phg->moveToPenalty(u, p) + ( it != _move_to_penalty_delta.cend() ? it->second : 0 );
  }

  Gain km1Gain(const HypernodeID u, const PartitionID from, const PartitionID to) const {
    unused(from);
    ASSERT(from == partID(u), "While gain computation works for from != partID(u), such a query makes no sense");
    ASSERT(from != to, "The gain computation doesn't work for from = to");
    return moveFromBenefit(u) - moveToPenalty(u, to);
  }

  void initializeGainCacheEntry(const HypernodeID u, vec<Gain>& penalty_aggregator) {
    _phg->initializeGainCacheEntry(u, penalty_aggregator);
  }

  // ! Clears all deltas applied to the partitioned hypergraph
  void clear() {
    // O(k)
    _part_weights_delta.assign(_k, 0);
    // Constant Time
    _part_ids_delta.clear();
    _pins_in_part_delta.clear();
    _move_to_penalty_delta.clear();
    _move_from_benefit_delta.clear();
  }

  void dropMemory() {
    if (!_memory_dropped) {
      _memory_dropped = true;
      _part_ids_delta.clear(); _part_ids_delta.rehash(0);
      _pins_in_part_delta.clear(); _pins_in_part_delta.rehash(0);
      _move_to_penalty_delta.clear(); _move_to_penalty_delta.rehash(0);
      _move_from_benefit_delta.clear(); _move_from_benefit_delta.rehash(0);
    }
  }

  size_t combinedMemoryConsumption() const {
    return _part_ids_delta.calcNumBytesTotal(_part_ids_delta.size())
            + _pins_in_part_delta.calcNumBytesTotal(_pins_in_part_delta.size())
            + _move_to_penalty_delta.calcNumBytesTotal(_move_to_penalty_delta.size())
            + _move_from_benefit_delta.calcNumBytesTotal(_move_from_benefit_delta.size());
  }

  PartitionID k() const {
    return _k;
  }

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode* delta_phg_node = parent->addChild("Delta Partitioned Hypergraph");
    utils::MemoryTreeNode* part_weights_node = delta_phg_node->addChild("Delta Part Weights");
    part_weights_node->updateSize(_part_weights_delta.capacity() * sizeof(HypernodeWeight));
    utils::MemoryTreeNode* part_ids_node = delta_phg_node->addChild("Delta Part IDs");
    part_ids_node->updateSize(_part_ids_delta.calcNumBytesTotal(_part_ids_delta.size()));
    utils::MemoryTreeNode* pins_in_part_node = delta_phg_node->addChild("Delta Pins In Part");
    pins_in_part_node->updateSize(_pins_in_part_delta.calcNumBytesTotal(_pins_in_part_delta.size()));
    utils::MemoryTreeNode* move_from_benefit_node = delta_phg_node->addChild("Delta Move From Benefit");
    move_from_benefit_node->updateSize(_move_from_benefit_delta.calcNumBytesTotal(_move_from_benefit_delta.size()));
    utils::MemoryTreeNode* move_to_penalty_node = delta_phg_node->addChild("Delta Move To Penalty");
    move_to_penalty_node->updateSize(_move_to_penalty_delta.calcNumBytesTotal(_move_to_penalty_delta.size()));
  }

 private:

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  size_t penalty_index(const HypernodeID u, const PartitionID p) const {
    return size_t(u) * _k + p;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID decrementPinCountInPart(const HyperedgeID e, const PartitionID p) {
    return std::max(static_cast<int32_t>(
      _phg->pinCountInPart(e, p)) + --_pins_in_part_delta[e * _k + p], static_cast<int32_t>(0));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID incrementPinCountInPart(const HyperedgeID e, const PartitionID p) {
    return std::max(static_cast<int32_t>(
      _phg->pinCountInPart(e, p)) + ++_pins_in_part_delta[e * _k + p], static_cast<int32_t>(0));
  }

  bool _memory_dropped = false;

  // ! Number of blocks
  const PartitionID _k;

  // ! Partitioned hypergraph where all deltas are stored relative to
  PartitionedHypergraph* _phg;

  // ! Delta for block weights
  vec< HypernodeWeight > _part_weights_delta;

  // ! Stores for each locally moved node, its new block id
  robin_hood::unordered_map<HypernodeID, PartitionID> _part_ids_delta;

  // ! Stores the delta of each locally touched pin count entry
  // ! relative to the _pins_in_part member in '_phg'
  robin_hood::unordered_map<size_t, int32_t> _pins_in_part_delta;

  // ! Stores the delta of each locally touched move to penalty entry
  // ! relative to the _move_to_penalty member in '_phg'
  robin_hood::unordered_map<size_t, HyperedgeWeight> _move_to_penalty_delta;

  // ! Stores the delta of each locally touched move from benefit entry
  // ! relative to the _move_from_benefit member in '_phg'
  robin_hood::unordered_map<HypernodeID, HyperedgeWeight> _move_from_benefit_delta;
};

} // namespace ds
} // namespace mt_kahypar