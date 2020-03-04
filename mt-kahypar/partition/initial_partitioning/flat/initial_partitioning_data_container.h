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

#include <sstream>

#include "tbb/enumerable_thread_specific.h"

#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/i_refiner.h"
#include "mt-kahypar/partition/refinement/label_propagation_refiner.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/utils/initial_partitioning_stats.h"

namespace mt_kahypar {

template <typename TypeTraits>
class InitialPartitioningDataContainerT {
  using HyperGraph = typename TypeTraits::HyperGraph;
  using PartitionedHyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
  using PartitionedHyperGraphWithoutBorderVertices = typename TypeTraits::template PartitionedHyperGraph<false>;
  using TBB = typename TypeTraits::TBB;
  using Refiner = IRefinerT<TypeTraits, false>;

  static constexpr bool debug = false;
  static PartitionID kInvalidPartition;
  static HypernodeID kInvalidHypernode;

  struct PartitioningResult {
    explicit PartitioningResult(InitialPartitioningAlgorithm algorithm,
                                HyperedgeWeight objective,
                                double imbalance) :
      _algorithm(algorithm),
      _objective(objective),
      _imbalance(imbalance) { }

    bool is_other_better(const PartitioningResult& other, const double epsilon) {
      bool equal_metric = other._objective == _objective;
      bool improved_metric = other._objective < _objective;
      bool improved_imbalance = other._imbalance < _imbalance;
      bool is_feasible = _imbalance <= epsilon;
      bool is_other_feasible = other._imbalance <= epsilon;
      return ( improved_metric && (is_other_feasible || improved_imbalance) ) ||
             ( equal_metric && improved_imbalance ) ||
             ( is_other_feasible && !is_feasible ) ||
             ( improved_imbalance && !is_other_feasible && !is_feasible );
    }

    std::string str() const {
      std::stringstream ss;
      ss << "Algorithm = " << _algorithm << ", "
         << "Objective = " << _objective << ", "
         << "Imbalance = " << _imbalance;
      return ss.str();
    }

    InitialPartitioningAlgorithm _algorithm;
    HyperedgeWeight _objective;
    double _imbalance;
  };

  struct LocalInitialPartitioningHypergraph {
    using LabelPropagationKm1Refiner = LabelPropagationRefinerT<TypeTraits, Km1Policy, false>;
    using LabelPropagationCutRefiner = LabelPropagationRefinerT<TypeTraits, CutPolicy, false>;

    LocalInitialPartitioningHypergraph(HyperGraph& hypergraph,
                                       const Context& context,
                                       const TaskGroupID task_group_id) :
      _partitioned_hypergraph(context.partition.k, hypergraph),
      _context(context),
      _partition(hypergraph.initialNumNodes(), kInvalidPartition),
      _result(InitialPartitioningAlgorithm::UNDEFINED,
              std::numeric_limits<HypernodeWeight>::max(),
              std::numeric_limits<double>::max()),
      _label_propagation(nullptr),
      _stats() {

      for ( uint8_t algo = 0; algo < static_cast<size_t>(InitialPartitioningAlgorithm::UNDEFINED); ++algo ) {
        _stats.emplace_back(static_cast<InitialPartitioningAlgorithm>(algo));
      }

      if ( _context.refinement.label_propagation.algorithm != LabelPropagationAlgorithm::do_nothing ) {
        if ( _context.partition.objective == kahypar::Objective::km1 ) {
          _label_propagation = std::make_unique<LabelPropagationKm1Refiner>(
            _partitioned_hypergraph, _context, task_group_id);
        } else if ( _context.partition.objective == kahypar::Objective::cut ) {
          _label_propagation = std::make_unique<LabelPropagationCutRefiner>(
            _partitioned_hypergraph, _context, task_group_id);
        }
      }
    }

    void commit(const InitialPartitioningAlgorithm algorithm, const double time = 0.0) {
      ASSERT([&]() {
          for (const HypernodeID& hn : _partitioned_hypergraph.nodes()) {
            if (_partitioned_hypergraph.partID(hn) == kInvalidPartition) {
              return false;
            }
          }
          return true;
        } (), "There are unassigned hypernodes!");

      _partitioned_hypergraph.initializeNumCutHyperedges();

      kahypar::Metrics current_metric = {
        metrics::hyperedgeCut(_partitioned_hypergraph),
        metrics::km1(_partitioned_hypergraph),
        metrics::imbalance(_partitioned_hypergraph, _context) };

      if ( _label_propagation ) {
        _label_propagation->initialize(_partitioned_hypergraph);
        _label_propagation->refine(_partitioned_hypergraph, {}, current_metric);
      }

      PartitioningResult result(algorithm,
        current_metric.getMetric(kahypar::Mode::direct_kway, _context.partition.objective),
        current_metric.imbalance);

      DBG << "Thread ID:" << std::this_thread::get_id()
          << "- CPU ID:" << sched_getcpu()
          << "[" << result.str() << "]";

      // Aggregate Stats
      uint8_t algorithm_index = static_cast<uint8_t>(algorithm);
      _stats[algorithm_index].total_sum_quality += result._objective;
      _stats[algorithm_index].total_time += time;
      ++_stats[algorithm_index].total_calls;

      if ( _result.is_other_better(result, _context.partition.epsilon) ) {
        for ( const HypernodeID& hn : _partitioned_hypergraph.nodes() ) {
          const HypernodeID original_id = _partitioned_hypergraph.originalNodeID(hn);
          const PartitionID part_id = _partitioned_hypergraph.partID(hn);
          ASSERT(original_id < _partition.size());
          ASSERT(part_id != kInvalidPartition);
          _partition[original_id] = part_id;
        }
        _result = std::move(result);
      }

      _partitioned_hypergraph.resetPartition();
    }

    void aggregate_stats(parallel::scalable_vector<utils::InitialPartitionerSummary>& main_stats) const {
      ASSERT(main_stats.size() == _stats.size());
      for ( size_t i = 0; i < _stats.size(); ++i ) {
        main_stats[i].add(_stats[i]);
      }
    }

    PartitionedHyperGraphWithoutBorderVertices _partitioned_hypergraph;
    const Context& _context;
    parallel::scalable_vector<PartitionID> _partition;
    PartitioningResult _result;
    std::unique_ptr<Refiner> _label_propagation;
    parallel::scalable_vector<utils::InitialPartitionerSummary> _stats;
  };

  using ThreadLocalHypergraph = tbb::enumerable_thread_specific<LocalInitialPartitioningHypergraph>;
  using ThreadLocalUnassignedHypernodes = tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeID>>;

 public:
  InitialPartitioningDataContainerT(PartitionedHyperGraph& hypergraph,
                                    const Context& context,
                                    const TaskGroupID task_group_id) :
    _partitioned_hg(hypergraph),
    _context(context),
    _task_group_id(task_group_id),
    _local_hg([&] {
      return construct_local_partitioned_hypergraph();
    }),
    _local_kway_pq(_context.partition.k),
    _is_local_pq_initialized(false),
    _local_hn_visited(_context.partition.k * hypergraph.initialNumNodes()),
    _local_he_visited(_context.partition.k * hypergraph.initialNumEdges()),
    _local_unassigned_hypernodes(),
    _local_unassigned_hypernode_pointer(std::numeric_limits<size_t>::max())  {
    // Setup Label Propagation Refiner Config for Initial Partitioning
    _context.refinement = _context.initial_partitioning.refinement;
    _context.refinement.label_propagation.execute_sequential = true;
  }

  InitialPartitioningDataContainerT(const InitialPartitioningDataContainerT&) = delete;
  InitialPartitioningDataContainerT & operator= (const InitialPartitioningDataContainerT &) = delete;

  InitialPartitioningDataContainerT(InitialPartitioningDataContainerT&&) = delete;
  InitialPartitioningDataContainerT & operator= (InitialPartitioningDataContainerT &&) = delete;

  PartitionedHyperGraph& global_partitioned_hypergraph() {
    return _partitioned_hg;
  }

  PartitionedHyperGraphWithoutBorderVertices& local_partitioned_hypergraph() {
    return _local_hg.local()._partitioned_hypergraph;
  }

  KWayPriorityQueue& local_kway_priority_queue() {
    bool& is_local_pq_initialized = _is_local_pq_initialized.local();
    KWayPriorityQueue& local_kway_pq = _local_kway_pq.local();
    if ( !is_local_pq_initialized ) {
      local_kway_pq.initialize(local_partitioned_hypergraph().initialNumNodes());
      is_local_pq_initialized = true;
    }
    return local_kway_pq;
  }

  kahypar::ds::FastResetFlagArray<>& local_hypernode_fast_reset_flag_array() {
    return _local_hn_visited.local();
  }

  kahypar::ds::FastResetFlagArray<>& local_hyperedge_fast_reset_flag_array() {
    return _local_he_visited.local();
  }

  void reset_unassigned_hypernodes() {
    parallel::scalable_vector<HypernodeID>& unassigned_hypernodes =
      _local_unassigned_hypernodes.local();
    size_t& unassigned_hypernode_pointer = _local_unassigned_hypernode_pointer.local();
    if ( unassigned_hypernode_pointer == std::numeric_limits<size_t>::max() ) {
      // In case the local unassigned hypernode vector was not initialized before
      // we initialize it here
      const PartitionedHyperGraphWithoutBorderVertices& hypergraph = local_partitioned_hypergraph();
      for ( const HypernodeID& hn : hypergraph.nodes() ) {
        unassigned_hypernodes.push_back(hn);
      }
    }
    unassigned_hypernode_pointer = unassigned_hypernodes.size();
  }

  HypernodeID get_unassigned_hypernode(const PartitionID unassigned_block = kInvalidPartition) {
    const PartitionedHyperGraphWithoutBorderVertices& hypergraph = local_partitioned_hypergraph();
    parallel::scalable_vector<HypernodeID>& unassigned_hypernodes =
      _local_unassigned_hypernodes.local();
    size_t& unassigned_hypernode_pointer = _local_unassigned_hypernode_pointer.local();
    ASSERT(unassigned_hypernodes.size() > 0);
    ASSERT(unassigned_hypernode_pointer <= unassigned_hypernodes.size());

    while ( unassigned_hypernode_pointer > 0 ) {
      const HypernodeID current_hn = unassigned_hypernodes[0];
      // In case the current hypernode is unassigned we return it
      if ( hypergraph.partID(current_hn) == unassigned_block ) {
        return current_hn;
      }
      // In case the hypernode on the first position is already assigned,
      // we swap it to end of the unassigned hypernode vector and decrement
      // the pointer such that we will not visit it again
      std::swap(unassigned_hypernodes[0], unassigned_hypernodes[--unassigned_hypernode_pointer]);
    }
    return kInvalidHypernode;
  }

  /*!
   * Commits the current partition computed on the local hypergraph. Partition replaces
   * the best local partition, if it has a better quality (or better imbalance).
   * Partition on the local hypergraph is resetted afterwards.
   */
  void commit(const InitialPartitioningAlgorithm algorithm, const double time = 0.0) {
    _local_hg.local().commit(algorithm, time);
  }

  /*!
   * Determines the best partition computed by all threads and applies it to
   * the hypergraph. Note, this function is not thread-safe and should be called
   * if no other thread using that object operates on it.
   */
  void apply() {
    // Initialize Stats
    parallel::scalable_vector<utils::InitialPartitionerSummary> stats;
    size_t number_of_threads = 0;
    for ( uint8_t algo = 0; algo < static_cast<size_t>(InitialPartitioningAlgorithm::UNDEFINED); ++algo ) {
      stats.emplace_back(static_cast<InitialPartitioningAlgorithm>(algo));
    }

    // Determine best partition
    LocalInitialPartitioningHypergraph* best = nullptr;
    LocalInitialPartitioningHypergraph* worst = nullptr;
    LocalInitialPartitioningHypergraph* best_imbalance = nullptr;
    LocalInitialPartitioningHypergraph* best_objective = nullptr;
    for ( LocalInitialPartitioningHypergraph& partition : _local_hg ) {
      ++number_of_threads;
      partition.aggregate_stats(stats);
      if ( !best || best->_result.is_other_better(partition._result, _context.partition.epsilon) ) {
        best = &partition;
      }
      if ( !worst || !worst->_result.is_other_better(partition._result, _context.partition.epsilon) ) {
        worst = &partition;
      }
      if ( !best_imbalance || best_imbalance->_result._imbalance > partition._result._imbalance ||
           (best_imbalance->_result._imbalance == partition._result._imbalance &&
            best_objective->_result._objective > partition._result._objective)) {
        best_imbalance = &partition;
      }
      if ( !best_objective || best_objective->_result._objective > partition._result._objective ) {
        best_objective = &partition;
      }
    }

    ASSERT(best);
    ASSERT(worst);
    ASSERT(best_imbalance);
    ASSERT(best_objective);
    DBG << "Num Vertices =" << _partitioned_hg.initialNumNodes()
        << ", Num Edges =" << _partitioned_hg.initialNumEdges()
        << ", k =" << _context.partition.k << ", epsilon =" << _context.partition.epsilon;
    DBG << "Best Partition                [" << best->_result.str() << "]";
    DBG << "Worst Partition               [" << worst->_result.str() << "]";
    DBG << "Best Balanced Partition       [" << best_imbalance->_result.str() << "]";
    DBG << "Partition with Best Objective [" << best_objective->_result.str() << "]";

    // Applies best partition to hypergraph
    _partitioned_hg.doParallelForAllNodes(_task_group_id, [&](const HypernodeID hn) {
      const HypernodeID original_id = _partitioned_hg.originalNodeID(hn);
      ASSERT(original_id < best->_partition.size());
      const PartitionID part_id = best->_partition[original_id];
      ASSERT(part_id != kInvalidPartition && part_id < _partitioned_hg.k());
      _partitioned_hg.setOnlyNodePart(hn, part_id);
    });
    _partitioned_hg.initializePartition(_task_group_id);

    utils::InitialPartitioningStats::instance().add_initial_partitioning_result(
      best->_result._algorithm, number_of_threads, stats);
    ASSERT(best->_result._objective == metrics::objective(_partitioned_hg, _context.partition.objective),
           V(best->_result._objective) << V(metrics::objective(_partitioned_hg, _context.partition.objective)));
  }

 private:
  LocalInitialPartitioningHypergraph construct_local_partitioned_hypergraph() {
    return LocalInitialPartitioningHypergraph(_partitioned_hg.hypergraph(), _context, _task_group_id);
  }

  PartitionedHyperGraph& _partitioned_hg;
  Context _context;
  const TaskGroupID _task_group_id;

  ThreadLocalHypergraph _local_hg;
  ThreadLocalKWayPriorityQueue _local_kway_pq;
  tbb::enumerable_thread_specific<bool> _is_local_pq_initialized;
  ThreadLocalFastResetFlagArray _local_hn_visited;
  ThreadLocalFastResetFlagArray _local_he_visited;
  ThreadLocalUnassignedHypernodes _local_unassigned_hypernodes;
  tbb::enumerable_thread_specific<size_t> _local_unassigned_hypernode_pointer;
};

template <typename TypeTraits>
PartitionID InitialPartitioningDataContainerT<TypeTraits>::kInvalidPartition = -1;
template <typename TypeTraits>
HypernodeID InitialPartitioningDataContainerT<TypeTraits>::kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

using InitialPartitioningDataContainer = InitialPartitioningDataContainerT<GlobalTypeTraits>;
} // namespace mt_kahypar