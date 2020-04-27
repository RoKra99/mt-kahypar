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

#include <algorithm>
#include <limits>
#include <vector>

#include "tbb/parallel_invoke.h"
#include "tbb/task_arena.h"
#include "tbb/task_group.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

/*!
 * RECURSIVE BISECTION INITIAL PARTITIONER
 * The recursive bisection initial partitioner starts by performing a parallel multilevel bisection.
 * Once the hypergraph is bisected both blocks are partitioned recursively in parallel until the
 * desired number of blocks are reached.
 * Note, the recursive bisection initial partitioner is written in TBB continuation style. The TBB
 * continuation style is especially useful for recursive patterns. Each task defines its continuation
 * task. A continuation task defines how computation should continue, if all its child tasks are completed.
 * As a consequence, tasks can be spawned without waiting for their completion, because the continuation
 * task is automatically invoked if all child tasks are terminated. Therefore, no thread will waste CPU
 * time while waiting for their recursive tasks to complete.
 *
 * Implementation Details
 * ----------------------
 * The recursive bisection initial partitioner starts by spawning the root RecursiveMultilevelBisectionTask. The RecursiveMultilevelBisectionTask
 * spawns a MultilevelBisectionTask that bisects the hypergraph (multilevel-fashion). Afterwards, the MultilevelBisectionContinuationTask continues
 * and applies the bisection to the hypergraph and spawns two RecursiveBisectionChildTasks. Both are responsible for exactly one block of
 * the partition. The RecursiveBisectionChildTask extracts its corresponding block as unpartitioned hypergraph and spawns
 * recursively a RecursiveMultilevelBisectionTask for that hypergraph. Once that RecursiveMultilevelBisectionTask is completed, a
 * RecursiveBisectionChildContinuationTask is started and the partition of the recursion is applied to the original hypergraph.
 */
template <typename TypeTraits>
class RecursiveBisectionInitialPartitionerT : public IInitialPartitioner {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using PartitionedHyperGraph = typename TypeTraits::template PartitionedHyperGraph<>;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;

  using BlockRange = std::pair<PartitionID, PartitionID>;

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  static PartitionID kInvalidPartition;
  static HypernodeID kInvalidHypernode;

  struct OriginalHypergraphInfo {

    double computeAdaptiveEpsilon(const HypernodeWeight current_hypergraph_weight,
                                  const PartitionID current_k) const {
      if ( current_hypergraph_weight == 0 ) {
        return 0.0;
      } else {
        double base = ceil(static_cast<double>(original_hypergraph_weight) / original_k)
                      / ceil(static_cast<double>(current_hypergraph_weight) / current_k)
                      * (1.0 + original_epsilon);
        double adaptive_epsilon = std::min(0.99, std::max(std::pow(base, 1.0 /
          ceil(log2(static_cast<double>(current_k)))) - 1.0,0.0));
        return adaptive_epsilon;
      }
    }

    const HypernodeWeight original_hypergraph_weight;
    const PartitionID original_k;
    const double original_epsilon;
  };

  /*!
   * A recursive bisection child task extracts a block of the partition
   * and recursively partition it into the desired number of blocks.
   */
  class RecursiveBisectionChildTask : public tbb::task {

   public:
    RecursiveBisectionChildTask(const OriginalHypergraphInfo original_hypergraph_info,
                                PartitionedHyperGraph& hypergraph,
                                const Context& context,
                                const PartitionID block,
                                const BlockRange range,
                                const TaskGroupID task_group_id) :
      _original_hypergraph_info(original_hypergraph_info),
      _hg(hypergraph),
      _context(context),
      _block(block),
      _range(range),
      _task_group_id(task_group_id) { }

    tbb::task* execute() override {
      const PartitionID k = _range.second - _range.first;
      Context rb_context = setupRecursiveBisectionContext(k);

      // Extracts the block of the hypergraph which we recursively want to partition as
      // seperate unpartitioned hypergraph.
      bool cut_net_splitting = _context.partition.objective == kahypar::Objective::km1;
      auto copy_hypergraph = _hg.extract(_task_group_id, _block, cut_net_splitting);
      HyperGraph& rb_hypergraph = copy_hypergraph.first;
      auto& mapping = copy_hypergraph.second;

      // Spawns a new recursive bisection task to partition the current block of the hypergraph
      // into the desired number of blocks
      if ( rb_hypergraph.initialNumNodes() > 0 ) {
        RecursiveBisectionChildContinuationTask& child_continuation = *new(allocate_continuation())
          RecursiveBisectionChildContinuationTask(_hg, std::move(rb_context),
            std::move(rb_hypergraph), std::move(mapping), k, _task_group_id, _block);
        RecursiveMultilevelBisectionTask& recursion = *new(child_continuation.allocate_child())
          RecursiveMultilevelBisectionTask(
            _original_hypergraph_info,
            child_continuation.recursivePartitionedHypergraph(),
            child_continuation.recursiveContext(), _task_group_id);
        child_continuation.set_ref_count(1);
        tbb::task::spawn(recursion);
      }

      return nullptr;
    }

   private:
    Context setupRecursiveBisectionContext(const PartitionID k) {
      ASSERT(k >= 2);
      Context rb_context(_context);
      rb_context.partition.k = k;

      rb_context.partition.perfect_balance_part_weights.assign(k, 0);
      rb_context.partition.max_part_weights.assign(k, 0);
      for ( PartitionID part_id = _range.first; part_id < _range.second; ++part_id ) {
        rb_context.partition.perfect_balance_part_weights[part_id - _range.first] =
          _context.partition.perfect_balance_part_weights[part_id];
        rb_context.partition.max_part_weights[part_id - _range.first] =
          _context.partition.max_part_weights[part_id];
      }

      return rb_context;
    }

    const OriginalHypergraphInfo _original_hypergraph_info;
    PartitionedHyperGraph& _hg;
    const Context _context;
    const PartitionID _block;
    const BlockRange _range;
    const TaskGroupID _task_group_id;
  };

  /*!
   * Continuation task for the recursive bisection child task. Applies the
   * partition obtained by a recursive bisection task to the original
   * hypergraph.
   */
  class RecursiveBisectionChildContinuationTask : public tbb::task {

    using DeltaFunction = std::function<void (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID)>;
    #define NOOP_FUNC [] (const HyperedgeID, const HyperedgeWeight, const HypernodeID, const HypernodeID, const HypernodeID) { }

   public:
    RecursiveBisectionChildContinuationTask(PartitionedHyperGraph& original_hypergraph,
                                            Context&& context,
                                            HyperGraph&& rb_hypergraph,
                                            parallel::scalable_vector<HypernodeID>&& mapping,
                                            const PartitionID k,
                                            const TaskGroupID task_group_id,
                                            const PartitionID part_id) :
      _original_hg(original_hypergraph),
      _context(std::move(context)),
      _rb_hg(std::move(rb_hypergraph)),
      _rb_partitioned_hg(k, task_group_id, _rb_hg),
      _mapping(std::move(mapping)),
      _task_group_id(task_group_id),
      _part_id(part_id) { }

    tbb::task* execute() override {
      // Applying partition of the recursively bisected hypergraph (_rb_partitioned_hg) to
      // original hypergraph (original_hg). All hypernodes that belong to block
      // 'part_id' in the original hypergraph are moved to the block defined in
      // _rb_partitioned_hg.
      ASSERT(_original_hg.initialNumNodes() == _mapping.size());
      _original_hg.doParallelForAllNodes(_task_group_id, [&](const HypernodeID& hn) {
        if ( _original_hg.partID(hn) == _part_id ) {
          const HypernodeID original_id = _original_hg.originalNodeID(hn);
          ASSERT(original_id < _mapping.size());
          PartitionID to = _part_id + _rb_partitioned_hg.partID(
            _rb_partitioned_hg.globalNodeID(_mapping[original_id]));
          ASSERT(to != kInvalidPartition && to < _original_hg.k());
          if ( _part_id != to ) {
            _original_hg.changeNodePart(hn, _part_id, to, NOOP_FUNC);
          }
        }
      });
      return nullptr;
    }

    PartitionedHyperGraph& recursivePartitionedHypergraph() {
      return _rb_partitioned_hg;
    }

    const Context& recursiveContext() const {
      return _context;
    }

   private:
    PartitionedHyperGraph& _original_hg;
    const Context _context;
    Hypergraph _rb_hg;
    PartitionedHyperGraph _rb_partitioned_hg;
    const parallel::scalable_vector<HypernodeID> _mapping;
    const TaskGroupID _task_group_id;
    const PartitionID _part_id;
  };

  /*!
   * A bisection task bisects the hypergraph into two block. Internally, it
   * calls our multilevel partitioner for k = 2.
   */
  class MultilevelBisectionTask : public tbb::task {
   public:
    MultilevelBisectionTask(PartitionedHyperGraph& original_hypergraph,
                            HyperGraph& bisection_hypergraph,
                            PartitionedHyperGraph& bisection_partitioned_hypergraph,
                            const Context& bisection_context,
                            const TaskGroupID task_group_id) :
      _original_hg(original_hypergraph),
      _bisection_hg(bisection_hypergraph),
      _bisection_partitioned_hg(bisection_partitioned_hypergraph),
      _bisection_context(bisection_context),
      _task_group_id(task_group_id) {}

    tbb::task* execute() override {
      // Bisect hypergraph with parallel multilevel bisection
      _bisection_hg = _original_hg.hypergraph().copy(_task_group_id);
      multilevel::partition_async(_bisection_hg, _bisection_partitioned_hg,
        _bisection_context, false, _task_group_id, this);
      return nullptr;
    }

   private:
    PartitionedHyperGraph& _original_hg;
    HyperGraph& _bisection_hg;
    PartitionedHyperGraph& _bisection_partitioned_hg;
    const Context& _bisection_context;
    const TaskGroupID _task_group_id;
  };

  /*!
   * Continuation task for the bisection task. Applies the bisection
   * obtained by the bisection task and spawns two recursive child
   * bisection task for each block, that further partitions the hypergraph
   * into the desired number of blocks
   */
  class MultilevelBisectionContinuationTask : public tbb::task {

   public:
    MultilevelBisectionContinuationTask(const OriginalHypergraphInfo original_hypergraph_info,
                                        PartitionedHyperGraph& hypergraph,
                                        const Context& context,
                                        const TaskGroupID task_group_id) :
      _bisection_hg(),
      _bisection_partitioned_hg(),
      _bisection_context(),
      _original_hypergraph_info(original_hypergraph_info),
      _hg(hypergraph),
      _context(context),
      _task_group_id(task_group_id) {
      _bisection_context = setupBisectionContext(_hg.hypergraph(), _context);
    }

    tbb::task* execute() override {
      ASSERT(_hg.initialNumNodes() == _bisection_hg.initialNumNodes());
      // Apply partition to hypergraph
      const PartitionID block_0 = 0;
      const PartitionID block_1 = _context.partition.k / 2 + (_context.partition.k % 2 != 0 ? 1 : 0);
      _hg.doParallelForAllNodes(_task_group_id, [&](const HypernodeID& hn) {
        const HypernodeID original_id = _hg.originalNodeID(hn);
        PartitionID part_id = _bisection_partitioned_hg.partID(
          _bisection_partitioned_hg.globalNodeID(original_id));
        ASSERT(part_id != kInvalidPartition && part_id < _hg.k());
        if ( part_id == 0 ) {
          _hg.setOnlyNodePart(hn, block_0);
        } else {
          _hg.setOnlyNodePart(hn, block_1);
        }
      });
      _hg.initializePartition(_task_group_id);

      ASSERT(metrics::objective(_bisection_partitioned_hg, _context.partition.objective) ==
        metrics::objective(_hg, _context.partition.objective));

      ASSERT(_context.partition.k >= 2);
      PartitionID num_blocks_part_0 = _context.partition.k / 2 + (_context.partition.k % 2 != 0 ? 1 : 0);
      PartitionID num_blocks_part_1 = _context.partition.k / 2;
      BlockRange range_0 = std::make_pair(0, num_blocks_part_0);
      BlockRange range_1 = std::make_pair(num_blocks_part_0, num_blocks_part_0 + num_blocks_part_1);

      if ( num_blocks_part_0 >= 2 && num_blocks_part_1 >= 2 ) {
        // In case we have to partition both blocks from the bisection further into
        // more than one block, we call the recursive bisection initial partitioner
        // recursively in parallel
        DBG << "Current k = " << _context.partition.k << "\n"
            << "Parallel Recursion 0: k =" << num_blocks_part_0 << "\n"
            << "Parallel Recursion 1: k =" << num_blocks_part_1;

        auto tbb_recursion_task_groups = TBB::instance().create_tbb_task_groups_for_recursion();
        DoNothingContinuation& recursive_continuation = *new(allocate_continuation()) DoNothingContinuation();
        RecursiveBisectionChildTask& recursion_0 = *new(recursive_continuation.allocate_child()) RecursiveBisectionChildTask(
          _original_hypergraph_info, _hg, _context, 0, range_0, tbb_recursion_task_groups.first);
        RecursiveBisectionChildTask& recursion_1 = *new(recursive_continuation.allocate_child()) RecursiveBisectionChildTask(
          _original_hypergraph_info, _hg, _context, num_blocks_part_0, range_1, tbb_recursion_task_groups.second);
        recursive_continuation.set_ref_count(2);
        tbb::task::spawn(recursion_1);
        tbb::task::spawn(recursion_0);
      } else if ( num_blocks_part_0 >= 2 ) {
        ASSERT(num_blocks_part_1 < 2);
        // In case only the first block has to be partitioned into more than one block, we call
        // the recursive bisection initial partitioner recusively on the block 0
        DBG << "Current k = " << _context.partition.k << ","
            << "Recursion 0: k =" << num_blocks_part_0;
        DoNothingContinuation& recursive_continuation = *new(allocate_continuation()) DoNothingContinuation();
        RecursiveBisectionChildTask& recursion_0 = *new(recursive_continuation.allocate_child()) RecursiveBisectionChildTask(
          _original_hypergraph_info, _hg, _context, 0, range_0, _task_group_id);
        recursive_continuation.set_ref_count(1);
        tbb::task::spawn(recursion_0);
      }

      return nullptr;
    }

    HyperGraph _bisection_hg;
    PartitionedHyperGraph _bisection_partitioned_hg;
    Context _bisection_context;

   private:
    Context setupBisectionContext(const HyperGraph& hypergraph, const Context& context) {
      Context bisection_context(context);

      bisection_context.partition.k = 2;
      bisection_context.partition.verbose_output = debug;
      bisection_context.initial_partitioning.mode = InitialPartitioningMode::direct;
      bisection_context.refinement.flow.algorithm = FlowAlgorithm::do_nothing;
      bisection_context.refinement.label_propagation.numa_aware = false;

      // Setup Part Weights
      HypernodeWeight total_weight = hypergraph.totalWeight();
      if ( context.initial_partitioning.use_adaptive_epsilon ) {
        bisection_context.partition.epsilon = _original_hypergraph_info.computeAdaptiveEpsilon(
          total_weight, context.partition.k);

        bisection_context.partition.perfect_balance_part_weights.clear();
        bisection_context.partition.max_part_weights.clear();
        const PartitionID k = context.partition.k;
        const PartitionID k0 = k / 2 + (k % 2 != 0 ? 1 : 0);
        const PartitionID k1 = k / 2;
        bisection_context.partition.perfect_balance_part_weights.push_back(
          std::ceil(k0 / static_cast<double>(k) * static_cast<double>(total_weight)));
        bisection_context.partition.perfect_balance_part_weights.push_back(
          std::ceil(k1 / static_cast<double>(k) * static_cast<double>(total_weight)));
        bisection_context.partition.max_part_weights.push_back(
          (1 + bisection_context.partition.epsilon) * bisection_context.partition.perfect_balance_part_weights[0]);
        bisection_context.partition.max_part_weights.push_back(
          (1 + bisection_context.partition.epsilon) * bisection_context.partition.perfect_balance_part_weights[1]);
      } else {
        PartitionID num_blocks_part_0 = context.partition.k / 2 + (context.partition.k % 2 != 0 ? 1 : 0);
        ASSERT(num_blocks_part_0 +  context.partition.k / 2 == context.partition.k);
        bisection_context.partition.perfect_balance_part_weights.assign(2, 0);
        bisection_context.partition.max_part_weights.assign(2, 0);
        for ( PartitionID i = 0; i < num_blocks_part_0; ++i ) {
          bisection_context.partition.perfect_balance_part_weights[0] +=
            context.partition.perfect_balance_part_weights[i];
          bisection_context.partition.max_part_weights[0] +=
            context.partition.max_part_weights[i];
        }
        for ( PartitionID i = num_blocks_part_0; i < context.partition.k; ++i ) {
          bisection_context.partition.perfect_balance_part_weights[1] +=
            context.partition.perfect_balance_part_weights[i];
          bisection_context.partition.max_part_weights[1] +=
            context.partition.max_part_weights[i];
        }

        // Special case, if balance constraint will be violated with this bisection
        HypernodeWeight total_max_part_weight = bisection_context.partition.max_part_weights[0] +
          bisection_context.partition.max_part_weights[1];
        if (total_max_part_weight < total_weight) {
          HypernodeWeight delta = total_weight - total_max_part_weight;
          bisection_context.partition.max_part_weights[0] += std::ceil(((double)delta) / 2.0);
          bisection_context.partition.max_part_weights[1] += std::ceil(((double)delta) / 2.0);
        }
      }
      bisection_context.setupContractionLimit(total_weight);
      bisection_context.setupSparsificationParameters();


      return bisection_context;
    }

    const OriginalHypergraphInfo _original_hypergraph_info;
    PartitionedHyperGraph& _hg;
    const Context& _context;
    const TaskGroupID _task_group_id;
  };

  /*!
   * The recursive bisection task spawns the multilevel bisection
   * task and its continuation task.
   */
  class RecursiveMultilevelBisectionTask : public tbb::task {

   public:
    RecursiveMultilevelBisectionTask(const OriginalHypergraphInfo original_hypergraph_info,
                                     PartitionedHyperGraph& hypergraph,
                                     const Context& context,
                                     const TaskGroupID task_group_id) :
      _original_hypergraph_info(original_hypergraph_info),
      _hg(hypergraph),
      _context(context),
      _task_group_id(task_group_id) { }

    tbb::task* execute() override {
      ASSERT(_context.partition.k >= 2);
      MultilevelBisectionContinuationTask& bisection_continuation_task = *new(allocate_continuation())
        MultilevelBisectionContinuationTask(_original_hypergraph_info, _hg, _context, _task_group_id);
      bisection_continuation_task.set_ref_count(1);
      tbb::task::spawn(*new(bisection_continuation_task.allocate_child())
        MultilevelBisectionTask(
          _hg, bisection_continuation_task._bisection_hg,
          bisection_continuation_task._bisection_partitioned_hg,
          bisection_continuation_task._bisection_context,
          _task_group_id));
      return nullptr;
    }

   private:
    const OriginalHypergraphInfo _original_hypergraph_info;
    PartitionedHyperGraph& _hg;
    const Context& _context;
    const TaskGroupID _task_group_id;
  };

  class DoNothingContinuation : public tbb::task {
  public:
    tbb::task* execute() override {
      return nullptr;
    }
  };

 public:
  RecursiveBisectionInitialPartitionerT(PartitionedHyperGraph& hypergraph,
                                        const Context& context,
                                        const bool top_level,
                                        const TaskGroupID task_group_id) :
    _hg(hypergraph),
    _context(context),
    _top_level(top_level),
    _task_group_id(task_group_id) { }

  RecursiveBisectionInitialPartitionerT(const RecursiveBisectionInitialPartitionerT&) = delete;
  RecursiveBisectionInitialPartitionerT(RecursiveBisectionInitialPartitionerT&&) = delete;
  RecursiveBisectionInitialPartitionerT & operator= (const RecursiveBisectionInitialPartitionerT &) = delete;
  RecursiveBisectionInitialPartitionerT & operator= (RecursiveBisectionInitialPartitionerT &&) = delete;

 private:
  void initialPartitionImpl() override final {
    if (_top_level) {
      parallel::MemoryPool::instance().deactivate_unused_memory_allocations();
      utils::Timer::instance().disable();
      utils::Stats::instance().disable();
    }

    RecursiveMultilevelBisectionTask& root_bisection_task = *new(tbb::task::allocate_root()) RecursiveMultilevelBisectionTask(
      OriginalHypergraphInfo { _hg.totalWeight(), _context.partition.k, _context.partition.epsilon },
      _hg, _context, _task_group_id);
    tbb::task::spawn_root_and_wait(root_bisection_task);

    if (_top_level) {
      parallel::MemoryPool::instance().activate_unused_memory_allocations();
      utils::Timer::instance().enable();
      utils::Stats::instance().enable();
    }
  }

  PartitionedHyperGraph& _hg;
  const Context& _context;
  const bool _top_level;
  const TaskGroupID _task_group_id;
};

template <typename TypeTraits>
PartitionID RecursiveBisectionInitialPartitionerT<TypeTraits>::kInvalidPartition = -1;
template <typename TypeTraits>
HypernodeID RecursiveBisectionInitialPartitionerT<TypeTraits>::kInvalidHypernode = std::numeric_limits<HypernodeID>::max();

using RecursiveBisectionInitialPartitioner = RecursiveBisectionInitialPartitionerT<GlobalTypeTraits>;
}  // namespace mt_kahypar
