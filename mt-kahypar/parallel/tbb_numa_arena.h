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

#include <mutex>

#include "tbb/task_arena.h"
#include "tbb/task_group.h"
#include "tbb/task_scheduler_init.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/global_thread_pinning_observer.h"
#include "mt-kahypar/parallel/numa_thread_pinning_observer.h"

namespace mt_kahypar {
namespace parallel {
/**
 * Creates number of NUMA nodes TBB task arenas. Each task arena is pinned
 * to a unique NUMA node. Each task arena can then be used to execute tasks
 * on specific NUMA node.
 */
template <typename HwTopology>
class TBBNumaArena {
  static constexpr bool debug = false;

 private:
  using GlobalThreadPinning = mt_kahypar::parallel::GlobalThreadPinning<HwTopology>;
  using GlobalThreadPinningObserver = mt_kahypar::parallel::GlobalThreadPinningObserver<HwTopology>;
  using NumaThreadPinningObserver = mt_kahypar::parallel::NumaThreadPinningObserver<HwTopology>;

 public:
  TBBNumaArena(const TBBNumaArena&) = delete;
  TBBNumaArena & operator= (const TBBNumaArena &) = delete;

  TBBNumaArena(TBBNumaArena&&) = delete;
  TBBNumaArena & operator= (TBBNumaArena &&) = delete;

  static TBBNumaArena& instance(const size_t num_threads = 1) {
    static TBBNumaArena instance(num_threads);
    return instance;
  }

  int total_number_of_threads() const {
    return _num_threads;
  }

  int number_of_threads_on_numa_node(const int node) const {
    ASSERT(node < (int)_arenas.size());
    return _arenas[node].max_concurrency();
  }

  int num_used_numa_nodes() const {
    return _arenas.size();
  }

  tbb::task_arena& numa_task_arena(const int node) {
    ASSERT(node < (int)_arenas.size());
    return _arenas[node];
  }

  tbb::task_group& numa_task_group(const int node) {
    ASSERT(_arenas.size() <= _groups.size());
    ASSERT(node < (int)_arenas.size());
    return _groups[node];
  }

  template <typename F>
  void execute_sequential_on_all_numa_nodes(F&& func) {
    for (int node = 0; node < num_used_numa_nodes(); ++node) {
      numa_task_arena(node).execute([&] {
            numa_task_group(node).run([&, node] {
              func(node);
            });
          });
      wait(node, numa_task_group(node));
    }
  }

  template <typename F>
  void execute_parallel_on_all_numa_nodes(F&& func) {
    for (int node = 0; node < num_used_numa_nodes(); ++node) {
      numa_task_arena(node).execute([&] {
            numa_task_group(node).run([&, node] {
              func(node);
            });
          });
    }
    wait();
  }

  void wait() {
    for (int node = 0; node < num_used_numa_nodes(); ++node) {
      _arenas[node].execute([&, node] {
            _groups[node].wait();
          });
    }
  }

  void wait(const int node, tbb::task_group& group) {
    ASSERT(node < (int)_arenas.size());
    _arenas[node].execute([&] {
          group.wait();
        });
  }

  void terminate() {
    for (NumaThreadPinningObserver& observer : _observer) {
      observer.observe(false);
    }
    _global_observer.observe(false);

    for (tbb::task_arena& arena : _arenas) {
      arena.terminate();
    }
    _init.terminate();
  }

 private:
  explicit TBBNumaArena(const int num_threads) :
    _num_threads(num_threads),
    _init(num_threads),
    _arenas(),
    _groups(HwTopology::instance().num_numa_nodes()),
    _global_observer(),
    _observer() {
    HwTopology& topology = HwTopology::instance();
    int threads_left = num_threads;
    int num_numa_nodes = topology.num_numa_nodes();
    DBG << "Initialize TBB with" << num_threads << "threads";
    _arenas.reserve(num_numa_nodes);
    // TODO(heuer): fix copy constructor of observer
    _observer.reserve(num_numa_nodes);

    std::vector<int> used_cpus_on_numa_node(num_numa_nodes);
    // First use cores (to prevent using hyperthreads)
    for (int node = 0; node < num_numa_nodes && threads_left > 0; ++node) {
      used_cpus_on_numa_node[node] = std::min(threads_left, topology.num_cores_on_numa_node(node));
      threads_left -= used_cpus_on_numa_node[node];
    }

    // If there are still thread to assign left we use hyperthreading
    for (int node = 0; node < num_numa_nodes && threads_left > 0; ++node) {
      int num_hyperthreads = std::min(threads_left, topology.num_cpus_on_numa_node(node) - used_cpus_on_numa_node[node]);
      used_cpus_on_numa_node[node] += num_hyperthreads;
      threads_left -= num_hyperthreads;
    }

    for (int node = 0; node < num_numa_nodes; ++node) {
      if (used_cpus_on_numa_node[node] > 0) {
        int num_cpus = used_cpus_on_numa_node[node];
        DBG << "Initialize TBB task arena on numa node" << node
            << "with" << num_cpus << "threads";
        _arenas.emplace_back(num_cpus, num_cpus == 1 ? 1 : 0);
        _observer.emplace_back(_arenas.back(), node);
      }
    }
    _num_threads -= threads_left;

    // Initialize Global Thread Pinning
    GlobalThreadPinning::instance(num_threads);
    for (int node = 0; node < num_numa_nodes; ++node) {
      // Seems that there is one extra worker threads when num_threads is equal to,
      // but only one will participate in task scheduling at a time
      int num_cpus = std::max(used_cpus_on_numa_node[node], 2);
      topology.use_only_num_cpus_on_numa_node(node, num_cpus);
    }
    _global_observer.observe(true);
  }

  int _num_threads;
  tbb::task_scheduler_init _init;
  std::vector<tbb::task_arena> _arenas;
  std::vector<tbb::task_group> _groups;
  GlobalThreadPinningObserver _global_observer;
  std::vector<NumaThreadPinningObserver> _observer;
};
}  // namespace parallel
}  // namespace mt_kahypar
