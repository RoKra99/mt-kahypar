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

#define USE_HARDWARE_MOCK false

#include <chrono>

#include "tbb/enumerable_thread_specific.h"

#include "kahypar/datastructure/fast_reset_flag_array.h"

#include "mt-kahypar/datastructures/hypergraph.h"
#include "mt-kahypar/datastructures/streaming_hypergraph.h"
#include "mt-kahypar/parallel/hardware_topology.h"
#include "mt-kahypar/parallel/tbb_numa_arena.h"

#include "tests/parallel/topology_mock.h"

namespace mt_kahypar {
#if USE_HARDWARE_MOCK
static constexpr int NUM_NUMA_NODES = 2;
using TopoMock = mt_kahypar::parallel::TopologyMock<NUM_NUMA_NODES>;
using topology_t = mt_kahypar::parallel::topology_t;
using node_t = mt_kahypar::parallel::node_t;
using HardwareTopology = mt_kahypar::parallel::HardwareTopology<TopoMock, topology_t, node_t>;
#else
using HardwareTopology = mt_kahypar::parallel::HardwareTopology<>;
#endif
using TBBNumaArena = mt_kahypar::parallel::TBBNumaArena<HardwareTopology>;
using ThreadLocalFastResetFlagArray = tbb::enumerable_thread_specific<kahypar::ds::FastResetFlagArray<> >;

using RatingType = double;
using HypernodeID = uint64_t;
using HyperedgeID = uint64_t;
using HypernodeWeight = int32_t;
using HyperedgeWeight = int32_t;
using PartitionID = int32_t;
using Gain = HyperedgeWeight;

// Note(gottesbueren) : why not keep it 32 bits until we are actually ready for those types of instances. Especially for the comparison against sequential KaHyPar?
using NodeID = uint32_t;

// #########Graph-Definitions#############
using EdgeID = uint32_t;
using EdgeWeight = long double;
using ClusterID = PartitionID;
using Flow = int32_t;
using Capacity = int32_t;


struct Move {
  PartitionID from;
  PartitionID to;
  Gain gain;
};

using StreamingHypergraph = mt_kahypar::ds::StreamingHypergraph<HypernodeID,
                                                                HyperedgeID,
                                                                HypernodeWeight,
                                                                HyperedgeWeight,
                                                                PartitionID,
                                                                HardwareTopology,
                                                                TBBNumaArena>;

using Hypergraph = mt_kahypar::ds::Hypergraph<HypernodeID,
                                              HyperedgeID,
                                              HypernodeWeight,
                                              HyperedgeWeight,
                                              PartitionID,
                                              HardwareTopology,
                                              TBBNumaArena>;

struct GlobalTypeTraits {
  using HyperGraph = Hypergraph;
  using StreamingHyperGraph = StreamingHypergraph;
  using TBB = TBBNumaArena;
  using HwTopology = HardwareTopology;
};

using HighResClockTimepoint = std::chrono::time_point<std::chrono::high_resolution_clock>;
}  // namespace mt_kahypar
