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

#include <cstdint>

#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {

#define UI64(X) static_cast<uint64_t>(X)

using TaskGroupID = size_t;
using RatingType = double;
#if ( KAHYPAR_ENABLE_NUMA_AWARE_PARTITIONING || HYPERGRAPH_UNIT_TEST )
#define ID(X) static_cast<uint64_t>(X)
using HypernodeID = uint64_t;
using HyperedgeID = uint64_t;
#else
#define ID(X) static_cast<uint32_t>(X)
using HypernodeID = uint32_t;
using HyperedgeID = uint32_t;
#endif
using HypernodeWeight = int32_t;
using HyperedgeWeight = int32_t;
using PartitionID = int32_t;
using Gain = HyperedgeWeight;
using ClusterID = PartitionID;
using Flow = int32_t;
using Capacity = int32_t;

// Constant Declarations
static constexpr PartitionID kInvalidPartition = -1;
static constexpr HypernodeID kInvalidHypernode = std::numeric_limits<HypernodeID>::max();
static constexpr HypernodeID kInvalidHyperedge = std::numeric_limits<HyperedgeID>::max();
static constexpr size_t kEdgeHashSeed = 42;

struct Move {
  PartitionID from;
  PartitionID to;
  Gain gain;
};

/*!
* A memento stores all information necessary to undo the contraction operation
* of a vertex pair \f$(u,v)\f$.
*/
struct Memento {
  Memento() :
    u(kInvalidHypernode),
    v(kInvalidHypernode),
    community_id(kInvalidPartition),
    one_pin_hes_begin(0),
    one_pin_hes_size(0),
    parallel_hes_begin(0),
    parallel_hes_size(0) { }

  Memento(HypernodeID representative, HypernodeID contraction_partner) :
    u(representative),
    v(contraction_partner),
    community_id(kInvalidPartition),
    one_pin_hes_begin(0),
    one_pin_hes_size(0),
    parallel_hes_begin(0),
    parallel_hes_size(0) { }

  Memento(HypernodeID representative, HypernodeID contraction_partner, PartitionID community) :
    u(representative),
    v(contraction_partner),
    community_id(community),
    one_pin_hes_begin(0),
    one_pin_hes_size(0),
    parallel_hes_begin(0),
    parallel_hes_size(0) { }

  // ! The representative hypernode that remains in the hypergraph
  HypernodeID u;
  // ! The contraction partner of u that is removed from the hypergraph after the contraction.
  HypernodeID v;
  // ! Community id of hypernodes
  PartitionID community_id;
  // ! start of removed single pin hyperedges
  int one_pin_hes_begin;
  // ! # removed single pin hyperedges
  int one_pin_hes_size;
  // ! start of removed parallel hyperedges
  int parallel_hes_begin;
  // ! # removed parallel hyperedges
  int parallel_hes_size;
};

/*!
  * This struct is used during multilevel coarsening to efficiently
  * detect parallel hyperedges.
  */
struct ContractedHyperedge {
  HyperedgeID he;
  size_t hash;
  HyperedgeWeight weight;
  bool is_parallel; // Indicates if this hyperedges is already detected as parallel
  int node; // NUMA Node of hyperedge
  parallel::scalable_vector<HypernodeID> hyperedge;
  HyperedgeID he_idx; // Index in hyperedge vector
  HypernodeID pin_idx; // Index of pins in incidence array
};

namespace common {

#if ( KAHYPAR_ENABLE_NUMA_AWARE_PARTITIONING || HYPERGRAPH_UNIT_TEST )

static constexpr size_t NUMA_NODE_IDENTIFIER = 48;

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HypernodeID get_global_vertex_id(const int node, const size_t vertex_pos) {
  return ( static_cast<HypernodeID>(node) << NUMA_NODE_IDENTIFIER ) | vertex_pos;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HypernodeID get_local_position_of_vertex(const HypernodeID u) {
  return ((1UL << NUMA_NODE_IDENTIFIER) - 1) & u;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static int get_numa_node_of_vertex(const HypernodeID u) {
  return static_cast<int>(u >> NUMA_NODE_IDENTIFIER);
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HyperedgeID get_global_edge_id(const int node, const size_t edge_pos) {
  return ( static_cast<HypernodeID>(node) << NUMA_NODE_IDENTIFIER ) | edge_pos;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HyperedgeID get_local_position_of_edge(const HyperedgeID e) {
  return ((1UL << NUMA_NODE_IDENTIFIER) - 1) & e;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static int get_numa_node_of_edge(const HyperedgeID e) {
  return static_cast<int>(e >> NUMA_NODE_IDENTIFIER);
}
#else
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HypernodeID get_global_vertex_id(const int, const size_t vertex_pos) {
  return vertex_pos;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HypernodeID get_local_position_of_vertex(const HypernodeID u) {
  return u;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static int get_numa_node_of_vertex(const HypernodeID) {
  return 0;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HyperedgeID get_global_edge_id(const int, const size_t edge_pos) {
  return edge_pos;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static HyperedgeID get_local_position_of_edge(const HyperedgeID e) {
  return e;
}

MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static int get_numa_node_of_edge(const HyperedgeID) {
  return 0;
}
#endif

template<typename Hypergraph>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static const Hypergraph& hypergraph_of_vertex(
  const HypernodeID u,
  const parallel::scalable_vector<Hypergraph>& hypergraphs) {
  int node = get_numa_node_of_vertex(u);
  ASSERT(node < static_cast<int>(hypergraphs.size()));
  return hypergraphs[node];
}

template<typename Hypergraph>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static Hypergraph& hypergraph_of_vertex(
  const HypernodeID u,
  parallel::scalable_vector<Hypergraph>& hypergraphs) {
  int node = get_numa_node_of_vertex(u);
  ASSERT(node < static_cast<int>(hypergraphs.size()));
  return hypergraphs[node];
}

template<typename Hypergraph>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static const Hypergraph& hypergraph_of_edge(
  const HyperedgeID e,
  const parallel::scalable_vector<Hypergraph>& hypergraphs) {
  int node = get_numa_node_of_edge(e);
  ASSERT(node < static_cast<int>(hypergraphs.size()));
  return hypergraphs[node];
}

template<typename Hypergraph>
MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE static Hypergraph& hypergraph_of_edge(
  const HyperedgeID e,
  parallel::scalable_vector<Hypergraph>& hypergraphs) {
  int node = get_numa_node_of_edge(e);
  ASSERT(node < static_cast<int>(hypergraphs.size()));
  return hypergraphs[node];
}

} // namespace common

} // namespace mt_kahypar
