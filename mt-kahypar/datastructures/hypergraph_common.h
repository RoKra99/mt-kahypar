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
#if KAHYPAR_USE_64_BIT_IDS
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

// Flow Network Types
using ClusterID = PartitionID;
using Flow = int32_t;
using Capacity = int32_t;

// Graph Types
using NodeID = uint32_t;
using ArcWeight = double;

struct Arc {
  NodeID head;
  ArcWeight weight;

  Arc() :
    head(0),
    weight(0) { }

  Arc(NodeID head, ArcWeight weight) :
    head(head),
    weight(weight) { }
};

// Constant Declarations
static constexpr PartitionID kInvalidPartition = -1;
static constexpr HypernodeID kInvalidHypernode = std::numeric_limits<HypernodeID>::max();
static constexpr HypernodeID kInvalidHyperedge = std::numeric_limits<HyperedgeID>::max();
static constexpr NodeID kinvalidFlowNetworkNode = std::numeric_limits<NodeID>::max() / 2;
static constexpr Flow kInfty = std::numeric_limits<Flow>::max() / 2;
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
  struct HyperedgeHash {
    HyperedgeID he;
    size_t hash;
    size_t size;
    bool valid;
  };

} // namespace mt_kahypar
