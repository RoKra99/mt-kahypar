/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2018 Sebastian Schlag <sebastian.schlag@kit.edu>
 * Copyright (C) 2018 Tobias Heuer <tobias.heuer@live.com>
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
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

#include "external_maximum_flow/bk/graph.h"
#include "external_maximum_flow/ibfs/ibfs.h"

#include "kahypar/datastructure/fast_reset_array.h"
#include "kahypar/datastructure/fast_reset_flag_array.h"
#include "kahypar/datastructure/sparse_set.h"
#include "kahypar/meta/mandatory.h"
#include "kahypar/meta/typelist.h"

#include "mt-kahypar/datastructures/flow_network.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/randomize.h"
#include "mt-kahypar/partition/refinement/flow/most_balanced_minimum_cut.h"

namespace mt_kahypar {
template <typename TypeTraits, class Network = Mandatory>
class MaximumFlow {
using HyperGraph = typename TypeTraits::HyperGraph;
using FlowNetwork = ds::FlowNetwork<TypeTraits>;

 public:
  MaximumFlow(const HypernodeID initial_size, const HypernodeID initial_num_nodes) :
    _parent(initial_size, nullptr),
    _visited(initial_size),
    _Q(),
    _mbmc(initial_size),
    _original_part_id(initial_num_nodes, 0) { }

  virtual ~MaximumFlow() { }

  MaximumFlow(const MaximumFlow&) = delete;
  MaximumFlow& operator= (const MaximumFlow&) = delete;

  MaximumFlow(MaximumFlow&&) = delete;
  MaximumFlow& operator= (MaximumFlow&&) = delete;


  virtual Flow maximumFlow(HyperGraph& hypergraph, FlowNetwork& flow_network) = 0;

  HyperedgeWeight minimumSTCut(HyperGraph& hypergraph, FlowNetwork& flow_network,
                               const Context& context,
                               const PartitionID block_0, const PartitionID block_1) {
    if (flow_network.isTrivialFlow()) {
      return Network::kInfty;
    }

    const PartitionID default_part =
      context.refinement.flow.use_most_balanced_minimum_cut ? block_0 : block_1;
    for (const HypernodeID& ogHn : flow_network.hypernodes()) {
      const HypernodeID& hn = hypergraph.globalNodeID(ogHn);
      _original_part_id[ogHn] = hypergraph.partID(hn);
      moveHypernode(hypergraph, hn, default_part);
    }

    const HyperedgeWeight cut = maximumFlow(hypergraph, flow_network);

    if (context.refinement.flow.use_most_balanced_minimum_cut) {
      _mbmc.mostBalancedMinimumCut(hypergraph, flow_network, context, block_0, block_1);
    } else {
      bfs<true>(hypergraph, flow_network, block_0);
    }

    return cut;
  }

  void rollback(HyperGraph& hypergraph, FlowNetwork& flow_network, const bool store_part_id = false) {
    for (const HypernodeID& ogHn : flow_network.hypernodes()) {
      const HypernodeID& hn = hypergraph.globalNodeID(ogHn);
      const PartitionID from = hypergraph.partID(hn);
      moveHypernode(hypergraph, hn, _original_part_id[ogHn]);
      if (store_part_id) {
        _original_part_id[ogHn] = from;
      }
    }
  }

  PartitionID getOriginalPartition(const HypernodeID hn_og) const {
    return _original_part_id[hn_og];
  }

  template <bool assign_hypernodes = false>
  bool bfs(HyperGraph& hypergraph, FlowNetwork& flow_network, const PartitionID block = 0) {
    bool augmenting_path_exists = false;
    _parent.resetUsedEntries();
    _visited.reset();
    while (!_Q.empty()) {
      _Q.pop();
    }

    // Initialize queue with all source nodes
    for (const NodeID& s : flow_network.sources()) {
      _visited.set(s, true);
      _parent.set(s, nullptr);
      _Q.push(s);
    }

    while (!_Q.empty()) {
      NodeID u_og = _Q.front();
      _Q.pop();

      if (assign_hypernodes) {
        if (flow_network.interpreteHypernode(u_og)) {
          
          moveHypernode(hypergraph, u_og, block);
        } else if (flow_network.interpreteHyperedge(u_og)) {
          const HyperedgeID he_og = flow_network.mapToHyperedgeID(u_og);
          const HyperedgeID he = hypergraph.globalNodeID(he_og);
          for (const HypernodeID& pin : hypergraph.pins(he)) {
            if (flow_network.containsHypernode(hypergraph, pin)) {
              moveHypernode(hypergraph, pin, block);
            }
          }
        }
      }

      if (flow_network.isIdSink(u_og)) {
        augmenting_path_exists = true;
        continue;
      }

      for (ds::FlowEdge& e : flow_network.incidentEdges(u_og)) {
        const NodeID v = e.target;
        if (!_visited[v] && flow_network.residualCapacity(e)) {
          _parent.set(v, &e);
          _visited.set(v, true);
          _Q.push(v);
        }
      }
    }
    return augmenting_path_exists;
  }

 protected:
  template <typename T>
  FRIEND_TEST(AMaximumFlow, ChecksIfAugmentingPathExist);
  template <typename T>
  FRIEND_TEST(AMaximumFlow, AugmentAlongPath);

  Flow augment(FlowNetwork& flow_network, const NodeID cur, const Flow min_flow = Network::kInfty) {
    if (flow_network.isGlobalIdSource(cur) || min_flow == 0) {
      return min_flow;
    } else {
      ds::FlowEdge* e = _parent.get(cur);
      const Flow f = augment(e->source, std::min(min_flow, flow_network.residualCapacity(*e)));

      ASSERT([&]() {
          const Flow residual_forward_before = flow_network.residualCapacity(*e);
          const Flow residual_backward_before = flow_network.residualCapacity(flow_network.reverseEdge(*e));
          flow_network.increaseFlow(*e, f);
          Flow residual_forward_after = flow_network.residualCapacity(*e);
          Flow residual_backward_after = flow_network.residualCapacity(flow_network.reverseEdge(*e));
          if (residual_forward_before != Network::kInfty && residual_forward_before != residual_forward_after + f) {
            LOG << "Residual capacity should be " << (residual_forward_before - f) << "!";
            return false;
          }
          if (residual_backward_before != Network::kInfty && residual_backward_before != residual_backward_after - f) {
            LOG << "Residual capacity should be " << (residual_backward_before + f) << "!";
            return false;
          }
          flow_network.increaseFlow(flow_network.reverseEdge(*e), f);
          residual_forward_after = flow_network.residualCapacity(*e);
          residual_backward_after = flow_network.residualCapacity(flow_network.reverseEdge(*e));
          if (residual_forward_before != residual_forward_after ||
              residual_backward_before != residual_backward_after) {
            LOG << "Restoring original capacities failed!";
            return false;
          }
          return true;
        } (), "Flow is not increased correctly!");

      flow_network.increaseFlow(*e, f);
      return f;
    }
  }

  void moveHypernode(HyperGraph& hypergraph, const HypernodeID hn, const PartitionID to) {
    ASSERT(hypergraph.partID(hn) != -1, "Hypernode " << hn << " should be assigned to a part");
    const PartitionID from = hypergraph.partID(hn);
    if (from != to ) { //&& !hypergraph.isFixedVertex(hn)
      bool success = hypergraph.changeNodePart(hn, from, to);
      unused(success);
      ASSERT(success);
    }
  }

  // Datastructure for BFS
  kahypar::ds::FastResetArray<mt_kahypar::ds::FlowEdge*> _parent;
  kahypar::ds::FastResetFlagArray<> _visited;
  std::queue<NodeID> _Q;

  MostBalancedMinimumCut<TypeTraits, Network> _mbmc;

  std::vector<PartitionID> _original_part_id;
};

template <typename TypeTraits, class Network = Mandatory>
class BoykovKolmogorov : public MaximumFlow<TypeTraits, Network>{
  using Base = MaximumFlow<TypeTraits,Network>;
  using FlowGraph = maxflow::Graph<int, int, int>;
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using FlowNetwork = ds::FlowNetwork<TypeTraits>;

 public:
  BoykovKolmogorov(const HypernodeID initial_size, const HypernodeID initial_num_nodes) :
    Base(initial_size, initial_num_nodes),
    _flow_graph(initial_size, initial_size),
    _flow_network_mapping(initial_size, 0) { }

  ~BoykovKolmogorov() = default;

  BoykovKolmogorov(const BoykovKolmogorov&) = delete;
  BoykovKolmogorov& operator= (const BoykovKolmogorov&) = delete;

  BoykovKolmogorov(BoykovKolmogorov&&) = delete;
  BoykovKolmogorov& operator= (BoykovKolmogorov&&) = delete;

  Flow maximumFlow(HyperGraph& hypergraph, FlowNetwork& flow_network) {
    unused(hypergraph);
    mapToExternalFlowNetwork();

    const Flow max_flow = _flow_graph.maxflow();

    FlowGraph::arc* a = _flow_graph.get_first_arc();
    while (a != _flow_graph.arc_last) {
      const Flow flow = a->flowEdge->capacity - _flow_graph.get_rcap(a);
      if (flow != 0) {
        a->flowEdge->increaseFlow(flow);
      }
      a = _flow_graph.get_next_arc(a);
    }

    ASSERT(!Base::bfs(hypergraph, flow_network), "Found augmenting path after flow computation finished!");
    return max_flow;
  }

 private:
  template <typename T>
  FRIEND_TEST(AMaximumFlow, ChecksIfAugmentingPathExist);
  template <typename T>
  FRIEND_TEST(AMaximumFlow, AugmentAlongPath);

  void mapToExternalFlowNetwork(FlowNetwork& flow_network) {
    _flow_graph.reset();
    _visited.reset();
    const Flow infty = flow_network.totalWeightHyperedges();

    for (const NodeID& node : flow_network.nodes()) {
      const NodeID id = _flow_graph.add_node();
      _flow_network_mapping[node] = id;
      if (flow_network.isGlobalIdSource(node)) {
        _flow_graph.add_tweights(id, infty, 0);
      }
      if (flow_network.isGlobalIdSink(node)) {
        _flow_graph.add_tweights(id, 0, infty);
      }
    }

    for (const NodeID node : flow_network.nodes()) {
      const NodeID u = _flow_network_mapping[node];
      for (ds::FlowEdge& edge : flow_network.incidentEdges(node)) {
        const NodeID v = _flow_network_mapping[edge.target];
        const Capacity c = edge.capacity;
        ds::FlowEdge& rev_edge = flow_network.reverseEdge(edge);
        const Capacity rev_c = rev_edge.capacity;
        if (!_visited[edge.target]) {
          FlowGraph::arc* a = _flow_graph.add_edge(u, v, c, rev_c);
          a->flowEdge = &edge;
          a->sister->flowEdge = &rev_edge;
        }
      }
      _visited.set(node, true);
    }
  }

  using Base::_parent;
  using Base::_visited;

  FlowGraph _flow_graph;
  std::vector<NodeID> _flow_network_mapping;
};


template <typename TypeTraits, class Network = Mandatory>
class IBFS : public MaximumFlow<TypeTraits, Network>{
  using Base = MaximumFlow<TypeTraits, Network>;
  using FlowGraph = maxflow::IBFSGraph;
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using FlowNetwork = ds::FlowNetwork<TypeTraits>;

 public:
  IBFS(const HypernodeID initial_size, const HypernodeID initial_num_nodes) :
    Base(initial_size, initial_num_nodes),
    _flow_graph(FlowGraph::IB_INIT_COMPACT),
    _flow_network_mapping(initial_size, 0) { }

  ~IBFS() = default;

  IBFS(const IBFS&) = delete;
  IBFS& operator= (const IBFS&) = delete;

  IBFS(IBFS&&) = delete;
  IBFS& operator= (IBFS&&) = delete;

  Flow maximumFlow(HyperGraph& hypergraph, FlowNetwork& flow_network) {
    unused(hypergraph);
    mapToExternalFlowNetwork(flow_network);

    _flow_graph.computeMaxFlow();
    const Flow max_flow = _flow_graph.getFlow();

    FlowGraph::Arc* a = _flow_graph.arcs;
    while (a != _flow_graph.arcEnd) {
      const Flow flow = a->flowEdge->capacity - a->rCap;
      if (flow != 0) {
        a->flowEdge->increaseFlow(flow);
      }
      a++;
    }

    ASSERT(!Base::bfs(hypergraph, flow_network), "Found augmenting path after flow computation finished!");
    return max_flow;
  }

 private:
  template <typename T>
  FRIEND_TEST(AMaximumFlow, ChecksIfAugmentingPathExist);
  template <typename T>
  FRIEND_TEST(AMaximumFlow, AugmentAlongPath);

  void mapToExternalFlowNetwork(FlowNetwork& flow_network) {
    _flow_graph.initSize(flow_network.numNodes(), flow_network.numEdges() - flow_network.numUndirectedEdges());
    _visited.reset();
    const Flow infty = flow_network.totalWeightHyperedges();
    NodeID cur_id = 0;

    for (const NodeID& node : flow_network.nodes()) {
      _flow_graph.addNode(cur_id,
                          flow_network.isIdSource(node) ? infty : 0,
                          flow_network.isIdSink(node) ? infty : 0);
      _flow_network_mapping[node] = cur_id;
      cur_id++;
    }

    for (const NodeID node : flow_network.nodes()) {
      const NodeID u = _flow_network_mapping[node];
      for (ds::FlowEdge& edge : flow_network.incidentEdges(node)) {
        const NodeID v = _flow_network_mapping[edge.target];
        const Capacity c = edge.capacity;
        ds::FlowEdge& rev_edge = flow_network.reverseEdge(edge);
        const Capacity rev_c = rev_edge.capacity;
        if (!_visited[edge.target]) {
          _flow_graph.addEdge(u, v, c, rev_c, &edge, &rev_edge);
        }
      }
      _visited.set(node, true);
    }

    _flow_graph.initGraph();
  }

  using Base::_parent;
  using Base::_visited;

  FlowGraph _flow_graph;
  std::vector<NodeID> _flow_network_mapping;
};
}  // namespace mt_kahypar