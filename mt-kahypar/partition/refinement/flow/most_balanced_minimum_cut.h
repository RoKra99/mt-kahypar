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
#include <array>
#include <queue>
#include <stack>
#include <string>
#include <utility>
#include <vector>

#include "external_tools/kahypar/kahypar/datastructure/fast_reset_array.h"
#include "external_tools/kahypar/kahypar/datastructure/fast_reset_flag_array.h"
#include "external_tools/kahypar/kahypar/datastructure/sparse_set.h"
#include "external_tools/kahypar/kahypar/datastructure/graph.h"

#include "mt-kahypar/datastructures/flow_network.h"


#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/refinement/flow/strongly_connected_components.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
using kahypar::ds::Graph;
using kahypar::ds::Edge;

template <typename TypeTraits, class Network = Mandatory>
class MostBalancedMinimumCut {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using FlowNetwork = ds::FlowNetwork<TypeTraits>;

 public:
  MostBalancedMinimumCut(const size_t initial_size) :
    _visited(initial_size),
    _graph_to_flow_network(initial_size, Network::kInvalidNode),
    _flow_network_to_graph(initial_size, Network::kInvalidNode),
    _scc_node_weight(initial_size, 0),
    _Q(),
    _sccs(initial_size) { }

  MostBalancedMinimumCut(const MostBalancedMinimumCut&) = delete;
  MostBalancedMinimumCut(MostBalancedMinimumCut&&) = delete;
  MostBalancedMinimumCut& operator= (const MostBalancedMinimumCut&) = delete;
  MostBalancedMinimumCut& operator= (MostBalancedMinimumCut&&) = delete;

  void mostBalancedMinimumCut(HyperGraph& hypergraph, FlowNetwork& flow_network,
                              const Context& context,
                              const PartitionID block_0, const PartitionID block_1) {
    reset();

    // Mark all reachable nodes from source and sink set as invalid
    markAllReachableNodesAsVisited<true>(hypergraph, flow_network, block_0, block_1);
    markAllReachableNodesAsVisited<false>(hypergraph, flow_network, block_0, block_1);

    // Build residual graph
    Graph residual_graph(std::move(buildResidualGraph(hypergraph, flow_network)));

    // Find strongly connected components
    findStronglyConnectedComponents(residual_graph);

    // Contract strongly connected components
    auto contraction = residual_graph.contractClusters();
    Graph dag(std::move(contraction.first));
    std::vector<NodeID> contraction_mapping = std::move(contraction.second);

    // Build mapping from contracted graph to flow network
    std::vector<std::vector<NodeID> > scc_to_flow_network(dag.numNodes(), std::vector<NodeID>());
    for (const NodeID& u : residual_graph.nodes()) {
      const NodeID flow_u = _graph_to_flow_network.get(u);
      if (flow_network.isHypernode(flow_u)) {
        scc_to_flow_network[contraction_mapping[u]].push_back(flow_u);
        _scc_node_weight.update(contraction_mapping[u], hypergraph.nodeWeight(flow_u));
      }
    }

    // Calculate in degrees of nodes in DAG graph for topological ordering
    std::vector<size_t> in_degree(dag.numNodes(), 0);
    for (const NodeID& u : dag.nodes()) {
      for (const Edge& e : dag.incidentEdges(u)) {
        const NodeID v = e.target_node;
        if (u != v) {
          in_degree[v]++;
        }
      }
    }

    // Find most balanced minimum cut
    std::vector<NodeID> topological_order(dag.numNodes(), 0);
    std::vector<PartitionID> best_partition_id(dag.numNodes(), block_0);
    double best_imbalance = metrics::localImbalance(hypergraph, context);;

    DBG << "Start Most Balanced Minimum Cut (Bipartition = {" << block_0 << "," << block_1 << "}";
    DBG << "Initial imbalance: " << V(metrics::localImbalance(hypergraph, context));

    for (size_t i = 0; i < 20; ++i) {
      // Compute random topological order
      topologicalSort(dag, in_degree, topological_order);

      // Sweep through topological order and find best imbalance
      std::vector<PartitionID> tmp_partition_id(dag.numNodes(), block_0);
      double tmp_best_imbalance = metrics::localImbalance(hypergraph, context);

      std::vector<HypernodeWeight> part_weight(context.partition.k, 0);
      for (PartitionID part = 0; part < context.partition.k; ++part) {
        part_weight[part] = hypergraph.localPartWeight(part);
      }
      for (size_t idx = 0; idx < topological_order.size(); ++idx) {
        const NodeID u = topological_order[idx];
        tmp_partition_id[u] = block_1;
        part_weight[block_0] -= _scc_node_weight.get(u);
        part_weight[block_1] += _scc_node_weight.get(u);
        double cur_imbalance = imbalance<true>(context, part_weight);
        if (context.partition.k > 2) {
          cur_imbalance = imbalance<false>(context, part_weight);
        }

        if (cur_imbalance > tmp_best_imbalance) {
          tmp_partition_id[u] = block_0;
          break;
        }
        tmp_best_imbalance = cur_imbalance;
      }

      if (tmp_best_imbalance < best_imbalance) {
        best_imbalance = tmp_best_imbalance;
        best_partition_id = tmp_partition_id;
      }
    }

    DBG << "Best imbalance: " << best_imbalance;

    ASSERT([&]() {
        //const HyperedgeWeight metric_before = metrics::objective(hypergraph, context.partition.objective);
        const double imbalance_before = metrics::localImbalance(hypergraph, context);
        std::vector<NodeID> topological_order(dag.numNodes(), 0);
        std::vector<NodeID> part_before(dag.numNodes(), block_0);
        topologicalSort(dag, in_degree, topological_order);
        for (const NodeID& u : topological_order) {
          for (const NodeID& v : scc_to_flow_network[u]) {
            const PartitionID from = hypergraph.partID(v);
            const PartitionID to = best_partition_id[u];
            if (from != to) {
              while(hypergraph.changeNodePart(v, from, to) == false);
              part_before[u] = from;
            }
          }
          // Check cut after assignment of an SCC
          // Should be the same as the starting cut
          // !does not hold in parallel environment
          /*const HyperedgeWeight metric_after = metrics::objective(hypergraph, context.partition.objective);
          if (metric_after != metric_before) {
            LOG << "Assignment of SCC leads to inconsistent hyperedge cut!";
            LOG << V(metric_before) << V(metric_after);
            return false;
          }*/
        }

        // Rollback hypernode assignment
        for (const NodeID& u : dag.nodes()) {
          for (const NodeID& v : scc_to_flow_network[u]) {
            const PartitionID from = hypergraph.partID(v);
            const PartitionID to = part_before[u];
            if (from != to) {
              while(hypergraph.changeNodePart(v, from, to) == false);
            }
          }
        }

        //const HyperedgeWeight metric = metrics::objective(hypergraph, context.partition.objective);
        //metric != metric_before ||
        if (metrics::localImbalance(hypergraph, context) != imbalance_before) {
          LOG << "Restoring original partition failed!";
          //LOG << V(metric_before) << V(metric);
          LOG << V(imbalance_before) << V(metrics::localImbalance(hypergraph, context));
          return false;
        }

        return true;
      } (), "Most balanced minimum cut failed!");

    // Assign most balanced minimum cut
    for (const NodeID& u : dag.nodes()) {
      for (const NodeID& v : scc_to_flow_network[u]) {
        const PartitionID from = hypergraph.partID(v);
        const PartitionID to = best_partition_id[u];
        if (from != to) {
          while(hypergraph.changeNodePart(v, from, to) == false);
        }
      }
    }

    ASSERT(best_imbalance == metrics::localImbalance(hypergraph, context),
           "Best imbalance didn't match with current imbalance"
           << V(best_imbalance) << V(metrics::localImbalance(hypergraph, context)));
  }

 private:
  static constexpr bool debug = false;

  void reset() {
    _visited.reset();
    _graph_to_flow_network.resetUsedEntries();
    _flow_network_to_graph.resetUsedEntries();
    _scc_node_weight.resetUsedEntries();
  }


  /**
   * Executes a BFS starting from the source (sourceSet = true)
   * or sink set (sourceSet = false). Touched nodes by BFS
   * are marked as visited and are not considered during
   * residual graph building step.
   *
   * @t_param sourceSet Indicates, if BFS start from source or sink set
   */
  template <bool sourceSet>
  void markAllReachableNodesAsVisited(Hypergraph& hypergraph, FlowNetwork& flow_network,
                                      const PartitionID block_0, const PartitionID block_1) {
    auto start_set_iterator = sourceSet ? flow_network.sources() : flow_network.sinks();
    for (const NodeID& node : start_set_iterator) {
      _Q.push(node);
      _visited.set(node, true);
    }

    while (!_Q.empty()) {
      const NodeID u = _Q.front();
      _Q.pop();

      if (flow_network.interpreteHypernode(u)) {
        if (!sourceSet) {
          const PartitionID from = hypergraph.partID(u);
          if (from == block_0) {
            while(hypergraph.changeNodePart(u, block_0, block_1) == false);
          }
        }
      } else if (flow_network.interpreteHyperedge(u, sourceSet)) {
        const HyperedgeID he = flow_network.mapToHyperedgeID(u);
        for (const HypernodeID& pin : hypergraph.pins(he)) {
          if (flow_network.containsHypernode(pin)) {
            if (!sourceSet) {
              PartitionID from = hypergraph.partID(pin);
              if (from == block_0) {
                while(hypergraph.changeNodePart(pin, block_0, block_1) == false);
              }
            }
            if (flow_network.isRemovedHypernode(pin)) {
              _visited.set(pin, true);
            }
          }
        }
      }

      for (FlowEdge& e : flow_network.incidentEdges(u)) {
        const FlowEdge& reverse_edge = flow_network.reverseEdge(e);
        const NodeID v = e.target;
        if (!_visited[v]) {
          if ((sourceSet && flow_network.residualCapacity(e)) ||
              (!sourceSet && flow_network.residualCapacity(reverse_edge))) {
            _Q.push(v);
            _visited.set(v, true);
          }
        }
      }
    }
  }

  Graph buildResidualGraph(HyperGraph& hypergraph, FlowNetwork& flow_network) {
    size_t cur_graph_node = 0;
    for (const NodeID& node : flow_network.nodes()) {
      if (!_visited[node]) {
        _graph_to_flow_network.set(cur_graph_node, node);
        _flow_network_to_graph.set(node, cur_graph_node++);
      }
    }

    for (const HypernodeID& hn : flow_network.removedHypernodes()) {
      if (!_visited[hn]) {
        _graph_to_flow_network.set(cur_graph_node, hn);
        _flow_network_to_graph.set(hn, cur_graph_node++);
      }
    }

    std::vector<std::vector<Edge> > adj_list(cur_graph_node, std::vector<Edge>());

    for (const NodeID& node : flow_network.nodes()) {
      if (!_visited[node]) {
        const NodeID source = _flow_network_to_graph.get(node);
        for (FlowEdge& flow_edge : flow_network.incidentEdges(node)) {
          const NodeID target = flow_edge.target;
          if (flow_network.residualCapacity(flow_edge) && !_visited[target]) {
            Edge e;
            e.target_node = _flow_network_to_graph.get(target);
            e.weight = 1.0;
            adj_list[source].push_back(e);
          }
        }
      }
    }

    for (const HypernodeID& hn : flow_network.removedHypernodes()) {
      if (!_visited[hn]) {
        const NodeID hn_node = _flow_network_to_graph.get(hn);
        for (const HyperedgeID& he : hypergraph.incidentEdges(hn)) {
          const NodeID in_he = _flow_network_to_graph.get(flow_network.mapToIncommingHyperedgeID(he));
          const NodeID out_he = _flow_network_to_graph.get(flow_network.mapToOutgoingHyperedgeID(he));
          if (in_he != Network::kInvalidNode) {
            Edge e;
            e.target_node = in_he;
            e.weight = 1.0;
            adj_list[hn_node].push_back(e);
          }
          if (out_he != Network::kInvalidNode) {
            Edge e;
            e.target_node = hn_node;
            e.weight = 1.0;
            adj_list[out_he].push_back(e);
          }
        }
      }
    }

    std::vector<NodeID> adj_array(1, 0);
    std::vector<Edge> edges;
    for (NodeID u = 0; u < cur_graph_node; ++u) {
      for (const Edge& e : adj_list[u]) {
        edges.push_back(e);
      }
      adj_array.push_back(adj_array[u] + adj_list[u].size());
    }

    return Graph(adj_array, edges);
  }

  void findStronglyConnectedComponents(Graph& g) {
    _sccs.compute(g);
  }


  void topologicalSort(const Graph& g,
                       std::vector<size_t> in_degree,
                       std::vector<NodeID>& topological_order) {
    std::vector<NodeID> start_nodes;
    for (const NodeID& u : g.nodes()) {
      if (in_degree[u] == 0) {
        start_nodes.push_back(u);
      }
    }
    utils::Randomize::instance().shuffleVector(start_nodes);
    for (const NodeID& u : start_nodes) {
      _Q.push(u);
    }

    size_t idx = 0;
    while (!_Q.empty()) {
      const NodeID u = _Q.front();
      _Q.pop();
      topological_order[idx++] = u;
      for (const Edge& e : g.incidentEdges(u)) {
        const NodeID v = e.target_node;
        if (u != v) {
          in_degree[v]--;
          if (in_degree[v] == 0) {
            _Q.push(v);
          }
        }
      }
    }

    ASSERT(idx == g.numNodes(), "Topological sort failed!" << V(idx) << V(g.numNodes()));
  }


  template <bool bipartition = true>
  double imbalance(const Context& context, const std::vector<HypernodeWeight>& part_weight) {
    if (bipartition) {
      const HypernodeWeight weight_part0 = part_weight[0];
      const HypernodeWeight weight_part1 = part_weight[1];
      const double imbalance_part0 = (weight_part0 /
                                      static_cast<double>(context.partition.perfect_balance_part_weights[0]));
      const double imbalance_part1 = (weight_part1 /
                                      static_cast<double>(context.partition.perfect_balance_part_weights[1]));
      return std::max(imbalance_part0, imbalance_part1) - 1.0;
    } else {
      double max_balance = (part_weight[0] /
                            static_cast<double>(context.partition.perfect_balance_part_weights[0]));
      for (PartitionID i = 1; i < context.partition.k; ++i) {
        const double balance_i =
          (part_weight[i] /
           static_cast<double>(context.partition.perfect_balance_part_weights[i]));
        max_balance = std::max(max_balance, balance_i);
      }
      return max_balance - 1.0;
    }
  }

  kahypar::ds::FastResetFlagArray<> _visited;
  kahypar::ds::FastResetArray<NodeID> _graph_to_flow_network;
  kahypar::ds::FastResetArray<NodeID> _flow_network_to_graph;
  kahypar::ds::FastResetArray<HypernodeWeight> _scc_node_weight;

  std::queue<NodeID> _Q;  // BFS queue
  StronglyConnectedComponents _sccs;
};
}  // namespace mt_kahypar