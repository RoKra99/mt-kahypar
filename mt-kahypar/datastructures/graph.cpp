/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019, 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "graph.h"


#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_invoke.h>
#include <tbb/enumerable_thread_specific.h>

#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar::ds {

  Graph::Graph(Hypergraph& hypergraph, const LouvainEdgeWeight edge_weight_type) :
          _num_nodes(0),
          _num_arcs(0),
          _total_volume(0),
          _max_degree(0),
          _indices(),
          _arcs(),
          _node_volumes(),
          _tmp_graph_buffer(nullptr) {

    switch( edge_weight_type ) {
      case LouvainEdgeWeight::uniform:
        construct(hypergraph,
                  [&](const HyperedgeWeight edge_weight,
                      const HypernodeID,
                      const HyperedgeID) {
                    return static_cast<ArcWeight>(edge_weight);
                  });
        break;
      case LouvainEdgeWeight::non_uniform:
        construct(hypergraph,
                  [&](const HyperedgeWeight edge_weight,
                      const HypernodeID edge_size,
                      const HyperedgeID) {
                    return static_cast<ArcWeight>(edge_weight) /
                           static_cast<ArcWeight>(edge_size);
                  });
        break;
      case LouvainEdgeWeight::degree:
        construct(hypergraph,
                  [&](const HyperedgeWeight edge_weight,
                      const HypernodeID edge_size,
                      const HyperedgeID node_degree) {
                    return static_cast<ArcWeight>(edge_weight) *
                           (static_cast<ArcWeight>(node_degree) /
                            static_cast<ArcWeight>(edge_size));
                  });
        break;
      case LouvainEdgeWeight::hybrid:
      case LouvainEdgeWeight::UNDEFINED:
      ERROR("No valid louvain edge weight");
    }
  }

  Graph::Graph(Graph&& other) :
    _num_nodes(other._num_nodes),
    _num_arcs(other._num_arcs),
    _total_volume(other._total_volume),
    _max_degree(other._max_degree),
    _indices(std::move(other._indices)),
    _arcs(std::move(other._arcs)),
    _node_volumes(std::move(other._node_volumes)),
    _tmp_graph_buffer(std::move(other._tmp_graph_buffer)) {
    other._num_nodes = 0;
    other._num_arcs = 0;
    other._total_volume = 0;
    other._max_degree = 0;
    other._tmp_graph_buffer = nullptr;
  }

  Graph& Graph::operator= (Graph&& other) {
    _num_nodes = other._num_nodes;
    _num_arcs = other._num_arcs;
    _total_volume = other._total_volume;
    _max_degree = other._max_degree;
    _indices = std::move(other._indices);
    _arcs = std::move(other._arcs);
    _node_volumes = std::move(other._node_volumes);
    _tmp_graph_buffer = std::move(other._tmp_graph_buffer);
    other._num_nodes = 0;
    other._num_arcs = 0;
    other._total_volume = 0;
    other._max_degree = 0;
    other._tmp_graph_buffer = nullptr;
    return *this;
  }

  Graph::~Graph() {
    if ( _tmp_graph_buffer ) {
      delete(_tmp_graph_buffer);
    }
  }


  /*!
 * Contracts the graph based on the community structure passed as argument.
 * In the first step the community ids are compactified (via parallel prefix sum)
 * which also determines the node ids in the coarse graph. Afterwards, we create
 * a temporary graph which contains all arcs that will not form a selfloop in the
 * coarse graph. Finally, the weights of each multiedge in that temporary graph
 * are aggregated and the result is written to the final contracted graph.
 */
  Graph Graph::contract(Clustering& communities) {
    ASSERT(canBeUsed());
    ASSERT(_num_nodes == communities.size());
    ASSERT(_tmp_graph_buffer);
    Graph coarse_graph;
    coarse_graph._total_volume = _total_volume;

    // #################### STAGE 1 ####################
    // Compute node ids of coarse graph with a parallel prefix sum
    utils::Timer::instance().start_timer("compute_cluster_mapping", "Compute Cluster Mapping");
    parallel::scalable_vector<size_t> mapping(_num_nodes, 0UL);
    ds::Array<parallel::IntegralAtomicWrapper<size_t>>& tmp_pos = _tmp_graph_buffer->tmp_pos;
    ds::Array<parallel::IntegralAtomicWrapper<size_t>>& tmp_indices = _tmp_graph_buffer->tmp_indices;
    ds::Array<parallel::AtomicWrapper<ArcWeight>>& coarse_node_volumes = _tmp_graph_buffer->tmp_node_volumes;
    tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
      ASSERT(static_cast<size_t>(communities[u]) < _num_nodes);
      mapping[communities[u]] = 1UL;
      tmp_pos[u] = 0;
      tmp_indices[u] = 0;
      coarse_node_volumes[u].store(0.0);
    });

    // Prefix sum determines node ids in coarse graph
    parallel::TBBPrefixSum<size_t> mapping_prefix_sum(mapping);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, _num_nodes), mapping_prefix_sum);

    // Remap community ids
    coarse_graph._num_nodes = mapping_prefix_sum.total_sum();
    tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
      communities[u] = mapping_prefix_sum[communities[u]];
    });
    utils::Timer::instance().stop_timer("compute_cluster_mapping");

    // #################### STAGE 2 ####################
    // Write all arcs, that will not form a selfloop in the coarse graph, into a tmp
    // adjacence array. For that, we compute a prefix sum over the sum of all arcs
    // in each community (which are no selfloop) and write them in parallel to
    // the tmp adjacence array.
    utils::Timer::instance().start_timer("construct_tmp_adjacent_array", "Construct Tmp Adjacent Array");
    // Compute number of arcs in tmp adjacence array with parallel prefix sum
    ASSERT(coarse_graph._num_nodes <= coarse_node_volumes.size());
    tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
      const NodeID coarse_u = communities[u];
      ASSERT(static_cast<size_t>(coarse_u) < coarse_graph._num_nodes);
      coarse_node_volumes[coarse_u] += nodeVolume(u);
      for ( const Arc& arc : arcsOf(u) ) {
        const NodeID coarse_v = communities[arc.head];
        if ( coarse_u != coarse_v ) {
          ++tmp_indices[coarse_u];
        }
      }
    });

    parallel::TBBPrefixSum<parallel::IntegralAtomicWrapper<size_t>, ds::Array> tmp_indices_prefix_sum(tmp_indices);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, _num_nodes), tmp_indices_prefix_sum);

    // Write all arcs into corresponding tmp adjacence array blocks
    ds::Array<Arc>& tmp_arcs = _tmp_graph_buffer->tmp_arcs;
    ds::Array<size_t>& valid_arcs = _tmp_graph_buffer->valid_arcs;
    tbb::parallel_for(0U, static_cast<NodeID>(_num_nodes), [&](const NodeID u) {
      const NodeID coarse_u = communities[u];
      ASSERT(static_cast<size_t>(coarse_u) < coarse_graph._num_nodes);
      for ( const Arc& arc : arcsOf(u) ) {
        const NodeID coarse_v = communities[arc.head];
        if ( coarse_u != coarse_v ) {
          const size_t tmp_arcs_pos = tmp_indices_prefix_sum[coarse_u] + tmp_pos[coarse_u]++;
          ASSERT(tmp_arcs_pos < tmp_indices_prefix_sum[coarse_u + 1]);
          tmp_arcs[tmp_arcs_pos] = Arc { coarse_v, arc.weight };
          valid_arcs[tmp_arcs_pos] = 1UL;
        }
      }
    });
    utils::Timer::instance().stop_timer("construct_tmp_adjacent_array");

    // #################### STAGE 3 ####################
    // Aggregate weights of arcs that are equal in each community.
    // Therefore, we sort the arcs according to their endpoints
    // and aggregate weight of arcs with equal endpoints.
    utils::Timer::instance().start_timer("contract_arcs", "Contract Arcs");
    tbb::enumerable_thread_specific<size_t> local_max_degree(0);
    tbb::parallel_for(0U, static_cast<NodeID>(coarse_graph._num_nodes), [&](const NodeID u) {
      const size_t tmp_arc_start = tmp_indices_prefix_sum[u];
      const size_t tmp_arc_end = tmp_indices_prefix_sum[u + 1];
      std::sort(tmp_arcs.begin() + tmp_arc_start, tmp_arcs.begin() + tmp_arc_end,
                [&](const Arc& lhs, const Arc& rhs) {
                  return lhs.head < rhs.head;
                });

      size_t arc_rep = tmp_arc_start;
      size_t degree = tmp_arc_start < tmp_arc_end ? 1 : 0;
      for ( size_t pos = tmp_arc_start + 1; pos < tmp_arc_end; ++pos ) {
        if ( tmp_arcs[arc_rep].head == tmp_arcs[pos].head ) {
          tmp_arcs[arc_rep].weight += tmp_arcs[pos].weight;
          valid_arcs[pos] = 0UL;
        } else {
          arc_rep = pos;
          ++degree;
        }
      }
      local_max_degree.local() = std::max(local_max_degree.local(), degree);
    });
    coarse_graph._max_degree = local_max_degree.combine(
            [&](const size_t& lhs, const size_t& rhs) {
              return std::max(lhs, rhs);
            });

    // Write all arcs to coarse graph
    parallel::TBBPrefixSum<size_t, ds::Array> valid_arcs_prefix_sum(valid_arcs);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL,
                                                  tmp_indices_prefix_sum.total_sum()), valid_arcs_prefix_sum);
    coarse_graph._num_arcs = valid_arcs_prefix_sum.total_sum();

    // Move memory down to coarse graph
    coarse_graph._indices = std::move(_indices);
    coarse_graph._arcs = std::move(_arcs);
    coarse_graph._node_volumes = std::move(_node_volumes);

    tbb::parallel_invoke([&] {
      const size_t tmp_num_arcs = tmp_indices_prefix_sum.total_sum();
      tbb::parallel_for(0UL, tmp_num_arcs, [&](const size_t i) {
        if ( valid_arcs_prefix_sum.value(i) ) {
          const size_t pos = valid_arcs_prefix_sum[i];
          ASSERT(pos < coarse_graph._num_arcs);
          coarse_graph._arcs[pos] = std::move(tmp_arcs[i]);
        }
      });
    }, [&] {
      tbb::parallel_for(0U, static_cast<NodeID>(coarse_graph._num_nodes), [&](const NodeID u) {
        const size_t start_index_pos = valid_arcs_prefix_sum[tmp_indices_prefix_sum[u]];
        ASSERT(start_index_pos <= coarse_graph._num_arcs);
        coarse_graph._indices[u] = start_index_pos;
        coarse_graph._node_volumes[u] = coarse_node_volumes[u];
      });
      coarse_graph._indices[coarse_graph._num_nodes] = coarse_graph._num_arcs;
    });
    coarse_graph._tmp_graph_buffer = _tmp_graph_buffer;
    _tmp_graph_buffer = nullptr;
    utils::Timer::instance().stop_timer("contract_arcs");

    return coarse_graph;
  }


  Graph::Graph() :
          _num_nodes(0),
          _num_arcs(0),
          _total_volume(0),
          _max_degree(0),
          _indices(),
          _arcs(),
          _node_volumes(),
          _tmp_graph_buffer(nullptr) {

  }




  /*!
   * Constructs a graph from a given hypergraph.
   */
  template<typename F>
  void Graph::construct(const Hypergraph& hypergraph,
                 const F& edge_weight_func) {
    // Test, if hypergraph is actually a graph
    const bool is_graph = tbb::parallel_reduce(tbb::blocked_range<HyperedgeID>(
            ID(0), hypergraph.initialNumEdges()), true, [&](const tbb::blocked_range<HyperedgeID>& range, bool isGraph) {
      if ( isGraph ) {
        bool tmp_is_graph = isGraph;
        for (HyperedgeID he = range.begin(); he < range.end(); ++he) {
          if ( hypergraph.edgeIsEnabled(he) ) {
            tmp_is_graph &= (hypergraph.edgeSize(he) == 2);
          }
        }
        return tmp_is_graph;
      }
      return false;
    }, [&](const bool lhs, const bool rhs) {
      return lhs && rhs;
    });

    if ( is_graph ) {
      _num_nodes = hypergraph.initialNumNodes();
      _num_arcs = 2 * hypergraph.initialNumEdges();
      constructGraph(hypergraph, edge_weight_func);
    } else {
      _num_nodes = hypergraph.initialNumNodes() + hypergraph.initialNumEdges();
      _num_arcs = 2 * hypergraph.initialNumPins();
      constructBipartiteGraph(hypergraph, edge_weight_func);
    }

    // Compute node volumes and total volume
    utils::Timer::instance().start_timer("compute_node_volumes", "Compute Node Volumes");
    _total_volume = 0.0;
    tbb::enumerable_thread_specific<ArcWeight> local_total_volume(0.0);
    tbb::parallel_for(0U, static_cast<NodeID>(numNodes()), [&](const NodeID u) {
      local_total_volume.local() += computeNodeVolume(u);
    });
    _total_volume = local_total_volume.combine(std::plus<ArcWeight>());
    utils::Timer::instance().stop_timer("compute_node_volumes");
  }

  template<typename F>
  void Graph::constructBipartiteGraph(const Hypergraph& hypergraph,
                               F& edge_weight_func) {
    _indices.resize("Preprocessing", "indices", _num_nodes + 1);
    _arcs.resize("Preprocessing", "arcs", _num_arcs);
    _node_volumes.resize("Preprocessing", "node_volumes", _num_nodes);
    _tmp_graph_buffer = new TmpGraphBuffer(_num_nodes, _num_arcs);

    // Initialize data structure
    utils::Timer::instance().start_timer("compute_node_degrees", "Compute Node Degrees");
    const HypernodeID num_hypernodes = hypergraph.initialNumNodes();
    const HypernodeID num_hyperedges = hypergraph.initialNumEdges();
    tbb::parallel_invoke([&] {
      tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID u) {
        ASSERT(u + 1 < _indices.size());
        _indices[u + 1] = hypergraph.nodeDegree(u);
      });
    }, [&] {
      tbb::parallel_for(num_hypernodes, num_hypernodes + num_hyperedges, [&](const HyperedgeID u) {
        ASSERT(u + 1 < _indices.size());
        const HyperedgeID he = u - num_hypernodes;
        _indices[u + 1] = hypergraph.edgeSize(he);
      });
    });

    parallel::TBBPrefixSum<size_t, ds::Array> indices_prefix_sum(_indices);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, _indices.size()), indices_prefix_sum);
    utils::Timer::instance().stop_timer("compute_node_degrees");

    utils::Timer::instance().start_timer("construct_arcs", "Construct Arcs");
    tbb::enumerable_thread_specific<size_t> local_max_degree(0);
    tbb::parallel_invoke([&] {
      tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID u) {
        ASSERT(u + 1 < _indices.size());
        size_t pos = _indices[u];
        const HypernodeID hn = u;
        const HyperedgeID node_degree = hypergraph.nodeDegree(hn);
        local_max_degree.local() = std::max(
                local_max_degree.local(), static_cast<size_t>(node_degree));
        for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
          const NodeID v = he + num_hypernodes;
          const HyperedgeWeight edge_weight = hypergraph.edgeWeight(he);
          const HypernodeID edge_size = hypergraph.edgeSize(he);
          ASSERT(pos < _indices[u + 1]);
          _arcs[pos++] = Arc(v, edge_weight_func(edge_weight, edge_size, node_degree));
        }
      });
    }, [&] {
      tbb::parallel_for(num_hypernodes, num_hypernodes + num_hyperedges, [&](const HyperedgeID u) {
        ASSERT(u + 1 < _indices.size());
        size_t pos = _indices[u];
        const HyperedgeID he = u - num_hypernodes;
        const HyperedgeWeight edge_weight = hypergraph.edgeWeight(he);
        const HypernodeID edge_size = hypergraph.edgeSize(he);
        local_max_degree.local() = std::max(
                local_max_degree.local(), static_cast<size_t>(edge_size));
        for ( const HypernodeID& pin : hypergraph.pins(he) ) {
          const NodeID v = pin;
          const HyperedgeID node_degree = hypergraph.nodeDegree(pin);
          ASSERT(pos < _indices[u + 1]);
          _arcs[pos++] = Arc(v, edge_weight_func(edge_weight, edge_size, node_degree));
        }
      });
    });
    _max_degree = local_max_degree.combine([&](const size_t& lhs, const size_t& rhs) {
      return std::max(lhs, rhs);
    });
    utils::Timer::instance().stop_timer("construct_arcs");
  }

  template<typename F>
  void Graph::constructGraph(const Hypergraph& hypergraph, const F& edge_weight_func) {
    _indices.resize("Preprocessing", "indices", _num_nodes + 1);
    _arcs.resize("Preprocessing", "arcs", _num_arcs);
    _node_volumes.resize("Preprocessing", "node_volumes", _num_nodes);
    _tmp_graph_buffer = new TmpGraphBuffer(_num_nodes, _num_arcs);

    // Initialize data structure
    utils::Timer::instance().start_timer("compute_node_degrees", "Compute Node Degrees");
    const HypernodeID num_hypernodes = hypergraph.initialNumNodes();
    tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID u) {
      ASSERT(u + 1 < _indices.size());
      _indices[u + 1] = hypergraph.nodeDegree(u);
    });

    parallel::TBBPrefixSum<size_t, ds::Array> indices_prefix_sum(_indices);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, num_hypernodes + 1), indices_prefix_sum);
    utils::Timer::instance().stop_timer("compute_node_degrees");

    utils::Timer::instance().start_timer("construct_arcs", "Construct Arcs");
    tbb::enumerable_thread_specific<size_t> local_max_degree(0);
    tbb::parallel_for(ID(0), num_hypernodes, [&](const HypernodeID u) {
      ASSERT(u + 1 < _indices.size());
      size_t pos = _indices[u];
      const HypernodeID hn = u;
      const HyperedgeID node_degree = hypergraph.nodeDegree(hn);
      local_max_degree.local() = std::max(
              local_max_degree.local(), static_cast<size_t>(node_degree));
      for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
        const HyperedgeWeight edge_weight = hypergraph.edgeWeight(he);
        NodeID v = std::numeric_limits<NodeID>::max();
        for ( const HypernodeID& pin : hypergraph.pins(he) ) {
          if ( pin != hn ) {
            v = pin;
            break;
          }
        }
        ASSERT(v != std::numeric_limits<NodeID>::max());
        ASSERT(pos < _indices[u + 1]);
        _arcs[pos++] = Arc(v, edge_weight_func(edge_weight, ID(2), node_degree));
      }
    });
    _max_degree = local_max_degree.combine([&](const size_t& lhs, const size_t& rhs) {
      return std::max(lhs, rhs);
    });
    utils::Timer::instance().stop_timer("construct_arcs");
  }

  bool Graph::canBeUsed(const bool verbose) const {
    const bool result = _indices.size() >= numNodes() + 1 && _arcs.size() >= numArcs() && _node_volumes.size() >= numNodes();
    if (verbose && !result) {
      LOG << "Some of the graph's members were stolen. For example the contract function does this. "
             "Make sure you're calling functions with a fresh graph or catch this condition and reinitialize."
             "If you do reinitialize, feel free to silence this warning by passing false to the canBeUsed function";
    }
    return result;
  }

} // namespace mt_kahypar::ds