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

#include "tbb/parallel_for.h"
#include "tbb/task_group.h"

#include "mt-kahypar/definitions.h"

namespace mt_kahypar {
namespace preprocessing {
template <typename TypeTraits>
class CommunityRedistributorT {
 private:
  using HyperGraph = typename TypeTraits::HyperGraph;
  using StreamingHyperGraph = typename TypeTraits::StreamingHyperGraph;
  using TBB = typename TypeTraits::TBB;
  using HwTopology = typename TypeTraits::HwTopology;

  static constexpr bool enable_heavy_assert = false;

 public:
  CommunityRedistributorT(const CommunityRedistributorT&) = delete;
  CommunityRedistributorT & operator= (const CommunityRedistributorT &) = delete;
  CommunityRedistributorT(CommunityRedistributorT&&) = delete;
  CommunityRedistributorT & operator= (CommunityRedistributorT &&) = delete;

  static HyperGraph redistribute(HyperGraph& hg,
                                 const PartitionID k,
                                 const std::vector<PartitionID>& community_assignment) {
    int used_numa_nodes = TBB::instance().num_used_numa_nodes();

    // Compute Node Mapping
    utils::Timer::instance().start_timer("compute_node_mapping", "Compute Node Mapping");
    std::vector<HypernodeID> node_mapping(hg.initialNumNodes(), -1);
    tbb::parallel_for(0UL, hg.initialNumNodes(), [&](const HypernodeID& hn) {
          node_mapping[hn] = community_assignment[hg.communityID(hg.globalNodeID(hn))];
        });
    utils::Timer::instance().stop_timer("compute_node_mapping");

    // Compute Hyperedge Mapping
    utils::Timer::instance().start_timer("compute_hyperedge_mapping", "Compute Hyperedge Mapping");
    std::vector<parallel::scalable_vector<PartitionID> > hyperedge_mapping(used_numa_nodes);
    for (int node = 0; node < used_numa_nodes; ++node) {
      TBB::instance().numa_task_arena(node).execute([&, node] {
            TBB::instance().numa_task_group(TBB::GLOBAL_TASK_GROUP, node).run([&, node] {
              hyperedge_mapping[node].assign(hg.initialNumEdges(node), -1);
              tbb::parallel_for(0UL, hg.initialNumEdges(node), [&](const HyperedgeID& local_he) {
                const HyperedgeID he = StreamingHyperGraph::get_global_edge_id(node, local_he);

                // Compute for each hyperedge the pin count on each node for the new assignment
                parallel::scalable_vector<size_t> pin_count(used_numa_nodes, 0);
                for (const HypernodeID& pin : hg.pins(he)) {
                  ASSERT(hg.communityID(pin) < (PartitionID)community_assignment.size());
                  ++pin_count[community_assignment[hg.communityID(pin)]];
                }

                // Compute assignment based maximum pin count
                size_t max_pin_count = 0;
                PartitionID assigned_node = -1;
                for (PartitionID current_node = 0; current_node < used_numa_nodes; ++current_node) {
                  if (pin_count[current_node] >= max_pin_count) {
                    max_pin_count = pin_count[current_node];
                    assigned_node = current_node;
                  }
                }
                ASSERT(assigned_node != -1);
                hyperedge_mapping[node][local_he] = assigned_node;
              });
            });
          });
    }
    TBB::instance().wait(TBB::GLOBAL_TASK_GROUP);
    utils::Timer::instance().stop_timer("compute_hyperedge_mapping");

    // Reset Pins to original node ids
    utils::Timer::instance().start_timer("reset_pins_to_original_ids", "Reset Pins to original IDs");
    hg.resetPinsToOriginalNodeIds(TBB::GLOBAL_TASK_GROUP);
    utils::Timer::instance().stop_timer("reset_pins_to_original_ids");

    utils::Timer::instance().start_timer("stream_hyperedges", "Stream Hyperedges");
    // Initialize Streaming Hypergraphs
    std::vector<StreamingHyperGraph> numa_hypergraphs;
    TBB::instance().execute_sequential_on_all_numa_nodes(TBB::GLOBAL_TASK_GROUP, [&](const int node) {
          numa_hypergraphs.emplace_back(node, k, TBB::instance().numa_task_arena(node));
        });

    // Stream hyperedges into hypergraphs
    for (int node = 0; node < used_numa_nodes; ++node) {
      for (int streaming_node = 0; streaming_node < used_numa_nodes; ++streaming_node) {
        TBB::instance().numa_task_arena(streaming_node).execute([&, node, streaming_node] {
              TBB::instance().numa_task_group(TBB::GLOBAL_TASK_GROUP, streaming_node).run([&, node, streaming_node] {
                tbb::parallel_for(0UL, hg.initialNumEdges(node), [&, node, streaming_node](const HyperedgeID& local_he) {
                  ASSERT(streaming_node == HwTopology::instance().numa_node_of_cpu(sched_getcpu()));
                  if (hyperedge_mapping[node][local_he] == streaming_node) {
                    const HyperedgeID he = StreamingHyperGraph::get_global_edge_id(node, local_he);
                    parallel::scalable_vector<HypernodeID> hyperedge;
                    for (const HypernodeID& pin : hg.pins(he)) {
                      hyperedge.emplace_back(pin);
                    }
                    numa_hypergraphs[streaming_node].streamHyperedge(hyperedge, hg.originalEdgeID(he), hg.edgeWeight(he));
                  }
                });
              });
            });
      }
      TBB::instance().wait(TBB::GLOBAL_TASK_GROUP);
    }

    // Initialize hyperedges in numa hypergraphs
    for (int node = 0; node < used_numa_nodes; ++node) {
      TBB::instance().numa_task_arena(node).execute([&] {
            TBB::instance().numa_task_group(TBB::GLOBAL_TASK_GROUP, node).run([&, node] {
              numa_hypergraphs[node].initializeHyperedges(hg.initialNumNodes());
            });
          });
    }
    TBB::instance().wait(TBB::GLOBAL_TASK_GROUP);
    utils::Timer::instance().stop_timer("stream_hyperedges");

    // Initialize hypergraph
    utils::Timer::instance().start_timer("initialize hypergraph", "Initialize Hypergraph");
    HyperGraph hypergraph(hg.initialNumNodes(), std::move(numa_hypergraphs),
      std::move(node_mapping), k, TBB::GLOBAL_TASK_GROUP);
    ASSERT(hypergraph.initialNumNodes() == hg.initialNumNodes());
    ASSERT(hypergraph.initialNumEdges() == hg.initialNumEdges());
    ASSERT(hypergraph.initialNumPins() == hg.initialNumPins());
    utils::Timer::instance().stop_timer("initialize hypergraph");

    // Initialize Communities
    utils::Timer::instance().start_timer("initialize_communities", "Initialize Communities");
    tbb::parallel_for(0UL, hypergraph.initialNumNodes(), [&](const HypernodeID& hn) {
          HypernodeID old_global_id = hg.globalNodeID(hn);
          HypernodeID new_global_id = hypergraph.globalNodeID(hn);
          hypergraph.setCommunityID(new_global_id, hg.communityID(old_global_id));
        });
    hypergraph.initializeCommunities();
    utils::Timer::instance().stop_timer("initialize_communities");

    HEAVY_PREPROCESSING_ASSERT([&] {
          for (const HypernodeID& hn : hypergraph.nodes()) {
            int node = StreamingHyperGraph::get_numa_node_of_vertex(hn);
            PartitionID community_id = hypergraph.communityID(hn);
            if (community_assignment[community_id] != node) {
              LOG << "Hypernode" << hn << "should be on numa node" << community_assignment[community_id]
                  << "but is on node" << node;
              return false;
            }
          }
          return true;
        } (), "There are verticies assigned to wrong numa node");

    return hypergraph;
  }

 protected:
  CommunityRedistributorT() = default;
};

using CommunityRedistributor = CommunityRedistributorT<GlobalTypeTraits>;
}  // namespace preprocessing
}  // namespace mt_kahypar
