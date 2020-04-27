/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"

#include "kahypar/meta/policy_registry.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/partitioning_output.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/factories.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/partition/multilevel.h"
#include "mt-kahypar/partition/preprocessing/sparsification/degree_zero_hn_remover.h"
#include "mt-kahypar/partition/preprocessing/community_detection/parallel_louvain.h"
#include "mt-kahypar/partition/preprocessing/community_reassignment/community_redistributor.h"
#include "mt-kahypar/utils/stats.h"

namespace mt_kahypar {
namespace partition {
class Partitioner {
 private:
  static constexpr bool debug = false;

 public:
  Partitioner(Context& context) :
    _context(context),
    _degree_zero_hn_remover(context) { }

  Partitioner(const Partitioner&) = delete;
  Partitioner & operator= (const Partitioner &) = delete;

  Partitioner(Partitioner&&) = delete;
  Partitioner & operator= (Partitioner &&) = delete;

  inline PartitionedHypergraph<> partition(Hypergraph& hypergraph);

 private:
  static inline void setupContext(Hypergraph& hypergraph, Context& context);

  static inline void configurePreprocessing(const Hypergraph& hypergraph, Context& context);

  inline void sanitize(Hypergraph& hypergraph);

  inline void preprocess(Hypergraph& hypergraph);

  inline void redistribution(Hypergraph& hypergraph);

  inline void postprocess(PartitionedHypergraph<>& hypergraph);

  Context& _context;
  DegreeZeroHypernodeRemover _degree_zero_hn_remover;
};

inline void Partitioner::setupContext(Hypergraph& hypergraph, Context& context) {
  context.setupPartWeights(hypergraph.totalWeight());
  context.setupContractionLimit(hypergraph.totalWeight());
  context.sanityCheck();
}

inline void Partitioner::configurePreprocessing(const Hypergraph& hypergraph, Context& context) {
  const double density = static_cast<double>(hypergraph.initialNumEdges()) /
                         static_cast<double>(hypergraph.initialNumNodes());
  if (context.preprocessing.community_detection.edge_weight_function == LouvainEdgeWeight::hybrid) {
    if (density < 0.75) {
      context.preprocessing.community_detection.edge_weight_function = LouvainEdgeWeight::degree;
    } else {
      context.preprocessing.community_detection.edge_weight_function = LouvainEdgeWeight::uniform;
    }
  }
}

inline void Partitioner::sanitize(Hypergraph& hypergraph) {

  utils::Timer::instance().start_timer("degree_zero_hypernode_removal", "Degree Zero Hypernode Removal");
  const HypernodeID num_removed_degree_zero_hypernodes =
    _degree_zero_hn_remover.contractDegreeZeroHypernodes(hypergraph);
  utils::Timer::instance().stop_timer("degree_zero_hypernode_removal");

  const HyperedgeID num_removed_single_node_hes = hypergraph.numRemovedHyperedges();
  if (_context.partition.verbose_output &&
      ( num_removed_single_node_hes > 0 || num_removed_degree_zero_hypernodes > 0 )) {
    LOG << "Performed single-node HE removal and degree-zero HN contractions:";
    LOG << "\033[1m\033[31m" << " # removed"
        << num_removed_single_node_hes << "hyperedges with |e| = 1 during hypergraph file parsing"
        << "\033[0m";
    LOG << "\033[1m\033[31m" << " # contracted"
        << num_removed_degree_zero_hypernodes << "hypernodes with d(v) = 0 to supervertices"
        << "\033[0m";
    io::printStripe();
  }
}

inline void Partitioner::preprocess(Hypergraph& hypergraph) {
  for (int node = 0; node < TBBNumaArena::instance().num_used_numa_nodes(); ++node) {
    utils::Stats::instance().add_stat("initial_hns_on_numa_node_" + std::to_string(node),
                                      (int64_t)hypergraph.initialNumNodes(node));
    utils::Stats::instance().add_stat("initial_hes_on_numa_node_" + std::to_string(node),
                                      (int64_t)hypergraph.initialNumEdges(node));
    utils::Stats::instance().add_stat("initial_pins_on_numa_node_" + std::to_string(node),
                                      (int64_t)hypergraph.initialNumPins(node));
  }

  if ( _context.preprocessing.use_community_detection ) {
    io::printTopLevelPreprocessingBanner(_context);

    utils::Timer::instance().start_timer("community_detection", "Community Detection");
    utils::Timer::instance().start_timer("perform_community_detection", "Perform Community Detection");
    ds::Clustering communities(0);
    utils::Timer::instance().start_timer("construct_graph", "Construct Graph");
    Graph graph(hypergraph, _context.preprocessing.community_detection.edge_weight_function);
    utils::Timer::instance().stop_timer("construct_graph");
    communities = ParallelModularityLouvain<Graph>::run(graph, _context,
      _context.shared_memory.num_threads);
    utils::Timer::instance().stop_timer("perform_community_detection");

    // Stream community ids into hypergraph
    utils::Timer::instance().start_timer("stream_community_ids", "Stream Community IDs");
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(0UL, hypergraph.initialNumNodes()),
                      [&](const tbb::blocked_range<HypernodeID>& range) {
          for (HypernodeID id = range.begin(); id < range.end(); ++id) {
            const HypernodeID hn = hypergraph.globalNodeID(id);
            hypergraph.setCommunityID(hn, communities[id]);
          }
        });
    utils::Timer::instance().stop_timer("stream_community_ids");

    // Initialize Communities
    utils::Timer::instance().start_timer("initialize_communities", "Initialize Communities");
    hypergraph.initializeCommunities(TBBNumaArena::GLOBAL_TASK_GROUP);
    utils::Timer::instance().stop_timer("initialize_communities");

    utils::Stats::instance().add_stat("num_communities", hypergraph.numCommunities());
    utils::Timer::instance().stop_timer("community_detection");

    if (_context.partition.verbose_output) {
      io::printCommunityInformation(hypergraph);
      io::printStripe();
    }

    // Redistribute Hypergraph based on communities
    utils::Timer::instance().start_timer("redistribution", "Redistribution");
    redistribution(hypergraph);
    utils::Timer::instance().stop_timer("redistribution");
  } else {
    // Per default all communities are assigned to community 0
    utils::Timer::instance().disable();
    hypergraph.initializeCommunities(TBBNumaArena::GLOBAL_TASK_GROUP);
    parallel::scalable_vector<PartitionID> community_node_mapping(1, 0);
    hypergraph.setCommunityNodeMapping(std::move(community_node_mapping));
    utils::Timer::instance().enable();
  }

  parallel::MemoryPool::instance().release_mem_group("Preprocessing");
}

inline void Partitioner::redistribution(Hypergraph& hypergraph) {
  if (_context.preprocessing.use_community_redistribution &&
      TBBNumaArena::instance().num_used_numa_nodes() > 1) {
    std::unique_ptr<preprocessing::ICommunityAssignment> community_assignment =
      RedistributionFactory::getInstance().createObject(
        _context.preprocessing.community_redistribution.assignment_strategy,
        hypergraph, _context);

    parallel::scalable_vector<PartitionID> community_node_mapping =
      community_assignment->computeAssignment();
    HyperedgeWeight remote_pin_count_before = metrics::remotePinCount(hypergraph);
    hypergraph = preprocessing::CommunityRedistributor::redistribute(
      TBBNumaArena::GLOBAL_TASK_GROUP, hypergraph, community_node_mapping);
    HyperedgeWeight remote_pin_count_after = metrics::remotePinCount(hypergraph);
    utils::Stats::instance().add_stat("remote_pin_count_before", remote_pin_count_before);
    utils::Stats::instance().add_stat("remote_pin_count_after", remote_pin_count_after);
    for (int node = 0; node < TBBNumaArena::instance().num_used_numa_nodes(); ++node) {
      utils::Stats::instance().add_stat("hns_on_numa_node_" + std::to_string(node) + "_after_redistribution",
                                        (int64_t)hypergraph.initialNumNodes(node));
      utils::Stats::instance().add_stat("hes_on_numa_node_" + std::to_string(node) + "_after_redistribution",
                                        (int64_t)hypergraph.initialNumEdges(node));
      utils::Stats::instance().add_stat("pins_on_numa_node_" + std::to_string(node) + "_after_redistribution",
                                        (int64_t)hypergraph.initialNumPins(node));
    }
    if (_context.partition.verbose_output) {
      LOG << "Hypergraph Redistribution Results:";
      LOG << " Remote Pin Count Before Redistribution   =" << remote_pin_count_before;
      LOG << " Remote Pin Count After Redistribution    =" << remote_pin_count_after;
      io::printStripe();
    }

    hypergraph.setCommunityNodeMapping(std::move(community_node_mapping));
  }
}

inline void Partitioner::postprocess(PartitionedHypergraph<>& hypergraph) {
  _degree_zero_hn_remover.restoreDegreeZeroHypernodes(hypergraph);
}

inline PartitionedHypergraph<> Partitioner::partition(Hypergraph& hypergraph) {
  configurePreprocessing(hypergraph, _context);
  setupContext(hypergraph, _context);

  io::printContext(_context);
  io::printMemoryPoolConsumption(_context);
  io::printInputInformation(_context, hypergraph);

  // ################## PREPROCESSING ##################
  utils::Profiler::instance().activate("Preprocessing");
  utils::Timer::instance().start_timer("preprocessing", "Preprocessing");
  preprocess(hypergraph);
  sanitize(hypergraph);
  utils::Timer::instance().stop_timer("preprocessing");
  utils::Profiler::instance().deactivate("Preprocessing");

  // ################## MULTILEVEL ##################
  PartitionedHypergraph<> partitioned_hypergraph = multilevel::partition(
    hypergraph, _context, true, TBBNumaArena::GLOBAL_TASK_GROUP);

  postprocess(partitioned_hypergraph);

  if (_context.partition.verbose_output) {
    io::printHypergraphInfo(partitioned_hypergraph, "Uncoarsened Hypergraph",
      _context.partition.show_memory_consumption);
    io::printStripe();
  }

  return partitioned_hypergraph;
}
}  // namespace partition
}  // namespace mt_kahypar
