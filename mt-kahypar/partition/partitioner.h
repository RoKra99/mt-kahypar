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
#include "mt-kahypar/partition/preprocessing/sparsification/large_he_remover.h"
#include "mt-kahypar/partition/preprocessing/community_detection/parallel_louvain.h"
#include "mt-kahypar/partition/preprocessing/community_detection/kahypar_louvain.h"
#include "mt-kahypar/utils/stats.h"

namespace mt_kahypar {
namespace partition {
class Partitioner {
 private:
  static constexpr bool debug = false;

 public:
  Partitioner(Context& context) :
    _context(context),
    _degree_zero_hn_remover(context),
    _large_he_remover(context) { }

  Partitioner(const Partitioner&) = delete;
  Partitioner & operator= (const Partitioner &) = delete;

  Partitioner(Partitioner&&) = delete;
  Partitioner & operator= (Partitioner &&) = delete;

  inline PartitionedHypergraph partition(Hypergraph& hypergraph);

 private:
  static inline void setupContext(Hypergraph& hypergraph, Context& context);

  static inline void configurePreprocessing(const Hypergraph& hypergraph, Context& context);

  inline void sanitize(Hypergraph& hypergraph);

  inline void preprocess(Hypergraph& hypergraph);

  inline void postprocess(PartitionedHypergraph& hypergraph);

  inline PartitionedHypergraph partitionVCycle(Hypergraph& hypergraph, PartitionedHypergraph&& partitioned_hypergraph);

  Context& _context;
  DegreeZeroHypernodeRemover _degree_zero_hn_remover;
  LargeHyperedgeRemover _large_he_remover;
};

inline void Partitioner::setupContext(Hypergraph& hypergraph, Context& context) {
  context.partition.large_hyperedge_size_threshold = std::max(hypergraph.initialNumNodes() *
    context.partition.large_hyperedge_size_threshold_factor, 100.0);
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

  utils::Timer::instance().start_timer("large_hyperedge_removal", "Large Hyperedge Removal");
  const HypernodeID num_removed_large_hyperedges =
    _large_he_remover.removeLargeHyperedges(hypergraph);
  utils::Timer::instance().stop_timer("large_hyperedge_removal");

  const HyperedgeID num_removed_single_node_hes = hypergraph.numRemovedHyperedges();
  if (_context.partition.verbose_output &&
      ( num_removed_single_node_hes > 0 ||
        num_removed_degree_zero_hypernodes > 0 ||
        num_removed_large_hyperedges > 0 )) {
    LOG << "Performed single-node/large HE removal and degree-zero HN contractions:";
    LOG << "\033[1m\033[31m" << " # removed"
        << num_removed_single_node_hes << "single-pin hyperedges during hypergraph file parsing"
        << "\033[0m";
    LOG << "\033[1m\033[31m" << " # removed"
        << num_removed_large_hyperedges << "large hyperedges with |e| >" << _large_he_remover.largeHyperedgeThreshold() << "\033[0m";
    LOG << "\033[1m\033[31m" << " # contracted"
        << num_removed_degree_zero_hypernodes << "hypernodes with d(v) = 0 to supervertices"
        << "\033[0m";
    io::printStripe();
  }
}

inline void Partitioner::preprocess(Hypergraph& hypergraph) {
  if ( _context.preprocessing.use_community_detection ) {
    io::printTopLevelPreprocessingBanner(_context);

    utils::Timer::instance().start_timer("community_detection", "Community Detection");
    utils::Timer::instance().start_timer("perform_community_detection", "Perform Community Detection");
    ds::Clustering communities(0);
    if ( _context.preprocessing.use_kahypar_community_detection ) {
      communities = KaHyParLouvain::run(hypergraph, _context);
    } else {
      utils::Timer::instance().start_timer("construct_graph", "Construct Graph");
      Graph graph(hypergraph, _context.preprocessing.community_detection.edge_weight_function);
      utils::Timer::instance().stop_timer("construct_graph");
      communities = ParallelModularityLouvain<Graph>::run(graph, _context,
        _context.shared_memory.num_threads);
    }
    utils::Timer::instance().stop_timer("perform_community_detection");

    // Stream community ids into hypergraph
    utils::Timer::instance().start_timer("stream_community_ids", "Stream Community IDs");
    tbb::parallel_for(tbb::blocked_range<HypernodeID>(0UL, hypergraph.initialNumNodes()),
                      [&](const tbb::blocked_range<HypernodeID>& range) {
          for (HypernodeID hn = range.begin(); hn < range.end(); ++hn) {
            hypergraph.setCommunityID(hn, communities[hn]);
          }
        });
    utils::Timer::instance().stop_timer("stream_community_ids");

    // Initialize Communities
    utils::Timer::instance().start_timer("initialize_communities", "Initialize Communities");
    hypergraph.initializeCommunities();
    utils::Timer::instance().stop_timer("initialize_communities");

    utils::Stats::instance().add_stat("num_communities", hypergraph.numCommunities());
    utils::Timer::instance().stop_timer("community_detection");

    if (_context.partition.verbose_output) {
      io::printCommunityInformation(hypergraph);
      io::printStripe();
    }
  } else {
    // Per default all communities are assigned to community 0
    utils::Timer::instance().disable();
    hypergraph.initializeCommunities();
    utils::Timer::instance().enable();
  }

  parallel::MemoryPool::instance().release_mem_group("Preprocessing");
}

inline void Partitioner::postprocess(PartitionedHypergraph& hypergraph) {
  _large_he_remover.restoreLargeHyperedges(hypergraph);
  _degree_zero_hn_remover.restoreDegreeZeroHypernodes(hypergraph);
}

inline PartitionedHypergraph Partitioner::partitionVCycle(Hypergraph& hypergraph,
                                                          PartitionedHypergraph&& partitioned_hypergraph) {
  ASSERT(_context.partition.num_vcycles > 0);
  parallel::scalable_vector<PartitionID> part_ids(hypergraph.initialNumNodes(), kInvalidPartition);

  for ( size_t i = 0; i < _context.partition.num_vcycles; ++i ) {
    // Reset memory pool
    parallel::MemoryPool::instance().reset();
    parallel::MemoryPool::instance().release_mem_group("Preprocessing");

    // Store partition and assign it as community ids in order to
    // restrict contractions in v-cycle to partition ids
    hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
      part_ids[hn] = partitioned_hypergraph.partID(hn);
    });
    hypergraph.setCommunityIDs(part_ids);

    // V-Cycle Multilevel Partitioning
    io::printVCycleBanner(_context, i + 1);
    partitioned_hypergraph = multilevel::partition(
      hypergraph, _context, true, TBBNumaArena::GLOBAL_TASK_GROUP, true /* vcycle */);
  }

  return std::move(partitioned_hypergraph);
}

inline PartitionedHypergraph Partitioner::partition(Hypergraph& hypergraph) {
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
  PartitionedHypergraph partitioned_hypergraph = multilevel::partition(
    hypergraph, _context, true, TBBNumaArena::GLOBAL_TASK_GROUP);

  // ################## V-Cycle s##################
  if ( _context.partition.num_vcycles > 0 ) {
    partitioned_hypergraph = partitionVCycle(
      hypergraph, std::move(partitioned_hypergraph));
  }

  // ################## POSTPROCESSING ##################
  utils::Timer::instance().start_timer("postprocessing", "Postprocessing");
  postprocess(partitioned_hypergraph);
  utils::Timer::instance().stop_timer("postprocessing");

  if (_context.partition.verbose_output) {
    io::printHypergraphInfo(partitioned_hypergraph, "Uncoarsened Hypergraph",
      _context.partition.show_memory_consumption);
    io::printStripe();
  }

  return partitioned_hypergraph;
}
}  // namespace partition
}  // namespace mt_kahypar
