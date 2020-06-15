/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#pragma once

#include <atomic>

#include <tbb/enumerable_thread_specific.h>

#include "mt-kahypar/datastructures/sparse_map.h"

#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/floating_point_comparisons.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/clustering.h"
#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/randomize.h"

namespace mt_kahypar {
template<typename G /* Graph */>
class PLM {
 private:
  static constexpr bool advancedGainAdjustment = false;

  using AtomicArcWeight = parallel::AtomicWrapper<ArcWeight>;
  using LargeIncidentClusterWeights = ds::FixedSizeSparseMap<PartitionID, ArcWeight>;
  using CacheEfficientIncidentClusterWeights = ds::FixedSizeSparseMap<PartitionID, ArcWeight>;

 public:
  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  explicit PLM(const Context& context,
               size_t numNodes,
               const bool disable_randomization = false) :
    _context(context),
    _max_degree(numNodes),
    _vertex_degree_sampling_threshold(context.preprocessing.community_detection.vertex_degree_sampling_threshold),
    _cluster_volumes(numNodes),
    _local_small_incident_cluster_weight(0),
    _local_large_incident_cluster_weight([&] {
      return construct_large_incident_cluster_weight_map();
    }),
    _disable_randomization(disable_randomization) { }

  ~PLM() {
    tbb::parallel_invoke([&] {
      parallel::parallel_free_thread_local_internal_data(
        _local_small_incident_cluster_weight, [&](CacheEfficientIncidentClusterWeights& data) {
          data.freeInternalData();
        });
    }, [&] {
      parallel::parallel_free_thread_local_internal_data(
        _local_large_incident_cluster_weight, [&](LargeIncidentClusterWeights& data) {
          data.freeInternalData();
        });
    }, [&] {
      parallel::free(_cluster_volumes);
    });
  }

  bool localMoving(G& graph, ds::Clustering& communities) {
    _max_degree = graph.max_degree();
    _reciprocal_total_volume = 1.0 / graph.totalVolume();
    _vol_multiplier_div_by_node_vol = _reciprocal_total_volume;

    parallel::scalable_vector<NodeID> nodes(graph.numNodes());
    tbb::parallel_for(0U, static_cast<NodeID>(graph.numNodes()), [&](const NodeID u) {
      nodes[u] = u;
      communities[u] = u;
      _cluster_volumes[u].store(graph.nodeVolume(u));
    });

    // local moving
    bool clustering_changed = false;
    size_t number_of_nodes_moved = graph.numNodes();
    for (size_t currentRound = 0;
         number_of_nodes_moved >= _context.preprocessing.community_detection.min_vertex_move_fraction * graph.numNodes()
         && currentRound < _context.preprocessing.community_detection.max_pass_iterations;
         currentRound++) {

      if ( !_disable_randomization ) {
        utils::Timer::instance().start_timer("random_shuffle", "Random Shuffle");
        utils::Randomize::instance().parallelShuffleVector(nodes, 0UL, nodes.size());
        utils::Timer::instance().stop_timer("random_shuffle");
      }

      tbb::enumerable_thread_specific<size_t> local_number_of_nodes_moved(0);
      auto moveNode = [&](const NodeID u) {
          const ArcWeight volU = graph.nodeVolume(u);
          const PartitionID from = communities[u];
          PartitionID best_cluster;

          if ( ratingsFitIntoSmallSparseMap(graph, u) ) {
            best_cluster = computeMaxGainCluster(graph, communities, u, _local_small_incident_cluster_weight.local());
          } else {
            LargeIncidentClusterWeights& large_incident_cluster_weight = _local_large_incident_cluster_weight.local();
            large_incident_cluster_weight.setMaxSize(3UL * std::min(_max_degree, _vertex_degree_sampling_threshold));
            best_cluster = computeMaxGainCluster(graph, communities, u, large_incident_cluster_weight);
          }

          if (best_cluster != from) {
            _cluster_volumes[best_cluster] += volU;
            _cluster_volumes[from] -= volU;
            communities[u] = best_cluster;
            ++local_number_of_nodes_moved.local();
          }
        };

      utils::Timer::instance().start_timer("local_moving_round", "Local Moving Round");

      #ifdef KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS
      std::for_each(nodes.begin(), nodes.end(), moveNode);
      #else
      tbb::parallel_for_each(nodes, moveNode);
      #endif

      number_of_nodes_moved = local_number_of_nodes_moved.combine(std::plus<size_t>());
      clustering_changed |= number_of_nodes_moved > 0;
      utils::Timer::instance().stop_timer("local_moving_round");

      DBG << "Louvain-Pass #" << currentRound << " - Modularity:" << metrics::modularity(graph, communities);
    }
    return clustering_changed;
  }

 private:
  FRIEND_TEST(ALouvain, ComputesMaxGainMove1);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove2);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove3);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove4);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove5);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove6);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove7);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove8);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove9);
  FRIEND_TEST(ALouvain, ComputesMaxGainMove10);

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool ratingsFitIntoSmallSparseMap(const G& graph,
                                                                    const HypernodeID u)  {
    static constexpr size_t cache_efficient_map_size = CacheEfficientIncidentClusterWeights::MAP_SIZE / 3UL;
    return std::min(_vertex_degree_sampling_threshold, _max_degree) > cache_efficient_map_size &&
           graph.degree(u) <= cache_efficient_map_size;
  }

  LargeIncidentClusterWeights construct_large_incident_cluster_weight_map() {
    return LargeIncidentClusterWeights(3UL * std::min(_max_degree, _vertex_degree_sampling_threshold), 0);
  }

  // ! Only for testing
  void initializeClusterVolumes(Graph& graph, ds::Clustering& communities) {
    _reciprocal_total_volume = 1.0 / graph.totalVolume();
    _vol_multiplier_div_by_node_vol =  _reciprocal_total_volume;
    tbb::parallel_for(0U, static_cast<NodeID>(graph.numNodes()), [&](const NodeID u) {
      const PartitionID community_id = communities[u];
      _cluster_volumes[community_id] += graph.nodeVolume(u);
    });
  }

  template<typename Map>
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE PartitionID computeMaxGainCluster(const G& graph,
                                                                    ds::Clustering& communities,
                                                                    const NodeID u,
                                                                    Map& incident_cluster_weights) {
    PartitionID from = communities[u];
    PartitionID bestCluster = communities[u];

    for (const Arc& arc : graph.arcsOf(u, _vertex_degree_sampling_threshold)) {
      incident_cluster_weights[communities[arc.head]] += arc.weight;
    }

    const ArcWeight volume_from = _cluster_volumes[from];
    const ArcWeight volU = graph.nodeVolume(u);
    const ArcWeight weight_from = incident_cluster_weights[from];

    const double volMultiplier = _vol_multiplier_div_by_node_vol * volU;
    double bestGain = weight_from - volMultiplier * (volume_from - volU);
    for (const auto& clusterWeight : incident_cluster_weights) {
      PartitionID to = clusterWeight.key;
      // if from == to, we would have to remove volU from volume_to as well.
      // just skip it. it has (adjusted) gain zero.
      if (from == to) {
        continue;
      }

      const ArcWeight volume_to = _cluster_volumes[to],
        weight_to = incident_cluster_weights.get(to);

      double gain = modularityGain(weight_to, volume_to, volMultiplier);

      if (gain > bestGain) {
        bestCluster = to;
        bestGain = gain;
      }
    }

    HEAVY_PREPROCESSING_ASSERT(verifyGain(graph, communities, u, bestCluster, bestGain, incident_cluster_weights));

    incident_cluster_weights.clear();

    return bestCluster;
  }


  inline double modularityGain(const ArcWeight weight_to,
                               const ArcWeight volume_to,
                               const double multiplier) {
    return weight_to - multiplier * volume_to;
    // missing term is - weight_from + multiplier * (volume_from - volume_node)
  }

  inline long double adjustAdvancedModGain(double gain,
                                           const ArcWeight weight_from,
                                           const ArcWeight volume_from,
                                           const ArcWeight volume_node) {
    return 2.0L * _reciprocal_total_volume *
      (gain - weight_from + _reciprocal_total_volume *
        volume_node * (volume_from - volume_node));
  }

  template<typename Map>
  bool verifyGain(const G& graph,
                  ds::Clustering& communities,
                  const NodeID u,
                  const PartitionID to,
                  double gain,
                  const Map& icw) {
    const PartitionID from = communities[u];

    long double adjustedGain = adjustAdvancedModGain(gain, icw.get(from), _cluster_volumes[from], graph.nodeVolume(u));
    const double volMultiplier = _vol_multiplier_div_by_node_vol * graph.nodeVolume(u);
    long double adjustedGainRecomputed = adjustAdvancedModGain(modularityGain(icw.get(to), _cluster_volumes[to], volMultiplier),
                                                               icw.get(from), _cluster_volumes[from], graph.nodeVolume(u));
    unused(adjustedGainRecomputed);

    if (from == to) {
      adjustedGainRecomputed = 0.0L;
      adjustedGain = 0.0L;
    }

    ASSERT(adjustedGain == adjustedGainRecomputed);

    long double dTotalVolumeSquared = static_cast<long double>(graph.totalVolume()) * static_cast<long double>(graph.totalVolume());

    auto accBeforeMove = intraClusterWeightsAndSumOfSquaredClusterVolumes(graph, communities);
    long double coverageBeforeMove = static_cast<long double>(accBeforeMove.first) / graph.totalVolume();
    long double expectedCoverageBeforeMove = accBeforeMove.second / dTotalVolumeSquared;
    long double modBeforeMove = coverageBeforeMove - expectedCoverageBeforeMove;

    // apply move
    communities[u] = to;
    _cluster_volumes[to] += graph.nodeVolume(u);
    _cluster_volumes[from] -= graph.nodeVolume(u);

    auto accAfterMove = intraClusterWeightsAndSumOfSquaredClusterVolumes(graph, communities);
    long double coverageAfterMove = static_cast<long double>(accAfterMove.first) / graph.totalVolume();
    long double expectedCoverageAfterMove = accAfterMove.second / dTotalVolumeSquared;
    long double modAfterMove = coverageAfterMove - expectedCoverageAfterMove;

    const bool result = math::are_almost_equal_ld(modBeforeMove + adjustedGain, modAfterMove, 1e-8);
    ASSERT(result,
           V(modBeforeMove + adjustedGain) << V(modAfterMove) << V(gain) << V(adjustedGain)
           << V(coverageBeforeMove) << V(expectedCoverageBeforeMove) << V(modBeforeMove)
           << V(coverageAfterMove) << V(expectedCoverageAfterMove) << V(modAfterMove));

    // revert move
    communities[u] = from;
    _cluster_volumes[to] -= graph.nodeVolume(u);
    _cluster_volumes[from] += graph.nodeVolume(u);

    return result;
  }

  static std::pair<ArcWeight, ArcWeight> intraClusterWeightsAndSumOfSquaredClusterVolumes(const G& graph, const ds::Clustering& communities) {
    ArcWeight intraClusterWeights = 0;
    ArcWeight sumOfSquaredClusterVolumes = 0;
    std::vector<ArcWeight> _cluster_volumes(graph.numNodes(), 0);

    for (NodeID u : graph.nodes()) {
      ArcWeight arcVol = 0;
      for (const Arc& arc : graph.arcsOf(u)) {
        if (communities[u] == communities[arc.head])
          intraClusterWeights += arc.weight;
        arcVol += arc.weight;
      }

      ArcWeight selfLoopWeight = graph.nodeVolume(u) - arcVol;          // already accounted for as twice!
      ASSERT(selfLoopWeight >= 0.0);
      intraClusterWeights += selfLoopWeight;
      _cluster_volumes[communities[u]] += graph.nodeVolume(u);
    }

    for (NodeID cluster : graph.nodes()) {      // unused cluster IDs have volume 0
      sumOfSquaredClusterVolumes += _cluster_volumes[cluster] * _cluster_volumes[cluster];
    }
    return std::make_pair(intraClusterWeights, sumOfSquaredClusterVolumes);
  }

  const Context& _context;
  size_t _max_degree;
  const size_t _vertex_degree_sampling_threshold;
  double _reciprocal_total_volume = 0.0;
  double _vol_multiplier_div_by_node_vol = 0.0;
  parallel::scalable_vector<AtomicArcWeight> _cluster_volumes;
  tbb::enumerable_thread_specific<CacheEfficientIncidentClusterWeights> _local_small_incident_cluster_weight;
  tbb::enumerable_thread_specific<LargeIncidentClusterWeights> _local_large_incident_cluster_weight;
  const bool _disable_randomization;
};
}
