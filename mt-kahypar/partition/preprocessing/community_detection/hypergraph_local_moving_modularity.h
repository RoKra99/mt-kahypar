#pragma once

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/community_hypergraph.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/exponentiations.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/datastructures/sparse_map.h"

#include "gtest/gtest_prod.h"


namespace mt_kahypar::community_detection {

class HypergraphLocalMovingModularity {
private:
    using AtomicHyperedgeWeight = parallel::AtomicWrapper<HyperedgeWeight>;
    using CommunityVolumes = parallel::scalable_vector<AtomicHyperedgeWeight>;
    using CommunityVolumeIterator = typename CommunityVolumes::const_iterator;
    using LargeTmpRatingMap = ds::SparseMap<PartitionID, HyperedgeWeight>;
    using CacheEfficientRatingMap = ds::FixedSizeSparseMap<PartitionID, HyperedgeWeight>;
    using ThreadLocalCacheEfficientRatingMap = tbb::enumerable_thread_specific<CacheEfficientRatingMap>;
    using ThreadLocalLargeTmpRatingMap = tbb::enumerable_thread_specific<LargeTmpRatingMap>;

public:
    static constexpr bool debug = false;
    static constexpr bool enable_heavy_assert = false;

    HypergraphLocalMovingModularity(ds::CommunityHypergraph& hypergraph, const Context& context, const bool deactivate_random = false) :
        positive_edgeContribution_count(0),
        //exp_makes_it_bad(0UL),
        //overall_checks(0UL),
        //pruned_by_old(0UL),
        edge_contribution_time(0.0L),
        exp_edge_contribution_time(0.0L),
        _community_neighbour_sampling_threshold(context.preprocessing.community_detection.community_neighbour_sampling_threshold),
        _hyperedge_size_caching_threshold(context.preprocessing.community_detection.hyperedge_size_caching_threshold),
        _vertex_degree_sampling_threshold(context.coarsening.vertex_degree_sampling_threshold),
        _context(context),
        _reciprocal_vol_total(1.0L / hypergraph.totalVolume()),
        _small_edge_contribution_map(0),
        _large_edge_contribution_map([&] {
        return construct_large_edge_contribution_map(hypergraph.initialNumNodes());
            }),
        _community_neighbours_of_edge(0),
                _powers_of_source_community(hypergraph.maxEdgeSize() + 1, 0.L),
                _community_volumes(hypergraph.initialNumNodes()),
                _deactivate_random(deactivate_random) {}


            ~HypergraphLocalMovingModularity() = default;

            // ! calculates the best modularity move for the given node
            template<typename Map>
            KAHYPAR_ATTRIBUTE_ALWAYS_INLINE PartitionID calculateBestMove(ds::CommunityHypergraph& chg, parallel::scalable_vector<HypernodeID>& communities, const HypernodeID v, Map& community_edge_contribution) {
                ASSERT(community_edge_contribution.size() == 0);
                //utils::Timer::instance().start_timer("calculate_best_move", "Calculate best move");
                auto t = tbb::tick_count::now();
                const PartitionID comm_v = communities[v];
                // sum of all edgeweights incident to v
                HyperedgeWeight sum_of_edgeweights = 0;
                //utils::Timer::instance().start_timer("edge_contribution", "EdgeContribution");
                for (const HyperedgeID& he : chg.incidentEdges(v)) {
                    const HyperedgeWeight edge_weight = chg.edgeWeight(he);
                    sum_of_edgeweights += edge_weight;
                    // cuts for this hyperedge are cached
                    if (chg.edgeSize(he) > _hyperedge_size_caching_threshold) {
                        for (const PartitionID community : chg.singleCuts(he)) {
                            //ASSERT(static_cast<HypernodeID>(community) < chg.initialNumNodes());
                            if (static_cast<HypernodeID>(community) < chg.initialNumNodes() && community >= 0) {
                                community_edge_contribution[community] -= edge_weight;
                            }
                        }

                        for (const auto& e : chg.multiCuts(he)) {
                            const PartitionID community = e.first;
                            const size_t count = e.second;
                            if (static_cast<HypernodeID>(community) < chg.initialNumNodes() && (community != comm_v || count == 1)) {
                                ASSERT(count > 0);
                                ASSERT(static_cast<HypernodeID>(community) < chg.initialNumNodes() && community >= 0);
                                community_edge_contribution[community] -= edge_weight;
                            }
                        }
                    } else { // hyperedge is not cached
                        CacheEfficientRatingMap& community_neighbours_of_edge = _community_neighbours_of_edge.local();
                        ASSERT(community_neighbours_of_edge.size() == 0);
                        for (const HypernodeID& hn : chg.pins(he)) {
                            const PartitionID comm_hn = communities[hn];
                            if (hn != v) {
                                community_neighbours_of_edge[comm_hn] += 1U;
                            }
                        }

                        if (!community_neighbours_of_edge.contains(comm_v)) {
                            community_edge_contribution[comm_v] -= edge_weight;
                        }

                        for (const auto& e : community_neighbours_of_edge) {
                            if (e.key != comm_v) {
                                community_edge_contribution[e.key] -= edge_weight;
                            }
                        }
                        community_neighbours_of_edge.clear();
                    }
                }
                const HyperedgeWeight edge_contribution_c = -community_edge_contribution[comm_v];
                //utils::Timer::instance().stop_timer("edge_contribution");
                edge_contribution_time += (tbb::tick_count::now() - t).seconds();
                t = tbb::tick_count::now();
                //utils::Timer::instance().start_timer("exp_edge_contribution", "ExpectedEdgeContribution");

                // ------------------------- Sampling --------------------------------------
                const auto end = community_edge_contribution.end();
                // const auto end = std::min(community_edge_contribution.begin() + _community_neighbour_sampling_threshold, community_edge_contribution.end());
                // if (end != community_edge_contribution.end()) {
                //     std::nth_element(community_edge_contribution.begin(), end, community_edge_contribution.end(), [&](const auto a, const auto b) {
                //         return a.value < b.value;
                //         });
                // }
                // -------------------------------------------------------------------------

                const HyperedgeWeight vol_v = chg.nodeVolume(v);
                const HyperedgeWeight vol_c = _community_volumes[comm_v];
                const HyperedgeWeight vol_c_minus_vol_v = vol_c - vol_v;

                const Volume source_fraction_minus = 1.0L - static_cast<Volume>(vol_c_minus_vol_v) * _reciprocal_vol_total;
                const Volume source_fraction = 1.0L - static_cast<Volume>(vol_c) * _reciprocal_vol_total;
                size_t biggest_d_yet = 1;
                Volume power_d_fraction_minus = source_fraction_minus;
                Volume power_d_fraction = source_fraction;
                bool calculated_c = false;

                PartitionID best_community = comm_v;
                Volume best_delta = 0.0L;
                const HyperedgeWeight sum_of_edgeweights_minus_edgecontribution_c = sum_of_edgeweights - edge_contribution_c;

                // expected edgecontribution starts here
                // for (const auto& e : community_edge_contribution) {
                for (auto it = community_edge_contribution.begin(); it != end; ++it) {
                    //++overall_checks;
                    const auto& e = *it;
                    const PartitionID community = e.key;

                    if (community == comm_v) {
                        continue;
                    }
                    const HyperedgeWeight vol_destination_minus = _community_volumes[community];
                    const HyperedgeWeight vol_destination = vol_destination_minus + vol_v;
                    const HyperedgeWeight destination_edge_contribution = e.value + sum_of_edgeweights_minus_edgecontribution_c;



                    // delta will not be < 0
                    if ((destination_edge_contribution >= 0 || best_delta < destination_edge_contribution)
                        && vol_c_minus_vol_v <= vol_destination_minus) {
                        //++pruned_by_old;
                        continue;
                    }

                    // #############################################################
                    ranking_after_km1.local().push_back(std::pair<PartitionID, Volume>(e.key, e.value));
                    // #############################################################

                    const Volume destination_fraction = 1.0L - static_cast<Volume>(vol_destination) * _reciprocal_vol_total;
                    const Volume destination_fraction_minus = 1.0L - static_cast<Volume>(vol_destination_minus) * _reciprocal_vol_total;
                    parallel::scalable_vector<Volume>& powers_of_source_community = _powers_of_source_community.local();

                    // precalculate the powers for the source community only once
                    // and only if not every possible move is pruned beforehand
                    if (!calculated_c) {
                        for (const size_t d : chg.edgeSizes()) {
                            const size_t remaining_d = d - biggest_d_yet;
                            power_d_fraction_minus *= math::fast_power(source_fraction_minus, remaining_d);
                            power_d_fraction *= math::fast_power(source_fraction, remaining_d);
                            powers_of_source_community[d] = power_d_fraction_minus - power_d_fraction;
                            biggest_d_yet = d;
                        }
                        calculated_c = true;
                    }

                    //Volume exp_edge_contribution = 0.0L;
                    Volume delta = destination_edge_contribution;
                    // if this is equal the expected_edge_contribution will be 0
                    if (vol_c_minus_vol_v != vol_destination_minus) {
                        biggest_d_yet = 1;
                        power_d_fraction_minus = destination_fraction_minus;
                        power_d_fraction = destination_fraction;
                        //actual calculation of the expected edge contribution for the given community
                        for (const size_t d : chg.edgeSizes()) {
                            const size_t remaining_d = d - biggest_d_yet;
                            power_d_fraction_minus *= math::fast_power(destination_fraction_minus, remaining_d);
                            power_d_fraction *= math::fast_power(destination_fraction, remaining_d);
                            delta += static_cast<Volume>(chg.edgeWeightBySize(d)) * (powers_of_source_community[d] + power_d_fraction - power_d_fraction_minus);
                            biggest_d_yet = d;
                            // if (delta > best_delta) {
                            //     break;
                            // }
                        }
                        ASSERT((vol_c_minus_vol_v > vol_destination_minus && exp_edge_contribution < 0.0L)
                            || (vol_c_minus_vol_v < vol_destination_minus&& exp_edge_contribution > 0.0L)
                            || (vol_c_minus_vol_v == vol_destination_minus));
                    }
                    if (delta < best_delta) {
                        best_delta = delta;
                        best_community = community;
                    }
                }
                //utils::Timer::instance().stop_timer("exp_edge_contribution");
                // #############################################################
                std::sort(ranking_after_km1.local().begin(), ranking_after_km1.local().end(), [&](const auto a, const auto b) {
                    return a.second < b.second;
                    });
                PartitionID i = 0;
                if (best_community != comm_v) {
                    while (ranking_after_km1.local()[i].first != best_community) {
                        ++i;
                    }
                    distance.push_back(i);
                    if (ranking_after_km1.local()[i].second > 0.0L) {
                        ++positive_edgeContribution_count;
                    }
                    com_neighbours.push_back(community_edge_contribution.size());
                }
                ranking_after_km1.local().clear();
                if (community_edge_contribution[best_community] > 0) {
                    ++positive_edgeContribution_count;
                }
                // #############################################################
                community_edge_contribution.clear();
                exp_edge_contribution_time += (tbb::tick_count::now() - t).seconds();
                //utils::Timer::instance().stop_timer("calculate_best_move");
                return best_community;
            }

            // ! executes the given move
            KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void makeMove(ds::CommunityHypergraph& chg, parallel::scalable_vector<HypernodeID>& communities, const HypernodeID node_to_move, const PartitionID destination_community) {
                const PartitionID source_community = communities[node_to_move];
                ASSERT(source_community != destination_community);
                _community_volumes[destination_community] += chg.nodeVolume(node_to_move);
                _community_volumes[source_community] -= chg.nodeVolume(node_to_move);
                communities[node_to_move] = destination_community;
                for (const HyperedgeID& he : chg.incidentEdges(node_to_move)) {
                    chg.addCommunityToHyperedge(he, destination_community);
                    chg.removeCommunityFromHyperedge(he, source_community);
                }
            }

            bool localMoving(ds::CommunityHypergraph& chg, parallel::scalable_vector<HypernodeID>& communities) {
                parallel::scalable_vector<HypernodeID> nodes(chg.initialNumNodes());

                tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID i) {
                    communities[i] = i;
                    nodes[i] = i;
                    _community_volumes[i].store(chg.nodeVolume(i));
                    });
                bool changed_clustering = false;
                size_t nr_nodes_moved = chg.initialNumNodes();
                for (size_t round = 0;
                    nr_nodes_moved >= _context.preprocessing.community_detection.min_vertex_move_fraction * chg.initialNumNodes()
                    && round < _context.preprocessing.community_detection.max_pass_iterations; ++round) {

                    if (!_deactivate_random) {
                        utils::Randomize::instance().parallelShuffleVector(nodes, 0UL, nodes.size());
                    }

                    tbb::enumerable_thread_specific<size_t> local_nr_nodes_moved(0);
                    tbb::parallel_for(0UL, nodes.size(), [&](const size_t i) {
                        PartitionID destination_community;
                        const HypernodeID node_to_move = nodes[i];
                        const PartitionID source_community = communities[node_to_move];
                        const size_t map_size = ratingsFitIntoSmallMap(chg, node_to_move);
                        if (!map_size) {
                            //LOG << "small";
                            destination_community = calculateBestMove(chg, communities, node_to_move, _small_edge_contribution_map.local());
                        } else {
                            //LOG << "large";
                            LargeTmpRatingMap& large_map = _large_edge_contribution_map.local();
                            large_map.setMaxSize(map_size);
                            destination_community = calculateBestMove(chg, communities, node_to_move, large_map);
                        }
                        if (destination_community != source_community) {
                            makeMove(chg, communities, node_to_move, destination_community);
                            ++local_nr_nodes_moved.local();
                        }
                        });
                    nr_nodes_moved = local_nr_nodes_moved.combine(std::plus<>());
                    changed_clustering |= nr_nodes_moved > 0;
                }
                return changed_clustering;
            }

            void initializeCommunityVolumes(const ds::CommunityHypergraph& chg, const parallel::scalable_vector<HypernodeID>& communities) {
                tbb::parallel_for(0UL, _community_volumes.size(), [&](const size_t i) {
                    _community_volumes[i].store(0);
                    });
                tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID i) {
                    _community_volumes[communities[i]] += chg.nodeVolume(i);
                    });
            }

            //! Iterator for the Community volumes
            IteratorRange<CommunityVolumeIterator> communityVolumes(const ds::CommunityHypergraph& chg) const {
                return IteratorRange<CommunityVolumeIterator>(_community_volumes.cbegin(), _community_volumes.cbegin() + chg.initialNumNodes());
            }

            parallel::AtomicWrapper<size_t> positive_edgeContribution_count;
            tbb::enumerable_thread_specific<parallel::scalable_vector<std::pair<PartitionID, Volume>>> ranking_after_km1;
            //tbb::enumerable_thread_specific<parallel::scalable_vector<Volume>> ranking_end;
            // parallel::AtomicWrapper<size_t> exp_makes_it_bad;
            tbb::concurrent_vector<int> distance;
            tbb::concurrent_vector<int> com_neighbours;
            //parallel::AtomicWrapper<size_t> overall_checks;
            //parallel::AtomicWrapper<size_t> pruned_by_old;
            parallel::AtomicWrapper<double> edge_contribution_time;
            parallel::AtomicWrapper<double> exp_edge_contribution_time;

private:

    KAHYPAR_ATTRIBUTE_ALWAYS_INLINE size_t ratingsFitIntoSmallMap(const ds::CommunityHypergraph& chg, const HypernodeID v) const {
        const bool use_vertex_degree_sampling =
            _context.coarsening.vertex_degree_sampling_threshold != std::numeric_limits<size_t>::max();
        const size_t vertex_degree_bounded_rating_map_size = use_vertex_degree_sampling ?
            3UL * _context.coarsening.vertex_degree_sampling_threshold : std::numeric_limits<size_t>::max();
        const size_t cache_efficient_rating_map_size = CacheEfficientRatingMap::MAP_SIZE;
        const size_t size_of_smaller_rating_map = std::min(
            vertex_degree_bounded_rating_map_size, cache_efficient_rating_map_size);

        // In case the current number of nodes is smaller than size
        // of the cache-efficient sparse map, the large tmp rating map
        // consumes less memory
        if (chg.initialNumNodes() < size_of_smaller_rating_map) {
            return chg.initialNumNodes();
        }

        // Compute estimation for the upper bound of neighbors of u
        HypernodeID ub_neighbors_v = 0;
        for (const HyperedgeID& he : chg.incidentEdges(v)) {
            const HypernodeID edge_size = chg.edgeSize(he);
            // Ignore large hyperedges
            ub_neighbors_v += edge_size;
            // If the number of estimated neighbors is greater than the size of the cache efficient rating map / 3, we
            // use the large sparse map. The division by 3 also ensures that the fill grade
            // of the cache efficient sparse map would be small enough such that linear probing
            // is fast.
            if (ub_neighbors_v > cache_efficient_rating_map_size / 3UL) {
                return std::min(vertex_degree_bounded_rating_map_size, static_cast<size_t>(chg.initialNumNodes()));
            }
        }
        return 0;
    }

    LargeTmpRatingMap construct_large_edge_contribution_map(const size_t num_nodes) {
        return LargeTmpRatingMap(3UL * std::min(num_nodes, _vertex_degree_sampling_threshold));
    }

    const size_t _community_neighbour_sampling_threshold;
    const size_t _hyperedge_size_caching_threshold;
    const size_t _vertex_degree_sampling_threshold;
    // ! contains Hyperparameters for the algorithms
    const Context& _context;

    // ! reciprocal total volume
    Volume _reciprocal_vol_total = 0.0L;

    // contains the edge contribution for each neighbouring community
    ThreadLocalCacheEfficientRatingMap _small_edge_contribution_map;
    ThreadLocalLargeTmpRatingMap _large_edge_contribution_map;

    // map for finding the neighbours of an edge (without counting them twice)
    ThreadLocalCacheEfficientRatingMap _community_neighbours_of_edge;

    // ! contains (vol_V - vol(C)+vol(v))^d - (vol(V)-vol(C))^d for all valid edgesizes d
    // ! Note the values in here are not cleared after each call to calculateBestMove
    tbb::enumerable_thread_specific<parallel::scalable_vector<Volume>> _powers_of_source_community;

    // ! volumes of each community
    parallel::scalable_vector<AtomicHyperedgeWeight> _community_volumes;

    // ! deactivates random node order in local moving
    const bool _deactivate_random;


    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunity0);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunity1);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunity2);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunity3);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunity4);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunity5);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunity6);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges0);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges1);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges2);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges3);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges4);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges5);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges6);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges0);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges1);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges2);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges3);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges4);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges5);
    FRIEND_TEST(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges6);
};
}

namespace mt_kahypar::metrics {
Volume hyp_modularity(const ds::CommunityHypergraph& hypergraph, const parallel::scalable_vector<HypernodeID>& communities, const community_detection::HypergraphLocalMovingModularity& hlmm);
}