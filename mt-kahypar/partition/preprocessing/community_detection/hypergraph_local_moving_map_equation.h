#pragma once

#include "mt-kahypar/datastructures/community_hypergraph.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/datastructures/sparse_map.h"

#include "gtest/gtest_prod.h"
#include "mt-kahypar/utils/floating_point_comparisons.h"

namespace mt_kahypar::metrics {
double hyp_map_equation(const ds::CommunityHypergraph& chg, const parallel::scalable_vector<HypernodeID>& communities);
}

namespace mt_kahypar::community_detection {

class HypergraphLocalMovingMapEquation {

    struct Move {
        HypernodeID destination_community;
        size_t delta_source;
        size_t delta_destination;
        double delta;
    };

private:

    using AtomicHyperedgeWeight = parallel::AtomicWrapper<HyperedgeWeight>;
    using AtomicDouble = parallel::AtomicWrapper<double>;
    using LargeTmpRatingMap = ds::SparseMap<PartitionID, HyperedgeWeight>;
    using CacheEfficientRatingMap = ds::FixedSizeSparseMap<PartitionID, HyperedgeWeight>;
    using ThreadLocalCacheEfficientRatingMap = tbb::enumerable_thread_specific<CacheEfficientRatingMap>;
    using ThreadLocalLargeTmpRatingMap = tbb::enumerable_thread_specific<LargeTmpRatingMap>;

public:

    HypergraphLocalMovingMapEquation(const ds::CommunityHypergraph& chg, const Context& context, const bool deactivate_random = false) :
        _vertex_degree_sampling_threshold(context.coarsening.vertex_degree_sampling_threshold),
        _hyperedge_size_caching_threshold(context.preprocessing.community_detection.hyperedge_size_caching_threshold),
        _context(context),
        _reciprocal_vol_total(1.0 / chg.totalVolume()),
        _community_volumes(chg.initialNumNodes()),
        _community_exit_probability_mul_vol_total(chg.initialNumNodes()),
        _sum_exit_probability_mul_vol_total(0.0),
        _small_overlap_map(0),
        _large_overlap_map([&] {
        return construct_large_overlap_map(chg.initialNumNodes());
            }),
        _community_neighbours_of_edge(0),
                _deactivate_random(deactivate_random) {
        initializeExitProbabilities(chg);
    }


            bool localMoving(ds::CommunityHypergraph& chg, parallel::scalable_vector<HypernodeID>& communities) {
                parallel::scalable_vector<HypernodeID> nodes(chg.initialNumNodes());

                for (size_t i = 0; i < chg.initialNumNodes(); ++i) {
                    communities[i] = i;
                    nodes[i] = i;
                    _community_volumes[i].store(chg.nodeVolume(i));
                    _community_exit_probability_mul_vol_total[i].store(0);
                }

                tbb::enumerable_thread_specific<parallel::scalable_vector<size_t>> overlap_local(chg.initialNumNodes(), 0);
                tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeID>> neighbouring_communities;

                // TODO: recalculate exit probabilities (later do contractions)
                // tbb::parallel_for(0U, chg.initialNumEdges(), [&](const HyperedgeID he) {
                //     const size_t edge_size = chg.edgeSize(he);
                //     const size_t edge_weight = chg.edgeWeight(he);

                //     for (const auto& mp : chg.multipins(he)) {
                //         const HypernodeID community = communities[mp.id];
                //         const size_t pincount = mp.multiplicity;
                //         if (overlap_local.local()[community] == 0) {
                //             neighbouring_communities.local().push_back(community);
                //         }
                //         overlap_local.local()[community] += pincount;
                //     }
                //     for (const HypernodeID community : neighbouring_communities.local()) {
                //         _community_exit_probability_mul_vol_total[community] += (edge_size - overlap_local.local()[community]) * edge_weight;
                //         overlap_local.local()[community] = 0;
                //     }
                //     neighbouring_communities.local().clear();
                //     });

                bool changed_clustering = false;
                size_t nr_nodes_moved = chg.initialNumNodes();
                for (size_t round = 0;
                    nr_nodes_moved >= _context.preprocessing.community_detection.min_vertex_move_fraction * chg.initialNumNodes()
                    && round < _context.preprocessing.community_detection.max_pass_iterations; ++round) {

                    nr_nodes_moved = 0;
                    if (!_deactivate_random) {
                        utils::Randomize::instance().parallelShuffleVector(nodes, 0UL, nodes.size());
                    }

                    for (size_t i = 0; i < nodes.size(); ++i) {
                        Move move;
                        const HypernodeID node_to_move = nodes[i];
                        const HypernodeID source_community = communities[node_to_move];
                        const size_t map_size = ratingsFitIntoSmallMap(chg, node_to_move);
                        if (!map_size) {
                            move = calculateBestMove(chg, communities, node_to_move, _small_overlap_map.local());
                        } else {
                            LargeTmpRatingMap& large_map = _large_overlap_map.local();
                            large_map.setMaxSize(map_size);
                            move = calculateBestMove(chg, communities, node_to_move, large_map);
                        }
                        if (move.destination_community != source_community) {
                            makeMove(chg, communities, node_to_move, move);
                            ++nr_nodes_moved;
                        }
                    }
                    changed_clustering |= nr_nodes_moved > 0;
                }
                return changed_clustering;
            }

            // ! initializes communityVolumes just for testing purposes
            void initializeCommunityVolumes(const ds::CommunityHypergraph& chg, const parallel::scalable_vector<HypernodeID>& communities) {
                tbb::parallel_for(0UL, _community_volumes.size(), [&](const size_t i) {
                    _community_volumes[i].store(0);
                    });
                tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID i) {
                    _community_volumes[communities[i]] += chg.nodeVolume(i);
                    });
            }

private:

    template<typename Map>
    Move calculateBestMove(ds::CommunityHypergraph& chg, parallel::scalable_vector<HypernodeID>& communities, const HypernodeID v, Map& overlap) {
        ASSERT(overlap.size() == 0);
        const PartitionID comm_v = communities[v];
        // sum of all edgeweights incident to v
        HyperedgeWeight sum_of_edgeweights = 0;


        const HyperedgeWeight vol_v = chg.nodeVolume(v);
        const HyperedgeWeight delta_source = -overlap[comm_v] + vol_v;

        const auto plogp_rel = [this](double p) -> double {
            if (p > 0.0) {
                p *= _reciprocal_vol_total;
                return p * log2(p);
            }
            return 0.0;
        };

        const double remain_plogp_sum_exit_prob = plogp_rel(_sum_exit_probability_mul_vol_total);
        // summands of map equation if we don't move the node at all and only depend on the source node/community
        // TODO: store this?
        const double remain_plogp_exit_prob_source = plogp_rel(_community_exit_probability_mul_vol_total[comm_v]);
        const double remain_plogp_exit_prob_plus_vol_source = plogp_rel(_community_exit_probability_mul_vol_total[comm_v] + _community_volumes[comm_v]);

        // summands if we move the node and only depend on the source node/community
        const double move_plogp_exit_prob_source = plogp_rel(_community_exit_probability_mul_vol_total[comm_v] + delta_source);
        const double move_plogp_exit_prob_plus_vol_source = plogp_rel(_community_exit_probability_mul_vol_total[comm_v] + delta_source + _community_volumes[comm_v] - vol_v);

        Move move;
        move.destination_community = comm_v;
        move.delta = 0.0;
        move.delta_source = delta_source;

        for (const auto& e : overlap) {
            if (e.key == comm_v) continue;

            const HypernodeID destination_community = e.key;
            const size_t delta_destination = sum_of_edgeweights - e.value - vol_v;
            // summands if we don't move the node at all and depend on destination
            //LOG << "testing" << _community_exit_probability_mul_vol_total[destination_community];
            const double remain_plogp_exit_prob_destination = plogp_rel(_community_exit_probability_mul_vol_total[destination_community]);
            const double remain_plogp_exit_plus_vol_destination = plogp_rel(_community_exit_probability_mul_vol_total[destination_community] + _community_volumes[destination_community]);

            //summands if we move the node and only depend on the destination node/community
            const double move_plogp_exit_prob_destination = plogp_rel(_community_exit_probability_mul_vol_total[destination_community] + delta_destination);
            const double move_plogp_exit_prob_plus_vol_destination = plogp_rel(_community_exit_probability_mul_vol_total[destination_community] + delta_destination + _community_volumes[destination_community] + vol_v);
            // summands if we move and depend on source and destination
            const double move_plogp_sum_exit_prob = plogp_rel(_sum_exit_probability_mul_vol_total + delta_source + delta_destination);

            const double delta = remain_plogp_sum_exit_prob - move_plogp_sum_exit_prob
                - 2 * (remain_plogp_exit_prob_source + remain_plogp_exit_prob_destination - move_plogp_exit_prob_source - move_plogp_exit_prob_destination)
                + (remain_plogp_exit_prob_plus_vol_source + remain_plogp_exit_plus_vol_destination - move_plogp_exit_prob_plus_vol_source - move_plogp_exit_prob_plus_vol_destination);

            // remain - delta, to minimize map equation choose large deltas > 0
            if (delta > move.delta) {
                //LOG << "improving move";
                //LOG << remain_plogp_sum_exit_prob - move_plogp_sum_exit_prob;
                //LOG << remain_plogp_exit_prob_source + remain_plogp_exit_prob_destination - move_plogp_exit_prob_source - move_plogp_exit_prob_destination;
                //LOG << remain_plogp_exit_prob_plus_vol_source + remain_plogp_exit_plus_vol_destination - move_plogp_exit_prob_plus_vol_source - move_plogp_exit_prob_plus_vol_destination;
                //LOG << delta;
                move.delta = delta;
                move.destination_community = destination_community;
                move.delta_destination = delta_destination;
            }
        }
        overlap.clear();
        return move;
    }

    KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void makeMove(ds::CommunityHypergraph& chg, parallel::scalable_vector<HypernodeID>& communitties, const HypernodeID node_to_move, const Move& move) {
        ASSERT(communitties[node_to_move] != move.destination_community);
#ifdef KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS
        const double before = metrics::hyp_map_equation(chg, communitties);
        LOG << "before" << before;

#endif

        const HypernodeID source_community = communitties[node_to_move];
        const HyperedgeWeight vol_v = chg.nodeVolume(node_to_move);
        communitties[node_to_move] = move.destination_community;
        _community_volumes[source_community] -= vol_v;
        _community_volumes[move.destination_community] += vol_v;
        for (const HyperedgeID& he : chg.incidentEdges(node_to_move)) {
            // remove has to be before add to ensure the amount of distinct communities per edge < edgeSize
            chg.removeCommunityFromHyperedge(he, source_community);
            chg.addCommunityToHyperedge(he, move.destination_community);
        }

        _community_exit_probability_mul_vol_total[source_community] += move.delta_source;
        _community_exit_probability_mul_vol_total[move.destination_community] += move.delta_destination;
        _sum_exit_probability_mul_vol_total += move.delta_destination + move.delta_source;
#ifdef KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS
        LOG << "delta" << move.delta;
        const double after = metrics::hyp_map_equation(chg, communitties);
        LOG << "after" << after;
        ASSERT(mt_kahypar::math::are_almost_equal_ld(move.delta, before - after, 1e-8L));
#endif
    }

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

    // ! initializes the exit probabilities for each community (expects each edge to only contain one pin of each community)
    void initializeExitProbabilities(const ds::CommunityHypergraph& chg) {
        tbb::parallel_for(0U, chg.initialNumEdges(), [&](const HyperedgeID he) {
            for (const auto& mp : chg.multipins(he)) {
                _community_neighbours_of_edge.local()[mp.id] += 1;
            }

            for (const auto& e : _community_neighbours_of_edge.local()) {
                const HypernodeID comm = e.key;
                const HypernodeWeight edge_size = static_cast<HypernodeWeight>(chg.edgeSize(he));
                _community_exit_probability_mul_vol_total[comm] += static_cast<double>(chg.edgeWeight(he) * (edge_size - 1)) / edge_size;
            }
            _community_neighbours_of_edge.local().clear();
            });
        tbb::enumerable_thread_specific<double> sum_exit_prob_local(0.0);
        tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID hn) {
            sum_exit_prob_local.local() += _community_exit_probability_mul_vol_total[hn];
            });
        _sum_exit_probability_mul_vol_total.store(sum_exit_prob_local.combine(std::plus<>()));
    }

    LargeTmpRatingMap construct_large_overlap_map(const size_t num_nodes) {
        return LargeTmpRatingMap(3UL * std::min(num_nodes, _vertex_degree_sampling_threshold));
    }

    const size_t _vertex_degree_sampling_threshold;

    const size_t _hyperedge_size_caching_threshold;

    // ! contains parameters for the algorithms
    const Context _context;

    // ! reciprocal total volume
    const double _reciprocal_vol_total;

    // ! volumes of each community
    parallel::scalable_vector<AtomicHyperedgeWeight> _community_volumes;

    // ! contains the exit probability of each community multiplied by vol(V) (qi_exit * vol(V))
    parallel::scalable_vector<AtomicDouble> _community_exit_probability_mul_vol_total;

    // ! the sum of all exit probabilities multiplied by vol(V)
    AtomicDouble _sum_exit_probability_mul_vol_total;

    // contains the overlap for each neighbouring community
    ThreadLocalCacheEfficientRatingMap _small_overlap_map;
    ThreadLocalLargeTmpRatingMap _large_overlap_map;

    // map for finding the neighbours of an edge (without counting them twice)
    ThreadLocalCacheEfficientRatingMap _community_neighbours_of_edge;

    // ! deactivates random node order in local moving
    const bool _deactivate_random;

    FRIEND_TEST(AHypergraphLocalMovingMapEquation, InitializesTheExitProbabilities);
    FRIEND_TEST(AHypergraphLocalMovingMapEquation, InitializesTheExitProbabilities0);
    FRIEND_TEST(AHypergraphLocalMovingMapEquation, InitializesTheExitProbabilities1);
    FRIEND_TEST(AHypergraphLocalMovingMapEquation, InitializesTheExitProbabilities2);
    FRIEND_TEST(AHypergraphLocalMovingMapEquation, InitializesTheExitProbabilities3);
    FRIEND_TEST(AHypergraphLocalMovingMapEquation, InitializesTheExitProbabilities4);
    FRIEND_TEST(AHypergraphLocalMovingMapEquation, test2);
};
}