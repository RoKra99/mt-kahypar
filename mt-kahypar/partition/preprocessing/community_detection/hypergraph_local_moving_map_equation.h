#pragma once

#include "mt-kahypar/datastructures/community_hypergraph.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/datastructures/sparse_map.h"

#include "gtest/gtest_prod.h"
#include "mt-kahypar/utils/floating_point_comparisons.h"

namespace mt_kahypar::metrics {
Probability hyp_map_equation(const ds::CommunityHypergraph& chg, const parallel::scalable_vector<HypernodeID>& communities);
}

namespace mt_kahypar::community_detection {

class HypergraphLocalMovingMapEquation {

    // ! represents a single map equation improving move and the corresponding changes in probabilities
    struct Move {
        HypernodeID destination_community;
        Probability delta_source;
        Probability delta_destination;
        Probability delta;
    };

private:

    using AtomicHyperedgeWeight = parallel::AtomicWrapper<HyperedgeWeight>;
    using AtomicProbability = parallel::AtomicWrapper<Probability>;
    using LargeTmpRatingMap = ds::SparseMap<HypernodeID, HyperedgeWeight>;
    using CacheEfficientRatingMap = ds::FixedSizeSparseMap<HypernodeID, HyperedgeWeight>;
    using ThreadLocalCacheEfficientRatingMap = tbb::enumerable_thread_specific<CacheEfficientRatingMap>;
    using ThreadLocalLargeTmpRatingMap = tbb::enumerable_thread_specific<LargeTmpRatingMap>;
    //TODO: Temporary
    using LargeTmpDoubleMap = ds::SparseMap<HypernodeID, Probability>;
    using CacheEfficentDoubleMap = ds::FixedSizeSparseMap<HypernodeID, Probability>;
    using ThreadLocalCacheEfficientDoubleMap = tbb::enumerable_thread_specific<CacheEfficentDoubleMap>;
    using ThreadLocalLargeTmpDoubleMap = tbb::enumerable_thread_specific<LargeTmpDoubleMap>;

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
        _small_deltas_mul_vol_total(0),
        _large_deltas_mul_vol_total([&] {
        return construct_large_delta_map(chg.initialNumNodes());
    }),
        _deactivate_random(deactivate_random) {}


    // ! executes the local moving for map equation
    bool localMoving(ds::CommunityHypergraph& chg, parallel::scalable_vector<HypernodeID>& communities) {
        parallel::scalable_vector<HypernodeID> nodes(chg.initialNumNodes());

        tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID i) {
            communities[i] = i;
            nodes[i] = i;
            _community_volumes[i].store(chg.nodeVolume(i));
            _community_exit_probability_mul_vol_total[i].store(0.0L);
        });

        tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeWeight>> overlap_local(chg.initialNumNodes(), 0);
        tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeID>> neighbouring_communities;

        tbb::parallel_for(0U, chg.initialNumEdges(), [&](const HyperedgeID he) {
            //for (HyperedgeID he = 0; he < chg.initialNumEdges(); ++he) {
            for (const auto& mp : chg.multipins(he)) {
                const size_t community_mp = communities[mp.id];
                if (overlap_local.local()[community_mp] == 0) {
                    neighbouring_communities.local().push_back(community_mp);
                }
                overlap_local.local()[community_mp] += mp.multiplicity;
            }

            for (const HypernodeID comm : neighbouring_communities.local()) {
                const HypernodeWeight pincount_in_edge = overlap_local.local()[comm];
                const HypernodeWeight edge_size = static_cast<HypernodeWeight>(chg.edgeSize(he));
                //TODO: Not sure if / (edge_size - 1) is better (since that results in an equal model to the original map equation)
                _community_exit_probability_mul_vol_total[comm] += static_cast<Probability>(pincount_in_edge * chg.edgeWeight(he) * (edge_size - pincount_in_edge)) / edge_size;//(edge_size - 1);
                overlap_local.local()[comm] = 0;
            }
            neighbouring_communities.local().clear();
            //}
            // _sum_exit_probability_mul_vol_total = 0.0L;
            // for (const auto& e : _community_exit_probability_mul_vol_total) {
            //     _sum_exit_probability_mul_vol_total += e;
            // }
        });
        _sum_exit_probability_mul_vol_total.store(0.0L);
        tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID hn) {
            _sum_exit_probability_mul_vol_total += _community_exit_probability_mul_vol_total[hn];
        });

        bool changed_clustering = false;
        size_t nr_nodes_moved = chg.initialNumNodes();
        for (size_t round = 0;
            nr_nodes_moved >= _context.preprocessing.community_detection.min_vertex_move_fraction * chg.initialNumNodes()
            && round < _context.preprocessing.community_detection.max_pass_iterations; ++round) {

            nr_nodes_moved = 0;
            if (!_deactivate_random) {
                utils::Randomize::instance().parallelShuffleVector(nodes, 0UL, nodes.size());
            }

            //for (size_t i = 0; i < nodes.size(); ++i) {
            tbb::parallel_for(0UL, nodes.size(), [&](const size_t i) {
                Move move;
                const HypernodeID node_to_move = nodes[i];
                const HypernodeID source_community = communities[node_to_move];
                const size_t map_size = ratingsFitIntoSmallMap(chg, node_to_move);
                if (!map_size) {
                    move = calculateBestMove(chg, communities, node_to_move, _small_overlap_map.local(), _small_deltas_mul_vol_total.local());
                } else {
                    LargeTmpRatingMap& large_map = _large_overlap_map.local();
                    large_map.setMaxSize(map_size);
                    LargeTmpDoubleMap& large_deltas = _large_deltas_mul_vol_total.local();
                    large_deltas.setMaxSize(map_size);
                    move = calculateBestMove(chg, communities, node_to_move, large_map, large_deltas);
                }
                if (move.destination_community != source_community) {
                    makeMove(chg, communities, node_to_move, move);
                    ++nr_nodes_moved;
                }
            });
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

    // ! recomputes all exit probabilities
    void recomputeExitProbability(const ds::CommunityHypergraph& chg, const parallel::scalable_vector<HypernodeID>& communities) {
        tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID hn) {
            _community_exit_probability_mul_vol_total[hn].store(0.0L);
        });
        tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeWeight>> overlap_local(chg.initialNumNodes(), 0);
        tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeID>> neighbouring_communities;
        for (HyperedgeID he = 0; he < chg.initialNumEdges(); ++he) {
            for (const auto& mp : chg.multipins(he)) {
                const size_t community_mp = communities[mp.id];
                if (overlap_local.local()[community_mp] == 0) {
                    neighbouring_communities.local().push_back(community_mp);
                }
                overlap_local.local()[community_mp] += mp.multiplicity;
            }

            for (const HypernodeID comm : neighbouring_communities.local()) {
                const HypernodeWeight pincount_in_edge = overlap_local.local()[comm];
                const HypernodeWeight edge_size = static_cast<HypernodeWeight>(chg.edgeSize(he));
                //TODO: Not sure if / (edge_size - 1) is better (since that results in an equal model to the original map equation)
                _community_exit_probability_mul_vol_total[comm] += static_cast<Probability>(pincount_in_edge * chg.edgeWeight(he) * (edge_size - pincount_in_edge)) / edge_size;//(edge_size - 1);
                overlap_local.local()[comm] = 0;
            }
            neighbouring_communities.local().clear();
        }
        _sum_exit_probability_mul_vol_total.store(0.0L);
        for (const auto& e : _community_exit_probability_mul_vol_total) {
            _sum_exit_probability_mul_vol_total += e;
        }
    }

    // ! recalculates and sets the exit probabilities of the source and destination community and then recomputes the sum
    void recomputeExitProbability(const ds::CommunityHypergraph& chg, const parallel::scalable_vector<HypernodeID>& communities, const HypernodeID source, const HypernodeID destination) {
        parallel::scalable_vector<Probability> exit_probs(chg.initialNumNodes(), 0.0);
        // tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID hn) {
        //     exit_probs[hn].store(0.0L);
        // });
        tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeWeight>> overlap_local(chg.initialNumNodes(), 0);
        tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeID>> neighbouring_communities;
        for (HyperedgeID he = 0; he < chg.initialNumEdges(); ++he) {
            for (const auto& mp : chg.multipins(he)) {
                const size_t community_mp = communities[mp.id];
                if (overlap_local.local()[community_mp] == 0) {
                    neighbouring_communities.local().push_back(community_mp);
                }
                overlap_local.local()[community_mp] += mp.multiplicity;
            }

            for (const HypernodeID comm : neighbouring_communities.local()) {
                if (comm == source || comm == destination) {
                    const HypernodeWeight pincount_in_edge = overlap_local.local()[comm];
                    const HypernodeWeight edge_size = static_cast<HypernodeWeight>(chg.edgeSize(he));
                    //TODO: Not sure if / (edge_size - 1) is better (since that results in an equal model to the original map equation)
                    exit_probs[comm] += static_cast<Probability>(pincount_in_edge * chg.edgeWeight(he) * (edge_size - pincount_in_edge)) / edge_size; //(edge_size - 1);
                }
                overlap_local.local()[comm] = 0;
            }
            neighbouring_communities.local().clear();
        }

        _community_exit_probability_mul_vol_total[source].store(exit_probs[source]);
        _community_exit_probability_mul_vol_total[destination].store(exit_probs[destination]);

        tbb::enumerable_thread_specific<Probability> sum_exit_prob_local(0.0);
        tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID com) {
            sum_exit_prob_local.local() += _community_exit_probability_mul_vol_total[com];
        });
        _sum_exit_probability_mul_vol_total.store(sum_exit_prob_local.combine(std::plus<>()));
    }

    // ! calculate the move wich improves the map equation the most and resturn this move
    template<typename Map, typename Map2>
    Move calculateBestMove(ds::CommunityHypergraph& chg, parallel::scalable_vector<HypernodeID>& communities, const HypernodeID v, Map& overlap, Map2& deltas) {
        ASSERT(overlap.size() == 0);
        ASSERT(deltas.size() == 0);
        Probability community_independent_part = 0.0L;
        const HypernodeID comm_v = communities[v];
        const HyperedgeWeight vol_v = chg.nodeVolume(v);
        HypernodeWeight pin_count_v = 0;
        //Probability sum_of_reciprocals = 0.0;
        for (const HyperedgeID& he : chg.incidentEdges(v)) {
            ASSERT(overlap.size() == 0);
            const Probability reciprocal_edge_size = 1.0L / chg.edgeSize(he); // (chg.edgeSize(he) - 1);
            const HypernodeWeight edge_weight = chg.edgeWeight(he);
            for (const auto& mp : chg.multipins(he)) {
                const HypernodeID community = communities[mp.id];
                if (mp.id == v) {
                    pin_count_v = mp.multiplicity;
                }
                overlap[community] += mp.multiplicity;
            }

            for (const auto& e : overlap) {
                if (e.key != comm_v) {
                    deltas[e.key] += static_cast<Probability>(edge_weight * pin_count_v * -2 * e.value) * reciprocal_edge_size;
                } else {
                    deltas[e.key] += static_cast<Probability>(edge_weight * pin_count_v * (2 * e.value - pin_count_v)) * reciprocal_edge_size;
                }
            }
            overlap.clear();
            community_independent_part -= static_cast<Probability>(edge_weight * pin_count_v * pin_count_v) * reciprocal_edge_size;
            //sum_of_reciprocals += edge_weight * pin_count_v * reciprocal_edge_size;
        }

        // the formula only works for communities with volume != 0 therefore we have to handle this edge case explicitly
        Probability delta_source = deltas[comm_v] - vol_v; //- sum_of_reciprocals;
        //if (vol_v == _community_volumes[comm_v]) {
            //Probability should = -_community_exit_probability_mul_vol_total[comm_v];
            //ASSERT(mt_kahypar::math::are_almost_equal_ld(should, delta_source, 1e-10L));
            //const int i = 1 + 1;
        //    delta_source = -_community_exit_probability_mul_vol_total[comm_v];
        //}

        const auto plogp_rel = [this](Probability p) -> Probability {
            if (p > 0.0L) {
                p *= _reciprocal_vol_total;
                return p * log2(p);
            }
            return 0.0L;
        };

        const Probability remain_plogp_sum_exit_prob = plogp_rel(_sum_exit_probability_mul_vol_total);
        // summands of map equation if we don't move the node at all and only depend on the source node/community
        // TODO: store this?
        const Probability remain_plogp_exit_prob_source = plogp_rel(_community_exit_probability_mul_vol_total[comm_v]);
        const Probability remain_plogp_exit_prob_plus_vol_source = plogp_rel(_community_exit_probability_mul_vol_total[comm_v] + _community_volumes[comm_v]);

        // summands if we move the node and only depend on the source node/community
        const Probability move_plogp_exit_prob_source = plogp_rel(_community_exit_probability_mul_vol_total[comm_v] + delta_source);
        const Probability move_plogp_exit_prob_plus_vol_source = plogp_rel(_community_exit_probability_mul_vol_total[comm_v] + delta_source + _community_volumes[comm_v] - vol_v);

        Move move;
        move.destination_community = comm_v;
        move.delta = 0.0L;
        move.delta_source = delta_source;

        for (const auto& e : deltas) {
            if (e.key == comm_v) { continue; }

            const HypernodeID destination_community = e.key;
            const Probability delta_destination = community_independent_part + e.value + vol_v; //+ sum_of_reciprocals;
            // summands if we don't move the node at all and depend on destination
            const Probability remain_plogp_exit_prob_destination = plogp_rel(_community_exit_probability_mul_vol_total[destination_community]);
            const Probability remain_plogp_exit_plus_vol_destination = plogp_rel(_community_exit_probability_mul_vol_total[destination_community] + _community_volumes[destination_community]);

            //summands if we move the node and only depend on the destination node/community
            const Probability move_plogp_exit_prob_destination = plogp_rel(_community_exit_probability_mul_vol_total[destination_community] + delta_destination);
            const Probability move_plogp_exit_prob_plus_vol_destination = plogp_rel(_community_exit_probability_mul_vol_total[destination_community] + delta_destination + _community_volumes[destination_community] + vol_v);
            // summands if we move and depend on source and destination
            const Probability move_plogp_sum_exit_prob = plogp_rel(_sum_exit_probability_mul_vol_total + delta_source + delta_destination);

            const Probability delta = remain_plogp_sum_exit_prob - move_plogp_sum_exit_prob
                - 2 * (remain_plogp_exit_prob_source + remain_plogp_exit_prob_destination - move_plogp_exit_prob_source - move_plogp_exit_prob_destination)
                + (remain_plogp_exit_prob_plus_vol_source + remain_plogp_exit_plus_vol_destination - move_plogp_exit_prob_plus_vol_source - move_plogp_exit_prob_plus_vol_destination);

            // remain - delta, to minimize map equation choose large deltas > 0
            if (delta > move.delta) {
                move.delta = delta;
                move.destination_community = destination_community;
                move.delta_destination = delta_destination;
            }
        }
        deltas.clear();
        return move;
    }

    // ! executes a the given move and also updates the data structures accordingly
    KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void makeMove(ds::CommunityHypergraph& chg, parallel::scalable_vector<HypernodeID>& communities, const HypernodeID node_to_move, const Move& move) {
        ASSERT(communities[node_to_move] != move.destination_community);
        ASSERT(move.delta > 0);
#ifdef KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS
        const Probability before = metrics::hyp_map_equation(chg, communities);
        LOG << "before" << before;

#endif

        const HypernodeID source_community = communities[node_to_move];
        const HyperedgeWeight vol_v = chg.nodeVolume(node_to_move);
        communities[node_to_move] = move.destination_community;
        _community_volumes[source_community] -= vol_v;
        _community_volumes[move.destination_community] += vol_v;
        // for (const HyperedgeID& he : chg.incidentEdges(node_to_move)) {
        //     // remove has to be before add to ensure the amount of distinct communities per edge < edgeSize
        //     chg.removeCommunityFromHyperedge(he, source_community);
        //     chg.addCommunityToHyperedge(he, move.destination_community);
        // }

        //recomputeExitProbability(chg, communities, source_community, move.destination_community);
        _community_exit_probability_mul_vol_total[source_community] += move.delta_source;
        _community_exit_probability_mul_vol_total[move.destination_community] += move.delta_destination;
        //ASSERT(_community_exit_probability_mul_vol_total[source_community] >= 0.0L, V(_community_exit_probability_mul_vol_total[source_community]));
        //ASSERT(_community_exit_probability_mul_vol_total[move.destination_community] >= 0.0L, V(_community_exit_probability_mul_vol_total[move.destination_community]));
        _sum_exit_probability_mul_vol_total += move.delta_destination + move.delta_source;
        ASSERT(_sum_exit_probability_mul_vol_total <= chg.totalVolume());
        ASSERT(_sum_exit_probability_mul_vol_total >= 0.0L);
#ifdef KAHYPAR_ENABLE_HEAVY_PREPROCESSING_ASSERTIONS
        LOG << std::fixed << "delta" << std::setprecision(15) << move.delta;
        const Probability after = metrics::hyp_map_equation(chg, communities);

        LOG << std::fixed << std::setprecision(15) << V(before - after);
        ASSERT(mt_kahypar::math::are_almost_equal_ld(move.delta, before - after, 1e-10L));
        //ASSERT(before - after > 0.0L);
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

    LargeTmpRatingMap construct_large_overlap_map(const size_t num_nodes) {
        return LargeTmpRatingMap(3UL * std::min(num_nodes, _vertex_degree_sampling_threshold));
    }

    LargeTmpDoubleMap construct_large_delta_map(const size_t num_nodes) {
        return LargeTmpDoubleMap(3UL * std::min(num_nodes, _vertex_degree_sampling_threshold));
    }

    const size_t _vertex_degree_sampling_threshold;

    const size_t _hyperedge_size_caching_threshold;

    // ! contains parameters for the algorithms
    const Context _context;

    // ! reciprocal total volume
    const Probability _reciprocal_vol_total;

    // ! volumes of each community
    parallel::scalable_vector<AtomicHyperedgeWeight> _community_volumes;

    // ! contains the exit probability of each community multiplied by vol(V) (qi_exit * vol(V))
    parallel::scalable_vector<AtomicProbability> _community_exit_probability_mul_vol_total;

    // ! the sum of all exit probabilities multiplied by vol(V)
    AtomicProbability _sum_exit_probability_mul_vol_total;

    // contains the overlap for each neighbouring community
    ThreadLocalCacheEfficientRatingMap _small_overlap_map;
    ThreadLocalLargeTmpRatingMap _large_overlap_map;

    // map for finding the neighbours of an edge (without counting them twice)
    ThreadLocalCacheEfficientRatingMap _community_neighbours_of_edge;

    //TODO: This is a temporary solution while there is no datastructure for the pincounts of a community in an edge
    ThreadLocalCacheEfficientDoubleMap _small_deltas_mul_vol_total;
    ThreadLocalLargeTmpDoubleMap _large_deltas_mul_vol_total;

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