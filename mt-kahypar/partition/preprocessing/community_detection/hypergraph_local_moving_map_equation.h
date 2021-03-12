#include "mt-kahypar/datastructures/community_hypergraph.h"
#include "mt-kahypar/partition/context.h"

#include "gtest/gtest_prod.h"

namespace mt_kahypar::community_detection {

class HypergraphLocalMovingMapEquation {
private:

    using AtomicHyperedgeWeight = parallel::AtomicWrapper<HyperedgeWeight>;
    using AtomicPinCOunt = parallel::AtomicWrapper<uint32_t>;

public:

    HypergraphLocalMovingMapEquation(const ds::CommunityHypergraph& chg, const Context& context, const bool deactivate_random = false) :
        _vertex_degree_sampling_threshold(context.coarsening.vertex_degree_sampling_threshold),
        _context(context),
        _reciprocal_vol_total(1.0 / chg.totalVolume()),
        _community_volumes(chg.initialNumNodes()),
        _community_exit_probability_mul_vol_total(chg.initialNumNodes()),
        _sum_exit_probability_mul_vol_total(0.0),
        _deactivate_random(deactivate_random) {
        initializeExitProbabilities(chg);
    }


    bool localMoving(ds::CommunityHypergraph& chg, parallel::scalable_vector<HypernodeID>& communities) {
        parallel::scalable_vector<HypernodeID> nodes(chg.initialNumNodes());

        for (size_t i = 0; i < chg.initialNumNodes(); ++i) {
            communities[i] = i;
            nodes[i] = i;
            _community_volumes[i].store(chg.nodeVolume(i));
        }

        bool changed_clustering = false;
        size_t nr_nodes_moved = chg.initialNumNodes();
        for (size_t round = 0;
            nr_nodes_moved >= _context.preprocessing.community_detection.min_vertex_move_fraction * chg.initialNumNodes()
            && round < _context.preprocessing.community_detection.max_pass_iterations; ++round) {

            if (!_deactivate_random) {
                utils::Randomize::instance().parallelShuffleVector(nodes, 0UL, nodes.size());
            }

            for (size_t i = 0; i < nodes.size(); ++i) {
                PartitionID destination_community;
                const HypernodeID node_to_move = nodes[i];
                const PartitionID source_community = communities[node_to_move];
                const size_t map_size = ratingsFitIntoSmallMap(chg, node_to_move);
                if (!map_size) {
                    //TODO: small Maps for what
                } else {
                    //TODO: large Maps for what 
                }
                if (destination_community != source_community) {
                    makeMove(chg, communities, node_to_move, destination_community);
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
    PartitionID calculateBestMove(ds::CommunityHypergraph& chg, parallel::scalable_vector<HypernodeID>& communities, const HypernodeID v, Map& overlap) {
        //TODO: calculate best move
    }

    KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void makeMove(ds::CommunityHypergraph& chg, parallel::scalable_vector<HypernodeID>& communitties, const HypernodeID node_to_move, const PartitionID destination_community) {
        ASSERT(communitties[node_to_move] != destination_community);
        const HypernodeID source_community = communitties[node_to_move];
        const HyperedgeWeight vol_v = chg.nodeVolume(node_to_move);
        communitties[node_to_move] = destination_community;
        _community_volumes[source_community] -= vol_v;
        _community_volumes[destination_community] += vol_v;
        //TODO: update data structures according to move
    }

    KAHYPAR_ATTRIBUTE_ALWAYS_INLINE size_t ratingsFitIntoSmallMap(const ds::CommunityHypergraph& chg, const HypernodeID v) const {
        // const bool use_vertex_degree_sampling =
        //     _context.coarsening.vertex_degree_sampling_threshold != std::numeric_limits<size_t>::max();
        // const size_t vertex_degree_bounded_rating_map_size = use_vertex_degree_sampling ?
        //     3UL * _context.coarsening.vertex_degree_sampling_threshold : std::numeric_limits<size_t>::max();
        // const size_t cache_efficient_rating_map_size = CacheEfficientRatingMap::MAP_SIZE;
        // const size_t size_of_smaller_rating_map = std::min(
        //     vertex_degree_bounded_rating_map_size, cache_efficient_rating_map_size);

        // // In case the current number of nodes is smaller than size
        // // of the cache-efficient sparse map, the large tmp rating map
        // // consumes less memory
        // if (chg.initialNumNodes() < size_of_smaller_rating_map) {
        //     return chg.initialNumNodes();
        // }

        // // Compute estimation for the upper bound of neighbors of u
        // HypernodeID ub_neighbors_v = 0;
        // for (const HyperedgeID& he : chg.incidentEdges(v)) {
        //     const HypernodeID edge_size = chg.edgeSize(he);
        //     // Ignore large hyperedges
        //     ub_neighbors_v += edge_size;
        //     // If the number of estimated neighbors is greater than the size of the cache efficient rating map / 3, we
        //     // use the large sparse map. The division by 3 also ensures that the fill grade
        //     // of the cache efficient sparse map would be small enough such that linear probing
        //     // is fast.
        //     if (ub_neighbors_v > cache_efficient_rating_map_size / 3UL) {
        //         return std::min(vertex_degree_bounded_rating_map_size, static_cast<size_t>(chg.initialNumNodes()));
        //     }
        // }
        return 0;
    }

    // ! initializes the exit probabilities for each communi
    void initializeExitProbabilities(const ds::CommunityHypergraph& chg) {
        for (size_t i = 0; i < chg.initialNumNodes(); ++i) {
            _community_exit_probability_mul_vol_total[i].store(0);
        }
        for (const HypernodeID hn : chg.nodes()) {
            for (const HyperedgeID he : chg.incidentEdges(hn)) {
                _community_exit_probability_mul_vol_total[hn] += chg.edgeSize(he) * chg.edgeWeight(he);
            }
            _community_exit_probability_mul_vol_total[hn] -= chg.nodeVolume(hn);
            _sum_exit_probability_mul_vol_total += _community_exit_probability_mul_vol_total[hn];
        }
    }

    const size_t _vertex_degree_sampling_threshold;

    // ! contains parameters for the algorithms
    const Context _context;

    // ! reciprocal total volume
    const double _reciprocal_vol_total;

    // ! volumes of each community
    parallel::scalable_vector<AtomicHyperedgeWeight> _community_volumes;

    // ! contains the exit probability of each community multiplied by vol(V) (qi_exit * vol(V))
    parallel::scalable_vector<AtomicPinCOunt> _community_exit_probability_mul_vol_total;

    // ! the sum of all exit probabilities multiplied by vol(V)
    AtomicPinCOunt _sum_exit_probability_mul_vol_total;

    // ! deactivates random node order in local moving
    const bool _deactivate_random;

    FRIEND_TEST(AHyperGraphLocalMovingMapEquation, InitializesTheExitProbabilities);
    FRIEND_TEST(AHyperGraphLocalMovingMapEquation, InitializesTheExitProbabilities0);
    FRIEND_TEST(AHyperGraphLocalMovingMapEquation, InitializesTheExitProbabilities1);
    FRIEND_TEST(AHyperGraphLocalMovingMapEquation, InitializesTheExitProbabilities2);
    FRIEND_TEST(AHyperGraphLocalMovingMapEquation, InitializesTheExitProbabilities3);
    FRIEND_TEST(AHyperGraphLocalMovingMapEquation, InitializesTheExitProbabilities4);
};
}

namespace mt_kahypar::metrics {
double hyp_map_equation(const ds::CommunityHypergraph& chg, const parallel::scalable_vector<HypernodeID>& communities);
}