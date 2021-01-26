#include "gmock/gmock.h"

#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_modularity.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/utils/floating_point_comparisons.h"

using ::testing::Test;

namespace mt_kahypar {
namespace community_detection {

using AHypergraphLocalMoving = ds::HypergraphFixture<ds::StaticHypergraph, ds::StaticHypergraphFactory>;

// precision with which the modularity change is compared (this is the highest precision where all tests pass)
static constexpr Volume PRECISION = 1e-18L;

TEST_F(AHypergraphLocalMoving, CalculatesModularityDelta0) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);
    Volume before = mt_kahypar::metrics::hyp_modularity(community_hypergraph, communities, hlmm);
    CommunityMove move = hlmm.calculateBestMove(community_hypergraph, communities, 0);
    hlmm.makeMove(community_hypergraph, communities, move);
    Volume after = mt_kahypar::metrics::hyp_modularity(community_hypergraph, communities, hlmm);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
    ASSERT_DOUBLE_EQ(move.delta, after - before);
}

TEST_F(AHypergraphLocalMoving, CalculatesModularityDelta1) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);
    Volume before = mt_kahypar::metrics::hyp_modularity(community_hypergraph, communities, hlmm);
    CommunityMove move = hlmm.calculateBestMove(community_hypergraph, communities, 1);
    hlmm.makeMove(community_hypergraph, communities, move);
    Volume after = mt_kahypar::metrics::hyp_modularity(community_hypergraph, communities, hlmm);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
    ASSERT_DOUBLE_EQ(move.delta, after - before);
}

TEST_F(AHypergraphLocalMoving, CalculatesModularityDelta2) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);
    Volume before = mt_kahypar::metrics::hyp_modularity(community_hypergraph, communities, hlmm);
    CommunityMove move = hlmm.calculateBestMove(community_hypergraph, communities, 2);
    hlmm.makeMove(community_hypergraph, communities, move);
    Volume after = mt_kahypar::metrics::hyp_modularity(community_hypergraph, communities, hlmm);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
    ASSERT_DOUBLE_EQ(move.delta, after - before);
}

TEST_F(AHypergraphLocalMoving, CalculatesModularityDelta3) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);
    Volume before = mt_kahypar::metrics::hyp_modularity(community_hypergraph, communities, hlmm);
    CommunityMove move = hlmm.calculateBestMove(community_hypergraph, communities, 3);
    hlmm.makeMove(community_hypergraph, communities, move);
    Volume after = mt_kahypar::metrics::hyp_modularity(community_hypergraph, communities, hlmm);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
    ASSERT_DOUBLE_EQ(move.delta, after - before);
}

TEST_F(AHypergraphLocalMoving, CalculatesModularityDelta4) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);
    Volume before = mt_kahypar::metrics::hyp_modularity(community_hypergraph, communities, hlmm);
    CommunityMove move = hlmm.calculateBestMove(community_hypergraph, communities, 4);
    hlmm.makeMove(community_hypergraph, communities, move);
    Volume after = mt_kahypar::metrics::hyp_modularity(community_hypergraph, communities, hlmm);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
    ASSERT_DOUBLE_EQ(move.delta, after - before);
}

TEST_F(AHypergraphLocalMoving, CalculatesModularityDelta5) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);
    Volume before = mt_kahypar::metrics::hyp_modularity(community_hypergraph, communities, hlmm);
    CommunityMove move = hlmm.calculateBestMove(community_hypergraph, communities, 5);
    hlmm.makeMove(community_hypergraph, communities, move);
    Volume after = mt_kahypar::metrics::hyp_modularity(community_hypergraph, communities, hlmm);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
    ASSERT_DOUBLE_EQ(move.delta, after - before);
}

TEST_F(AHypergraphLocalMoving, CalculatesModularityDelta6) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);
    Volume before = mt_kahypar::metrics::hyp_modularity(community_hypergraph, communities, hlmm);
    CommunityMove move = hlmm.calculateBestMove(community_hypergraph, communities, 6);
    hlmm.makeMove(community_hypergraph, communities, move);
    Volume after = mt_kahypar::metrics::hyp_modularity(community_hypergraph, communities, hlmm);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
    ASSERT_DOUBLE_EQ(move.delta, after - before);
}


TEST_F(AHypergraphLocalMoving, ExecutesAMoveCorrectly) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);
    CommunityMove cm;
    cm.node_to_move = 0;
    cm.delta = -10.0;
    cm.destination_community = 2;
    hlmm.makeMove(community_hypergraph, communities, cm);

    ASSERT_EQ(2, communities[0]);
    const std::vector<HypernodeID> exp_iterator = { 0,1,4,2,2,1,2 };
    size_t pos = 0;
    for (const HypernodeID id : hlmm.communityVolumes(community_hypergraph)) {
        ASSERT_EQ(exp_iterator[pos++], id);
    }
}

// // ############################### Fully Caching Hyperedges ###############################

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunity0) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 0);
    ASSERT_EQ(1, cm.destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunity1) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 1);
    std::set<size_t> possible_best_moves = { 0,3,4 };
    ASSERT_TRUE(possible_best_moves.end() != possible_best_moves.find(cm.destination_community));

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunity2) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 2);
    ASSERT_EQ(5, cm.destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunity3) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 3);
    ASSERT_EQ(4, cm.destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunity4) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 4);
    ASSERT_EQ(3, cm.destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunity5) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 5);
    std::set<size_t> possible_best_moves = { 2,6 };
    ASSERT_TRUE(possible_best_moves.end() != possible_best_moves.find(cm.destination_community));

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunity6) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 6);
    ASSERT_EQ(5, cm.destination_community);

}

// ############################### Partially Cached Hyperedges ###############################

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges0) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 2;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 0);
    ASSERT_EQ(1, cm.destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges1) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 2;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 1);
    std::set<size_t> possible_best_moves = { 0,3,4 };
    ASSERT_TRUE(possible_best_moves.end() != possible_best_moves.find(cm.destination_community));

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges2) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 2;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 2);
    ASSERT_EQ(5, cm.destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges3) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 2;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 3);
    ASSERT_EQ(4, cm.destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges4) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 2;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 4);
    ASSERT_EQ(3, cm.destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges5) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 2;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 5);
    std::set<size_t> possible_best_moves = { 2,6 };
    ASSERT_TRUE(possible_best_moves.end() != possible_best_moves.find(cm.destination_community));

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges6) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 2;
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 6);
    ASSERT_EQ(5, cm.destination_community);

}

// ############################### Without Caching Hyperedges ###############################

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges0) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = std::numeric_limits<size_t>::max();
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 0);
    ASSERT_EQ(1, cm.destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges1) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = std::numeric_limits<size_t>::max();
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 1);
    std::set<size_t> possible_best_moves = { 0,3,4 };
    ASSERT_TRUE(possible_best_moves.end() != possible_best_moves.find(cm.destination_community));

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges2) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = std::numeric_limits<size_t>::max();
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 2);
    ASSERT_EQ(5, cm.destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges3) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = std::numeric_limits<size_t>::max();
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 3);
    ASSERT_EQ(4, cm.destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges4) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = std::numeric_limits<size_t>::max();
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 4);
    ASSERT_EQ(3, cm.destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges5) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = std::numeric_limits<size_t>::max();
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 5);
    std::set<size_t> possible_best_moves = { 2,6 };
    ASSERT_TRUE(possible_best_moves.end() != possible_best_moves.find(cm.destination_community));

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges6) {
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = std::numeric_limits<size_t>::max();
    ds::CommunityHypergraph community_hypergraph(hypergraph, context);
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);

    CommunityMove cm = hlmm.calculateBestMove(community_hypergraph, communities, 6);
    ASSERT_EQ(5, cm.destination_community);

}
}
}