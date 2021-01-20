#include "gmock/gmock.h"

#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_modularity.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/utils/floating_point_comparisons.h"

using ::testing::Test;

namespace mt_kahypar {
namespace community_detection {

using AHypergraphLocalMoving = ds::HypergraphFixture;

// precision with which the modularity change is compared (this is the highest precision where all tests pass)
static constexpr Volume PRECISION = 1e-18L;

TEST_F(AHypergraphLocalMoving, CalculatesModularityDelta0) {
    ds::CommunityHypergraph community_hypergraph(hypergraph);
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
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
    ds::CommunityHypergraph community_hypergraph(hypergraph);
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
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
    ds::CommunityHypergraph community_hypergraph(hypergraph);
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
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
    ds::CommunityHypergraph community_hypergraph(hypergraph);
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
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
    ds::CommunityHypergraph community_hypergraph(hypergraph);
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
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
    ds::CommunityHypergraph community_hypergraph(hypergraph);
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
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
    ds::CommunityHypergraph community_hypergraph(hypergraph);
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
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
    ds::CommunityHypergraph community_hypergraph(hypergraph);
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    HypergraphLocalMovingModularity hlmm(community_hypergraph, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(community_hypergraph, communities);
    CommunityMove cm;
    cm.node_to_move = 0;
    cm.delta = -10.0;
    cm.destination_community = 2;
    hlmm.makeMove(community_hypergraph, communities, cm);

    ASSERT_EQ(2, communities[0]);
    const std::vector<HypernodeID> exp_iterator = {0,1,4,2,2,1,2};
    size_t pos = 0;
    for (const HypernodeID id : hlmm.communityVolumes(community_hypergraph)) {
        ASSERT_EQ(exp_iterator[pos++], id);
    }
}
}
}