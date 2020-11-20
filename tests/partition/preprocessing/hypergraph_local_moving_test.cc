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

// precision with the modularity change is compared (this is the highest precision where all tests pass)
static constexpr Volume PRECISION = 1e-15L;

TEST_F(AHypergraphLocalMoving, NewModularityDeltaNode0) {
    ds::CommunityHypergraph community_hypergraph(hypergraph);
    HypergraphLocalMovingModularity hlmm(community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(0);
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(community_hypergraph);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
}

TEST_F(AHypergraphLocalMoving, NewModularityDeltaNode1) {
    ds::CommunityHypergraph community_hypergraph(hypergraph);
    HypergraphLocalMovingModularity hlmm(community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(1);
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(community_hypergraph);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
}

TEST_F(AHypergraphLocalMoving, NewModularityDeltaNode2) {
    ds::CommunityHypergraph community_hypergraph(hypergraph);
    HypergraphLocalMovingModularity hlmm(community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(2);
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(community_hypergraph);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
}

TEST_F(AHypergraphLocalMoving, NewModularityDeltaNode3) {
    ds::CommunityHypergraph community_hypergraph(hypergraph);
    HypergraphLocalMovingModularity hlmm(community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(3);
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(community_hypergraph);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
}

TEST_F(AHypergraphLocalMoving, NewModularityDeltaNode4) {
    ds::CommunityHypergraph community_hypergraph(hypergraph);
    HypergraphLocalMovingModularity hlmm(community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(4);
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(community_hypergraph);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
}

TEST_F(AHypergraphLocalMoving, NewModularityDeltaNode5) {
    ds::CommunityHypergraph community_hypergraph(hypergraph);
    HypergraphLocalMovingModularity hlmm(community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(5);
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(community_hypergraph);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
}

TEST_F(AHypergraphLocalMoving, NewModularityDeltaNode6) {
    ds::CommunityHypergraph community_hypergraph(hypergraph);
    HypergraphLocalMovingModularity hlmm(community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(6);
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(community_hypergraph);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
}


TEST_F(AHypergraphLocalMoving, makeMoveTest) {
    ds::CommunityHypergraph community_hypergraph(hypergraph);
    HypergraphLocalMovingModularity hlmm(community_hypergraph);
    CommunityMove cm;
    cm.node_to_move = 0;
    cm.delta = -10.0;
    cm.destination_community = 2;
    hlmm.makeMove(cm);
    ASSERT_EQ(0, community_hypergraph.communityVolume(0));
    ASSERT_EQ(4, community_hypergraph.communityVolume(2));
    ASSERT_EQ(2, community_hypergraph.communityID(0));

}
}
}