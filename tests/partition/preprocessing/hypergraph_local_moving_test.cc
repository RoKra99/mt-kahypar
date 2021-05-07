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

class AHypergraphLocalMoving : public ds::HypergraphFixture<Hypergraph, HypergraphFactory> {
    using Base = ds::HypergraphFixture<Hypergraph, HypergraphFactory>;

public:
    AHypergraphLocalMoving() :
        Base(),
        chg(nullptr),
        context() {
        context.preprocessing.community_detection.max_pass_iterations = 100;
        context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
        context.shared_memory.num_threads = 1;
        context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
        context.preprocessing.community_detection.tie_breaking_rule = TieBreakingRule::random;

        chg = std::make_unique<ds::CommunityHypergraph>(hypergraph, context);
    }

    using Base::hypergraph;
    std::unique_ptr<ds::CommunityHypergraph> chg;
    Context context;
};

TEST_F(AHypergraphLocalMoving, ExecutesAMoveCorrectly) {
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);
    hlmm.makeMove(*chg, communities, 0, 2);

    ASSERT_EQ(2, communities[0]);
    const std::vector<HypernodeID> exp_iterator = { 0,1,4,2,2,1,2 };
    size_t pos = 0;
    for (const HypernodeID id : hlmm.communityVolumes(*chg)) {
        ASSERT_EQ(exp_iterator[pos++], id);
    }
}

// // ############################### Fully Caching Hyperedges ###############################

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunity0) {
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 0, hlmm._large_edge_contribution_map.local());
    ASSERT_EQ(1, destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunity1) {
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 1, hlmm._large_edge_contribution_map.local());
    std::set<size_t> possible_best_moves = { 0,3,4 };
    ASSERT_TRUE(possible_best_moves.end() != possible_best_moves.find(destination_community));

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunity2) {
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 2, hlmm._large_edge_contribution_map.local());
    ASSERT_EQ(5, destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunity3) {
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 3, hlmm._large_edge_contribution_map.local());
    ASSERT_EQ(4, destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunity4) {
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 4, hlmm._large_edge_contribution_map.local());
    ASSERT_EQ(3, destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunity5) {
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 5, hlmm._large_edge_contribution_map.local());
    std::set<size_t> possible_best_moves = { 2,6 };
    ASSERT_TRUE(possible_best_moves.end() != possible_best_moves.find(destination_community));

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunity6) {
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 6, hlmm._large_edge_contribution_map.local());
    ASSERT_EQ(5, destination_community);

}

// ############################### Partially Cached Hyperedges ###############################

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges0) {
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 2;
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 0, hlmm._large_edge_contribution_map.local());
    ASSERT_EQ(1, destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges1) {
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 2;
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 1, hlmm._large_edge_contribution_map.local());
    std::set<size_t> possible_best_moves = { 0,3,4 };
    ASSERT_TRUE(possible_best_moves.end() != possible_best_moves.find(destination_community));

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges2) {
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 2;
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 2, hlmm._large_edge_contribution_map.local());
    ASSERT_EQ(5, destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges3) {
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 2;
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 3, hlmm._large_edge_contribution_map.local());
    ASSERT_EQ(4, destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges4) {
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 2;
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 4, hlmm._large_edge_contribution_map.local());
    ASSERT_EQ(3, destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges5) {
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 2;
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 5, hlmm._large_edge_contribution_map.local());
    std::set<size_t> possible_best_moves = { 2,6 };
    ASSERT_TRUE(possible_best_moves.end() != possible_best_moves.find(destination_community));

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithPartiallyCachedHyperedges6) {
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 2;
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 6, hlmm._large_edge_contribution_map.local());
    ASSERT_EQ(5, destination_community);

}

// ############################### Without Caching Hyperedges ###############################

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges0) {
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = std::numeric_limits<size_t>::max();
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 0, hlmm._large_edge_contribution_map.local());
    ASSERT_EQ(1, destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges1) {
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = std::numeric_limits<size_t>::max();
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 1, hlmm._large_edge_contribution_map.local());
    std::set<size_t> possible_best_moves = { 0,3,4 };
    ASSERT_TRUE(possible_best_moves.end() != possible_best_moves.find(destination_community));

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges2) {
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = std::numeric_limits<size_t>::max();
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 2, hlmm._large_edge_contribution_map.local());
    ASSERT_EQ(5, destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges3) {
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = std::numeric_limits<size_t>::max();
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 3, hlmm._large_edge_contribution_map.local());
    ASSERT_EQ(4, destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges4) {
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = std::numeric_limits<size_t>::max();
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 4, hlmm._large_edge_contribution_map.local());
    ASSERT_EQ(3, destination_community);

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges5) {
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = std::numeric_limits<size_t>::max();
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 5, hlmm._large_edge_contribution_map.local());
    std::set<size_t> possible_best_moves = { 2,6 };
    ASSERT_TRUE(possible_best_moves.end() != possible_best_moves.find(destination_community));

}

TEST_F(AHypergraphLocalMoving, CalulatesBestDestinationCommunityWithoutCachedHyperedges6) {
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = std::numeric_limits<size_t>::max();
    HypergraphLocalMovingModularity hlmm(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmm.initializeCommunityVolumes(*chg, communities);

    PartitionID destination_community = hlmm.calculateBestMove(*chg, communities, 6, hlmm._large_edge_contribution_map.local());
    ASSERT_EQ(5, destination_community);

}
}
}