#include "gmock/gmock.h"

#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/datastructures/community_hypergraph.h"
#include "mt-kahypar/partition/context.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {

using ACommunityCount = ds::HypergraphFixture<StaticHypergraph, StaticHypergraphFactory>;

TEST_F(ACommunityCount, InitializesEdge0Properly) {
    Context context;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    CommunityHypergraph chg(hypergraph, context);
    const std::vector<size_t>  expected_iterator = { 0,2 };
    int pos = 0;
    for (const auto& e : chg.singleCuts(0)) {
        ASSERT_EQ(expected_iterator[pos++], e);
    }
    for (const auto& e : chg.multiCuts(0)) {
        LOG << e.first << e.second;
        ASSERT_TRUE(false);
    }
}

TEST_F(ACommunityCount, InitializesEdge1Properly) {
    Context context;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    CommunityHypergraph chg(hypergraph, context);
    const std::vector<size_t>  expected_iterator = { 0, 1, 3, 4 };
    int pos = 0;
    for (const auto& e : chg.singleCuts(1)) {
        ASSERT_EQ(expected_iterator[pos++], e);
    }
    for (const auto& e : chg.multiCuts(1)) {
        LOG << e.first << e.second;
        ASSERT_TRUE(false);
    }
}

TEST_F(ACommunityCount, InitializesEdge2Properly) {
    Context context;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    CommunityHypergraph chg(hypergraph, context);
    const std::vector<size_t>  expected_iterator = { 3,4,6 };
    int pos = 0;
    for (const auto& e : chg.singleCuts(2)) {
        ASSERT_EQ(expected_iterator[pos++], e);
    }
    for (const auto& e : chg.multiCuts(2)) {
        LOG << e.first << e.second;
        ASSERT_TRUE(false);
    }
}

TEST_F(ACommunityCount, InitializesEdge3Properly) {
    Context context;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    CommunityHypergraph chg(hypergraph, context);
    const std::vector<size_t>  expected_iterator = { 2,5,6 };
    int pos = 0;
    for (const auto& e : chg.singleCuts(3)) {
        ASSERT_EQ(expected_iterator[pos++], e);
    }
    for (const auto& e : chg.multiCuts(3)) {
        LOG << e.first << e.second;
        ASSERT_TRUE(false);
    }
}

TEST_F(ACommunityCount, AddsANewCommunity) {
    Context context;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    CommunityHypergraph chg(hypergraph, context);
    // we first a hvae to remove a community since we can't add more communities than the edgesize
    chg.removeCommunityFromHyperedge(0,2);
    chg.addCommunityToHyperedge(0, 1);
    const std::vector<size_t>  expected_iterator = { 0,1 };
    int pos = 0;
    for (const auto& e : chg.singleCuts(0)) {
        ASSERT_EQ(expected_iterator[pos++], e);
    }
}

TEST_F(ACommunityCount, MovesAcommunityToTheHashmap) {
    Context context;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    CommunityHypergraph chg(hypergraph, context);
    chg.addCommunityToHyperedge(0, 2);
    const std::vector<size_t>  expected_iterator = { 0 };
    const std::vector<size_t>  expected_iterator2 = { 2 };
    const std::vector<size_t>  expected_value2 = { 2 };
    int pos = 0;
    for (const auto& e : chg.singleCuts(0)) {
        ASSERT_EQ(expected_iterator[pos++], e);
    }
    pos = 0;
    for (const auto& e : chg.multiCuts(0)) {
        ASSERT_EQ(expected_iterator2[pos], e.first);
        ASSERT_EQ(expected_value2[pos], e.second);
        ++pos;
    }
}

TEST_F(ACommunityCount, RemovesACommunity) {
    Context context;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    CommunityHypergraph chg(hypergraph, context);
    chg.removeCommunityFromHyperedge(0, 0);
    const std::vector<size_t>  expected_iterator = { 2 };
    int pos = 0;
    for (const auto& e : chg.singleCuts(0)) {
        ASSERT_EQ(expected_iterator[pos++], e);
    }
    for (const auto& e : chg.multiCuts(0)) {
        LOG << e.first << e.second;
        ASSERT_TRUE(false);
    }
}

TEST_F(ACommunityCount, RemovesAcommunityFromTheHashmap) {
    Context context;
    context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
    CommunityHypergraph chg(hypergraph, context);
    chg.addCommunityToHyperedge(0, 0);
    chg.removeCommunityFromHyperedge(0, 0);
    chg.removeCommunityFromHyperedge(0, 0);
    const std::vector<size_t>  expected_iterator = { 2 };
    int pos = 0;
    for (const auto& e : chg.singleCuts(0)) {
        ASSERT_EQ(expected_iterator[pos++], e);
    }
    for (const auto& e : chg.multiCuts(0)) {
        LOG << e.first << e.second;
        ASSERT_TRUE(false);
    }
}

}
}