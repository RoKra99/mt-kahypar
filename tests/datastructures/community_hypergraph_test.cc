#include <mt-kahypar/parallel/tbb_numa_arena.h>

#include "gmock/gmock.h"

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/partition/preprocessing/community_detection/hyp_local_moving_modularity.h"
#include "mt-kahypar/utils/floating_point_comparisons.h"

using ::testing::Test;
using CommunityMove = mt_kahypar::community_detection::CommunityMove;

namespace mt_kahypar {
namespace ds {


template< typename CommunityHG,
    typename HG,
    typename HGFactory>
    struct CommunityHypergraphTypeTraits {
    using CommunityHyperGraph = CommunityHG;
    using Hypergraph = HG;
    using Factory = HGFactory;
};

template<typename TypeTraits>
class ACommunityHypergraph : public Test {

    using CommunityHyperGraph = typename TypeTraits::CommunityHyperGraph;
    using Factory = typename TypeTraits::Factory;

public:
    using Hypergraph = typename TypeTraits::Hypergraph;
    ACommunityHypergraph() :
        hypergraph(Factory::construct(TBBNumaArena::GLOBAL_TASK_GROUP,
                                      7 , 4, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} })),
        community_hypergraph(hypergraph) {
    }

    Hypergraph hypergraph;
    CommunityHypergraph community_hypergraph;
};

// precision with the modularity change is compared (this is the highest precision where all tests pass)
static constexpr Volume PRECISION = 1e-15;

using CommunityHypergraphTestTypes =
::testing::Types<
    CommunityHypergraphTypeTraits<
    CommunityHypergraph,
    StaticHypergraph,
    StaticHypergraphFactory>>;


TYPED_TEST_CASE(ACommunityHypergraph, CommunityHypergraphTestTypes);

TYPED_TEST(ACommunityHypergraph, HasCorrectInitialVolumes) {
    ASSERT_EQ(12, this->community_hypergraph.totalVolume());
    ASSERT_EQ(2, this->community_hypergraph.nodeVolume(0));
    ASSERT_EQ(1, this->community_hypergraph.nodeVolume(1));
    ASSERT_EQ(2, this->community_hypergraph.nodeVolume(2));
    ASSERT_EQ(2, this->community_hypergraph.nodeVolume(3));
    ASSERT_EQ(2, this->community_hypergraph.nodeVolume(4));
    ASSERT_EQ(1, this->community_hypergraph.nodeVolume(5));
    ASSERT_EQ(2, this->community_hypergraph.nodeVolume(6));
    ASSERT_EQ(2, this->community_hypergraph.communityVolume(0));
    ASSERT_EQ(1, this->community_hypergraph.communityVolume(1));
    ASSERT_EQ(2, this->community_hypergraph.communityVolume(2));
    ASSERT_EQ(2, this->community_hypergraph.communityVolume(3));
    ASSERT_EQ(2, this->community_hypergraph.communityVolume(4));
    ASSERT_EQ(1, this->community_hypergraph.communityVolume(5));
    ASSERT_EQ(2, this->community_hypergraph.communityVolume(6));
}

TYPED_TEST(ACommunityHypergraph, SingletonCommunityInitialization) {
    ASSERT_EQ(0, this->community_hypergraph.communityID(0));
    ASSERT_EQ(1, this->community_hypergraph.communityID(1));
    ASSERT_EQ(2, this->community_hypergraph.communityID(2));
    ASSERT_EQ(3, this->community_hypergraph.communityID(3));
    ASSERT_EQ(4, this->community_hypergraph.communityID(4));
    ASSERT_EQ(5, this->community_hypergraph.communityID(5));
    ASSERT_EQ(6, this->community_hypergraph.communityID(6));
}

TYPED_TEST(ACommunityHypergraph, SumOfEdgeWeights) {
    ASSERT_EQ(4, this->community_hypergraph.totalEdgeWeight());
    ASSERT_EQ(0, this->community_hypergraph.dEdgeWeight(0));
    ASSERT_EQ(0, this->community_hypergraph.dEdgeWeight(1));
    ASSERT_EQ(1, this->community_hypergraph.dEdgeWeight(2));
    ASSERT_EQ(2, this->community_hypergraph.dEdgeWeight(3));
    ASSERT_EQ(1, this->community_hypergraph.dEdgeWeight(4));
}

TYPED_TEST(ACommunityHypergraph, HasCorrectInitialCommunityVolumeIterator) {
    const std::vector<HyperedgeWeight> expected_iterator = { 2, 1, 2, 2, 2, 1, 2 };
    size_t pos = 0;
    for (const HyperedgeWeight& hw : this->community_hypergraph.communityVolumes()) {
        ASSERT_EQ(expected_iterator[pos++], hw);
    }
}

TYPED_TEST(ACommunityHypergraph, HasCorrectEdgeSizeIterator) {
    const std::vector<size_t> expected_iterator = { 2, 3, 4 };
    size_t pos = 0;
    for (const PartitionID& d : this->community_hypergraph.edgeSizes()) {
        ASSERT_EQ(expected_iterator[pos++], d);
    }
}


// TODO: Move Tests for HypLocalMovingModularity to it's own file


TYPED_TEST(ACommunityHypergraph, NewModularityDeltaNode0) {
    mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(0);
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
}

TYPED_TEST(ACommunityHypergraph, NewModularityDeltaNode1) {
    mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(1);
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
}

TYPED_TEST(ACommunityHypergraph, NewModularityDeltaNode2) {
    mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(2);
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
}

TYPED_TEST(ACommunityHypergraph, NewModularityDeltaNode3) {
    mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(3);
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
}

TYPED_TEST(ACommunityHypergraph, NewModularityDeltaNode4) {
    mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(4);
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
}

TYPED_TEST(ACommunityHypergraph, NewModularityDeltaNode5) {
    mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(5);
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
}

TYPED_TEST(ACommunityHypergraph, NewModularityDeltaNode6) {
    mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(6);
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    ASSERT_TRUE(mt_kahypar::math::are_almost_equal_ld(move.delta, after - before, PRECISION));
}


TYPED_TEST(ACommunityHypergraph, makeMoveTest) {
    mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
    CommunityMove cm;
    cm.node_to_move = 0;
    cm.delta = -10.0;
    cm.destination_community = 2;
    hlmm.makeMove(cm);
    ASSERT_EQ(0, this->community_hypergraph.communityVolume(0));
    ASSERT_EQ(4, this->community_hypergraph.communityVolume(2));
    ASSERT_EQ(2, this->community_hypergraph.communityID(0));

}
}
}