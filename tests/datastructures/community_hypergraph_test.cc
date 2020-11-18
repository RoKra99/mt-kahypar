#include <mt-kahypar/parallel/tbb_numa_arena.h>

#include "gmock/gmock.h"

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/datastructures/community_hypergraph.h"
#include "mt-kahypar/partition/preprocessing/community_detection/hyp_local_moving_modularity.h"

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
                                      7, 4, { {0, 2}, {0, 1, 3, 4}, {3, 4, 6}, {2, 5, 6} })),
        community_hypergraph(hypergraph) {
    }

    Hypergraph hypergraph;
    CommunityHypergraph community_hypergraph;
};

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

TYPED_TEST(ACommunityHypergraph, HasCorrectEdgeSizeIterator) {
    const std::vector<size_t> expected_iterator = {2, 3, 4};
    size_t pos = 0;
    for (const PartitionID& d : this->community_hypergraph.edgeSizes()) {
        ASSERT_EQ(expected_iterator[pos++], d);
    }
}


// TODO: Move Tests for HypLocalMovingModularity to it's own file
// TYPED_TEST(ACommunityHypergraph, ModularityDeltaNode0) {
//     mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
//     Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
//     Volume delta = hlmm.tryMove(0);
//     Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
//     ASSERT_DOUBLE_EQ(delta, after - before);
// }

// TYPED_TEST(ACommunityHypergraph, ModularityDeltaNode1) {
//     mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
//     Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
//     Volume delta = hlmm.tryMove(1);
//     Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
//     ASSERT_DOUBLE_EQ(delta, after - before);
// }

// TYPED_TEST(ACommunityHypergraph, ModularityDeltaNode2) {
//     mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
//     Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
//     Volume delta = hlmm.tryMove(2);
//     Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
//     ASSERT_DOUBLE_EQ(delta, after - before);
// }

// TYPED_TEST(ACommunityHypergraph, ModularityDeltaNode3) {
//     mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
//     Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
//     Volume delta = hlmm.tryMove(3);
//     Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
//     ASSERT_DOUBLE_EQ(delta, after - before);
// }

// TYPED_TEST(ACommunityHypergraph, ModularityDeltaNode4) {
//     mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
//     Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
//     Volume delta = hlmm.tryMove(4);
//     Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
//     ASSERT_DOUBLE_EQ(delta, after - before);
// }

// TYPED_TEST(ACommunityHypergraph, ModularityDeltaNode5) {
//     mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
//     Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
//     Volume delta = hlmm.tryMove(5);
//     Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
//     ASSERT_DOUBLE_EQ(delta, after - before);
// }

// TYPED_TEST(ACommunityHypergraph, ModularityDeltaNode6) {
//     mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
//     Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
//     Volume delta = hlmm.tryMove(6);
//     Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
//     ASSERT_DOUBLE_EQ(delta, after - before);
// }


TYPED_TEST(ACommunityHypergraph, NewModularityDeltaNode0) {
    mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(0);
    Volume delta = move.delta;
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    ASSERT_DOUBLE_EQ(delta, after - before);
}

TYPED_TEST(ACommunityHypergraph, NewModularityDeltaNode1) {
    mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(1);
    Volume delta = move.delta;
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    ASSERT_DOUBLE_EQ(delta, after - before);
}

TYPED_TEST(ACommunityHypergraph, NewModularityDeltaNode2) {
    mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(2);
    Volume delta = move.delta;
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    ASSERT_DOUBLE_EQ(delta, after - before);
}

TYPED_TEST(ACommunityHypergraph, NewModularityDeltaNode3) {
    mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(3);
    Volume delta = move.delta;
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    ASSERT_DOUBLE_EQ(delta, after - before);
}

TYPED_TEST(ACommunityHypergraph, NewModularityDeltaNode4) {
    mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(4);
    Volume delta = move.delta;
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    ASSERT_DOUBLE_EQ(delta, after - before);
}

TYPED_TEST(ACommunityHypergraph, NewModularityDeltaNode5) {
    mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(5);
    Volume delta = move.delta;
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    ASSERT_DOUBLE_EQ(delta, after - before);
}

TYPED_TEST(ACommunityHypergraph, NewModularityDeltaNode6) {
    mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
    Volume before = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    CommunityMove move = hlmm.calculateBestMove(6);
    Volume delta = move.delta;
    hlmm.makeMove(move);
    Volume after = mt_kahypar::metrics::hyp_modularity(this->community_hypergraph);
    ASSERT_DOUBLE_EQ(delta, after - before);
}


TYPED_TEST(ACommunityHypergraph, InitializingWeightToVector) {
    mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
    ASSERT_DOUBLE_EQ(std::numeric_limits<Volume>::max(), hlmm.weightTOCommunity(0));
    ASSERT_DOUBLE_EQ(std::numeric_limits<Volume>::max(), hlmm.weightTOCommunity(1));
    ASSERT_DOUBLE_EQ(std::numeric_limits<Volume>::max(), hlmm.weightTOCommunity(2));
    ASSERT_DOUBLE_EQ(std::numeric_limits<Volume>::max(), hlmm.weightTOCommunity(3));
    ASSERT_DOUBLE_EQ(std::numeric_limits<Volume>::max(), hlmm.weightTOCommunity(4));
    ASSERT_DOUBLE_EQ(std::numeric_limits<Volume>::max(), hlmm.weightTOCommunity(5));
    ASSERT_DOUBLE_EQ(std::numeric_limits<Volume>::max(), hlmm.weightTOCommunity(6));
}

TYPED_TEST(ACommunityHypergraph, makeMoveTest) {
    mt_kahypar::community_detection::HypLocalMovingModularity hlmm(this->community_hypergraph);
    HyperedgeWeight manuel = this->community_hypergraph.communityVolume(2) + this->community_hypergraph.nodeVolume(0);
    CommunityMove cm;
    cm.node_to_move = 0;
    cm.delta = -10.0;
    cm.destination_community = 2;
    hlmm.makeMove(cm);
    ASSERT_EQ(0,this->community_hypergraph.communityVolume(0));
    ASSERT_EQ(manuel, this->community_hypergraph.communityVolume(2));
    ASSERT_EQ(2, this->community_hypergraph.communityID(0));

}
}
}