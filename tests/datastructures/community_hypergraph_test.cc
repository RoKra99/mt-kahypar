#include <mt-kahypar/parallel/tbb_numa_arena.h>

#include "gmock/gmock.h"

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_modularity.h"


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

using CommunityHypergraphTestTypes =
::testing::Types<
    CommunityHypergraphTypeTraits<
    CommunityHypergraph,
    StaticHypergraph,
    StaticHypergraphFactory>>;


TYPED_TEST_CASE(ACommunityHypergraph, CommunityHypergraphTestTypes);

TYPED_TEST(ACommunityHypergraph, HasCorrectInitialNodeVolumes) {
    ASSERT_EQ(12, this->community_hypergraph.totalVolume());
    ASSERT_EQ(2, this->community_hypergraph.nodeVolume(0));
    ASSERT_EQ(1, this->community_hypergraph.nodeVolume(1));
    ASSERT_EQ(2, this->community_hypergraph.nodeVolume(2));
    ASSERT_EQ(2, this->community_hypergraph.nodeVolume(3));
    ASSERT_EQ(2, this->community_hypergraph.nodeVolume(4));
    ASSERT_EQ(1, this->community_hypergraph.nodeVolume(5));
    ASSERT_EQ(2, this->community_hypergraph.nodeVolume(6));
}

TYPED_TEST(ACommunityHypergraph, HasCorrectSumOfEdgeWeights) {
    ASSERT_EQ(4, this->community_hypergraph.totalEdgeWeight());
    ASSERT_EQ(0, this->community_hypergraph.edgeWeightBySize(0));
    ASSERT_EQ(0, this->community_hypergraph.edgeWeightBySize(1));
    ASSERT_EQ(1, this->community_hypergraph.edgeWeightBySize(2));
    ASSERT_EQ(2, this->community_hypergraph.edgeWeightBySize(3));
    ASSERT_EQ(1, this->community_hypergraph.edgeWeightBySize(4));
}

TYPED_TEST(ACommunityHypergraph, HasCorrectEdgeSizeIterator) {
    const std::vector<size_t> expected_iterator = { 2, 3, 4 };
    size_t pos = 0;
    for (const PartitionID& d : this->community_hypergraph.edgeSizes()) {
        ASSERT_EQ(expected_iterator[pos++], d);
    }
}

TYPED_TEST(ACommunityHypergraph, ReturnsTheMinimumEdgeSize) {
    ASSERT_EQ(2, this->community_hypergraph.minEdgeSize());
}
}
}