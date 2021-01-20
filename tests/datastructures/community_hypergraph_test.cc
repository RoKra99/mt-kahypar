#include <mt-kahypar/parallel/tbb_numa_arena.h>

#include "gmock/gmock.h"

#include "tests/datastructures/hypergraph_fixtures.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_modularity.h"
#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_louvain.h"


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
        community_hypergraph(hypergraph) {}

    void verifyIncidentNets(const Hypergraph& hg,
        const HypernodeID hn,
        const std::set<HypernodeID>& reference) {
        size_t count = 0;
        for (const HyperedgeID& he : hg.incidentEdges(hn)) {
            ASSERT_TRUE(reference.find(he) != reference.end()) << V(he);
            count++;
        }
        ASSERT_EQ(count, reference.size());
    }

    void verifyPins(const Hypergraph& hg,
        const std::vector<HyperedgeID> hyperedges,
        const std::vector< std::multiset<HypernodeID> >& references,
        bool log = false) {
        ASSERT(hyperedges.size() == references.size());
        for (size_t i = 0; i < hyperedges.size(); ++i) {
            const HyperedgeID he = hyperedges[i];
            const std::multiset<HypernodeID>& reference = references[i];
            size_t count = 0;
            for (const HypernodeID& pin : hg.pins(he)) {
                if (log) LOG << V(he) << V(pin);
                ASSERT_TRUE(reference.find(pin) != reference.end()) << V(he) << V(pin);
                count++;
            }
            ASSERT_EQ(count, reference.size());
        }
    }

    void verifyCommunityCounts(const CommunityHypergraph& chg, const HyperedgeID he, std::vector<std::set<HypernodeID>> exp_communities) {
        size_t count = 0;
        for (auto& e : chg.singleCuts(he)) {
            ASSERT_TRUE(exp_communities[0].find(e) != exp_communities[0].end());
            ++count;
        }
        ASSERT_EQ(count, exp_communities[0].size());
        count = 0;
        for (auto& e : chg.multiCuts(he)) {
            ASSERT_TRUE(exp_communities[1].find(e.first) != exp_communities[1].end());
            ++count;
        }
        ASSERT_EQ(count, exp_communities[1].size());
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

TYPED_TEST(ACommunityHypergraph, HasCorrectCommunityCounts) {
    this->verifyCommunityCounts(this->community_hypergraph, 0, { {0, 2},{} });
    this->verifyCommunityCounts(this->community_hypergraph, 1, { {0, 1, 3, 4},{} });
    this->verifyCommunityCounts(this->community_hypergraph, 2, { {3, 4, 6},{} });
    this->verifyCommunityCounts(this->community_hypergraph, 3, { {2, 5, 6},{} });
}

TYPED_TEST(ACommunityHypergraph, ContractsCommunities1) {
    parallel::scalable_vector<HypernodeID> c_communities = { 1,4,1,5,5,4,5 };
    StaticHypergraph hg;
    CommunityHypergraph cchg = this->community_hypergraph.contract(hg, c_communities);

    // community mapping
    ASSERT_EQ(0, c_communities[0]);
    ASSERT_EQ(1, c_communities[1]);
    ASSERT_EQ(0, c_communities[2]);
    ASSERT_EQ(2, c_communities[3]);
    ASSERT_EQ(2, c_communities[4]);
    ASSERT_EQ(1, c_communities[5]);
    ASSERT_EQ(2, c_communities[6]);

    // stats
    ASSERT_EQ(3, cchg.initialNumNodes());
    ASSERT_EQ(4, cchg.initialNumEdges());
    ASSERT_EQ(12, cchg.totalVolume());
    ASSERT_EQ(4, cchg.totalEdgeWeight());
    ASSERT_EQ(12, hg.initialNumPins());

    // edgeWeight by size
    ASSERT_EQ(1, cchg.edgeWeightBySize(2));
    ASSERT_EQ(2, cchg.edgeWeightBySize(3));
    ASSERT_EQ(1, cchg.edgeWeightBySize(4));

    // valid edge sizes
    std::vector<size_t> expected_valid_edge_sizes = { 2,3,4 };
    size_t i = 0;
    for (const size_t d : cchg.edgeSizes()) {
        ASSERT_EQ(d, expected_valid_edge_sizes[i++]);
    }

    // verify edge sizes
    ASSERT_EQ(2, cchg.edgeSize(0));
    ASSERT_EQ(4, cchg.edgeSize(1));
    ASSERT_EQ(3, cchg.edgeSize(2));
    ASSERT_EQ(3, cchg.edgeSize(3));

    // verify edge weights
    ASSERT_EQ(1, cchg.edgeWeight(0));
    ASSERT_EQ(1, cchg.edgeWeight(1));
    ASSERT_EQ(1, cchg.edgeWeight(2));
    ASSERT_EQ(1, cchg.edgeWeight(3));

    // verify community counts
    this->verifyCommunityCounts(cchg, 0, { {0}, {} });
    this->verifyCommunityCounts(cchg, 1, { {0,1,2}, {} });
    this->verifyCommunityCounts(cchg, 2, { {2}, {} });
    this->verifyCommunityCounts(cchg, 3, { {0,1,2}, {} });

    // verify node volumes
    ASSERT_EQ(4, cchg.nodeVolume(0));
    ASSERT_EQ(2, cchg.nodeVolume(1));
    ASSERT_EQ(6, cchg.nodeVolume(2));

    // verify hypergraph structure
    this->verifyIncidentNets(hg, 0, { 0, 1, 3 });
    this->verifyIncidentNets(hg, 1, { 1, 3 });
    this->verifyIncidentNets(hg, 2, { 1, 2, 3 });
    this->verifyPins(hg, { 0,1,2,3 }, { {0,0}, {0,1,2,2}, {2,2,2}, {0,1,2} });
}

TYPED_TEST(ACommunityHypergraph, ContractsCommunities2) {
    parallel::scalable_vector<HypernodeID> c_communities = { 1,1,5,4,4,5,5 };
    StaticHypergraph hg;
    CommunityHypergraph cchg = this->community_hypergraph.contract(hg, c_communities);

    // community mapping
    ASSERT_EQ(0, c_communities[0]);
    ASSERT_EQ(0, c_communities[1]);
    ASSERT_EQ(2, c_communities[2]);
    ASSERT_EQ(1, c_communities[3]);
    ASSERT_EQ(1, c_communities[4]);
    ASSERT_EQ(2, c_communities[5]);
    ASSERT_EQ(2, c_communities[6]);

    // stats
    ASSERT_EQ(3, cchg.initialNumNodes());
    ASSERT_EQ(4, cchg.initialNumEdges());
    ASSERT_EQ(12, cchg.totalVolume());
    ASSERT_EQ(4, cchg.totalEdgeWeight());
    ASSERT_EQ(12, hg.initialNumPins());

    // edgeWeight by size
    ASSERT_EQ(1, cchg.edgeWeightBySize(2));
    ASSERT_EQ(2, cchg.edgeWeightBySize(3));
    ASSERT_EQ(1, cchg.edgeWeightBySize(4));


    // valid edge sizes
    std::vector<size_t> expected_valid_edge_sizes = { 2,3,4 };
    size_t i = 0;
    for (const size_t d : cchg.edgeSizes()) {
        ASSERT_EQ(d, expected_valid_edge_sizes[i++]);
    }

    // verify edge sizes
    ASSERT_EQ(2, cchg.edgeSize(0));
    ASSERT_EQ(4, cchg.edgeSize(1));
    ASSERT_EQ(3, cchg.edgeSize(2));
    ASSERT_EQ(3, cchg.edgeSize(3));

    // verify edge weights
    ASSERT_EQ(1, cchg.edgeWeight(0));
    ASSERT_EQ(1, cchg.edgeWeight(1));
    ASSERT_EQ(1, cchg.edgeWeight(2));
    ASSERT_EQ(1, cchg.edgeWeight(3));

    // verify community counts
    this->verifyCommunityCounts(cchg, 0, { {0, 2}, {} });
    this->verifyCommunityCounts(cchg, 1, { {0,1}, {} });
    this->verifyCommunityCounts(cchg, 2, { {1,2}, {} });
    this->verifyCommunityCounts(cchg, 3, { {2}, {} });

    // verify node volumes
    ASSERT_EQ(3, cchg.nodeVolume(0));
    ASSERT_EQ(4, cchg.nodeVolume(1));
    ASSERT_EQ(5, cchg.nodeVolume(2));

    // verify hypergraph structure
    this->verifyIncidentNets(hg, 0, { 0, 1 });
    this->verifyIncidentNets(hg, 1, { 1, 2 });
    this->verifyIncidentNets(hg, 2, { 0, 2, 3 });
    this->verifyPins(hg, { 0,1,2,3 }, { {0,2}, {0,0,1,1}, {1,1,2}, {2,2,2} });
}

TYPED_TEST(ACommunityHypergraph, ContractsCommunities3) {
    parallel::scalable_vector<HypernodeID> c_communities = { 1,1,1,1,4,4,5 };
    StaticHypergraph hg;
    CommunityHypergraph cchg = this->community_hypergraph.contract(hg, c_communities);

    // community mapping
    ASSERT_EQ(0, c_communities[0]);
    ASSERT_EQ(0, c_communities[1]);
    ASSERT_EQ(0, c_communities[2]);
    ASSERT_EQ(0, c_communities[3]);
    ASSERT_EQ(1, c_communities[4]);
    ASSERT_EQ(1, c_communities[5]);
    ASSERT_EQ(2, c_communities[6]);

    // stats
    ASSERT_EQ(3, cchg.initialNumNodes());
    ASSERT_EQ(3, cchg.initialNumEdges());
    ASSERT_EQ(12, cchg.totalVolume());
    ASSERT_EQ(4, cchg.totalEdgeWeight());
    ASSERT_EQ(9, hg.initialNumPins());

    // edgeWeight by size
    ASSERT_EQ(1, cchg.edgeWeightBySize(2));
    ASSERT_EQ(2, cchg.edgeWeightBySize(3));
    ASSERT_EQ(1, cchg.edgeWeightBySize(4));

    // valid edge sizes
    std::vector<size_t> expected_valid_edge_sizes = { 2,3,4 };
    size_t i = 0;
    for (const size_t d : cchg.edgeSizes()) {
        ASSERT_EQ(d, expected_valid_edge_sizes[i++]);
    }

    // verify edge sizes
    ASSERT_EQ(2, cchg.edgeSize(0));
    ASSERT_EQ(4, cchg.edgeSize(1));
    ASSERT_EQ(3, cchg.edgeSize(2));

    // verify edge weights
    ASSERT_EQ(1, cchg.edgeWeight(0));
    ASSERT_EQ(1, cchg.edgeWeight(1));
    ASSERT_EQ(2, cchg.edgeWeight(2));

    // verify community counts
    this->verifyCommunityCounts(cchg, 0, { {0}, {} });
    this->verifyCommunityCounts(cchg, 1, { {0,1}, {} });
    this->verifyCommunityCounts(cchg, 2, { {0,1,2}, {} });

    // verify node volumes
    ASSERT_EQ(7, cchg.nodeVolume(0));
    ASSERT_EQ(3, cchg.nodeVolume(1));
    ASSERT_EQ(2, cchg.nodeVolume(2));

    // verify hypergraph structure
    this->verifyIncidentNets(hg, 0, { 0, 1, 2 });
    this->verifyIncidentNets(hg, 1, { 1, 2 });
    this->verifyIncidentNets(hg, 2, { 2 });
    this->verifyPins(hg, { 0,1,2 }, { {0,0}, {0,0,0,1}, {0,1,2} });
}

TYPED_TEST(ACommunityHypergraph, ContractsCommunities4) {
    parallel::scalable_vector<HypernodeID> c_communities = { 0,0,0,0,0,0,0 };
    StaticHypergraph hg;
    CommunityHypergraph cchg = this->community_hypergraph.contract(hg, c_communities);

    // community mapping
    ASSERT_EQ(0, c_communities[0]);
    ASSERT_EQ(0, c_communities[1]);
    ASSERT_EQ(0, c_communities[2]);
    ASSERT_EQ(0, c_communities[3]);
    ASSERT_EQ(0, c_communities[4]);
    ASSERT_EQ(0, c_communities[5]);
    ASSERT_EQ(0, c_communities[6]);

    // stats
    ASSERT_EQ(1, cchg.initialNumNodes());
    ASSERT_EQ(3, cchg.initialNumEdges());
    ASSERT_EQ(12, cchg.totalVolume());
    ASSERT_EQ(4, cchg.totalEdgeWeight());
    ASSERT_EQ(9, hg.initialNumPins());

    // edgeWeight by size
    ASSERT_EQ(1, cchg.edgeWeightBySize(2));
    ASSERT_EQ(2, cchg.edgeWeightBySize(3));
    ASSERT_EQ(1, cchg.edgeWeightBySize(4));

    // valid edge sizes
    std::vector<size_t> expected_valid_edge_sizes = { 2,3,4 };
    size_t i = 0;
    for (const size_t d : cchg.edgeSizes()) {
        ASSERT_EQ(d, expected_valid_edge_sizes[i++]);
    }

    // verify edge sizes
    ASSERT_EQ(2, cchg.edgeSize(0));
    ASSERT_EQ(4, cchg.edgeSize(1));
    ASSERT_EQ(3, cchg.edgeSize(2));

    // verify edge weights
    ASSERT_EQ(1, cchg.edgeWeight(0));
    ASSERT_EQ(1, cchg.edgeWeight(1));
    ASSERT_EQ(2, cchg.edgeWeight(2));

    // verify community counts
    this->verifyCommunityCounts(cchg, 0, { {0}, {} });
    this->verifyCommunityCounts(cchg, 1, { {0}, {} });
    this->verifyCommunityCounts(cchg, 2, { {0}, {} });

    // verify node volumes
    ASSERT_EQ(12, cchg.nodeVolume(0));

    // verify hypergraph structure
    this->verifyIncidentNets(hg, 0, { 0, 1, 2 });
    this->verifyPins(hg, { 0,1,2 }, { {0,0}, {0,0,0,0}, {0,0,0} });
}
}
}