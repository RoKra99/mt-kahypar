/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include "gmock/gmock.h"

#include "tests/datastructures/hypergraph_fixtures.h"

using ::testing::Test;

namespace mt_kahypar {
namespace ds {
using AHypergraphWithTwoStreamingHypergraphs = AHypergraph<2>;
using TestHypergraph = typename AHypergraphWithTwoStreamingHypergraphs::TestHypergraph;
using TestStreamingHypergraph = typename AHypergraphWithTwoStreamingHypergraphs::TestStreamingHypergraph;
using TBBArena = typename AHypergraphWithTwoStreamingHypergraphs::TBBArena;

TestHypergraph construct_test_hypergraph(const AHypergraphWithTwoStreamingHypergraphs& test) {
  return test.construct_hypergraph(7, { { 0, 2 }, { 0, 1, 3, 4 }, { 3, 4, 6 }, { 2, 5, 6 } },
                                   { 0, 0, 0, 1, 1, 1, 1 },
                                   { 0, 0, 1, 1 });
}

void assignPartitionIDs(TestHypergraph& hypergraph) {
  for (const HypernodeID& hn : hypergraph.nodes()) {
    PartitionID part_id = TestStreamingHypergraph::get_numa_node_of_vertex(hn);
    hypergraph.setNodePart(hn, part_id);
  }
  hypergraph.updateGlobalPartInfos();
  hypergraph.initializeNumCutHyperedges();
}

template <typename IDType>
auto identity = [](const IDType& id) {
                  return id;
                };

template <typename IDType, typename F, typename K = decltype(identity<IDType>)>
void verifyIterator(const std::set<IDType>& reference, F&& it_func, K map_func = identity<IDType>, bool log = false) {
  size_t count = 0;
  for (const IDType& id : it_func()) {
    if (log) LOG << V(id) << V(map_func(id));
    ASSERT_TRUE(reference.find(map_func(id)) != reference.end()) << V(map_func(id));
    count++;
  }
  ASSERT_EQ(count, reference.size());
}

void verifyPinIterators(const TestHypergraph& hypergraph,
                        const std::vector<HyperedgeID> hyperedges,
                        const std::vector<std::set<HypernodeID> >& references,
                        bool log = false) {
  ASSERT(hyperedges.size() == references.size());
  for (size_t i = 0; i < hyperedges.size(); ++i) {
    const HyperedgeID he = hyperedges[i];
    const std::set<HypernodeID>& reference = references[i];
    size_t count = 0;
    for (const HypernodeID& pin : hypergraph.pins(he)) {
      if (log) LOG << V(he) << V(pin);
      ASSERT_TRUE(reference.find(pin) != reference.end()) << V(he) << V(pin);
      count++;
    }
    ASSERT_EQ(count, reference.size());
  }
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ContainsCorrectNumberofNodesEdgesAndPins) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  ASSERT_EQ(7, hypergraph.initialNumNodes());
  ASSERT_EQ(4, hypergraph.initialNumEdges());
  ASSERT_EQ(12, hypergraph.initialNumPins());
  ASSERT_EQ(7, hypergraph.totalWeight());

  ASSERT_EQ(3, hypergraph.initialNumNodes(0));
  ASSERT_EQ(2, hypergraph.initialNumEdges(0));
  ASSERT_EQ(6, hypergraph.initialNumPins(0));

  ASSERT_EQ(4, hypergraph.initialNumNodes(1));
  ASSERT_EQ(2, hypergraph.initialNumEdges(1));
  ASSERT_EQ(6, hypergraph.initialNumPins(1));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, HasCorrectWeightAfterUpdatingNodeWeights) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  ASSERT_EQ(7, hypergraph.totalWeight());

  hypergraph.setNodeWeight(0, 3);
  hypergraph.setNodeWeight(1, 5);
  hypergraph.setNodeWeight(281474976710656, 10);
  hypergraph.setNodeWeight(281474976710658, 5);
  hypergraph.updateTotalWeight(TBBArena::GLOBAL_TASK_GROUP);

  ASSERT_EQ(26, hypergraph.totalWeight());
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalNodeIterators1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HypernodeID>({ 0, 1, 2 }, [&] {
        return hypergraph.nodes(0);
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalNodeIterators2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HypernodeID>({ 281474976710656, 281474976710657, 281474976710658, 281474976710659 }, [&] {
        return hypergraph.nodes(1);
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalNodeIteratorsWithDisabledHypernodes1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHypernode(1);

  verifyIterator<HypernodeID>({ 0, 2 }, [&] {
        return hypergraph.nodes(0);
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalNodeIteratorsWithDisabledHypernodes2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHypernode(281474976710657);
  hypergraph.disableHypernode(281474976710658);

  verifyIterator<HypernodeID>({ 281474976710656, 281474976710659 }, [&] {
        return hypergraph.nodes(1);
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksOriginalNodeIds1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HypernodeID>({ 0, 1, 2 }, [&] {
        return hypergraph.nodes(0);
      }, [&](const HypernodeID& hn) {
        return hypergraph.originalNodeID(hn);
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksOriginalNodeIds2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HypernodeID>({ 3, 4, 5, 6 }, [&] {
        return hypergraph.nodes(1);
      }, [&](const HypernodeID& hn) {
        return hypergraph.originalNodeID(hn);
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalNodeIterator) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HypernodeID>({ 0, 1, 2, 281474976710656, 281474976710657,
                                281474976710658, 281474976710659 }, [&] {
        return hypergraph.nodes();
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalNodeIteratorIfAllNodesAreAssignedToSecondNode) {
  TestHypergraph hypergraph = this->construct_hypergraph(7,
                                                         { { 0, 2 }, { 0, 1, 3, 4 }, { 3, 4, 6 }, { 2, 5, 6 } },
                                                         { 1, 1, 1, 1, 1, 1, 1 },
                                                         { 1, 1, 1, 1 });

  verifyIterator<HypernodeID>({ 281474976710656, 281474976710657,
                                281474976710658, 281474976710659, 281474976710660, 281474976710661,
                                281474976710662 }, [&] {
        return hypergraph.nodes();
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalNodeIteratorWithDisabledHypernodes1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHypernode(2);
  hypergraph.disableHypernode(281474976710657);
  hypergraph.disableHypernode(281474976710659);

  verifyIterator<HypernodeID>({ 0, 1, 281474976710656, 281474976710658 }, [&] {
        return hypergraph.nodes();
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalNodeIteratorWithDisabledHypernodes2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHypernode(0);
  hypergraph.disableHypernode(2);
  hypergraph.disableHypernode(281474976710656);
  hypergraph.disableHypernode(281474976710659);

  verifyIterator<HypernodeID>({ 1, 281474976710657, 281474976710658 }, [&] {
        return hypergraph.nodes();
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalEdgeIterators1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HyperedgeID>({ 0, 1 }, [&] {
        return hypergraph.edges(0);
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalEdgeIterators2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HyperedgeID>({ 281474976710656, 281474976710657 }, [&] {
        return hypergraph.edges(1);
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalEdgeIteratorsWithDisabledHyperedges1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHyperedge(1);

  verifyIterator<HyperedgeID>({ 0 }, [&] {
        return hypergraph.edges(0);
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksLocalEdgeIteratorsWithDisabledHyperedges2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHyperedge(281474976710656);

  verifyIterator<HyperedgeID>({ 281474976710657 }, [&] {
        return hypergraph.edges(1);
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalEdgeIterators) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HyperedgeID>({ 0, 1, 281474976710656, 281474976710657 }, [&] {
        return hypergraph.edges();
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalEdgeIteratorsWithDisabledHyperedges1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHyperedge(1);
  hypergraph.disableHyperedge(281474976710656);

  verifyIterator<HyperedgeID>({ 0, 281474976710657 }, [&] {
        return hypergraph.edges();
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksGlobalEdgeIteratorsWithDisabledHyperedges2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.disableHyperedge(0);
  hypergraph.disableHyperedge(281474976710657);

  verifyIterator<HyperedgeID>({ 1, 281474976710656 }, [&] {
        return hypergraph.edges();
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksIncidentEdgesOfHypernode1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HyperedgeID>({ 0, 1 }, [&] {
        return hypergraph.incidentEdges(GLOBAL_ID(hypergraph, 0));
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksIncidentEdgesOfHypernode2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HyperedgeID>({ 0, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(GLOBAL_ID(hypergraph, 2));
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksIncidentEdgesOfHypernode3) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HyperedgeID>({ 1, 281474976710656 }, [&] {
        return hypergraph.incidentEdges(GLOBAL_ID(hypergraph, 3));
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksIncidentEdgesOfHypernode4) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HyperedgeID>({ 281474976710657 }, [&] {
        return hypergraph.incidentEdges(GLOBAL_ID(hypergraph, 5));
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksIncidentEdgesOfHypernode5) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HyperedgeID>({ 281474976710656, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(GLOBAL_ID(hypergraph, 6));
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksPinsOfHyperedge1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HypernodeID>({ GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 2) }, [&] {
        return hypergraph.pins(0);
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksPinsOfHyperedge2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HypernodeID>({ GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1),
                                GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4) }, [&] {
        return hypergraph.pins(1);
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksPinsOfHyperedge3) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HypernodeID>({ GLOBAL_ID(hypergraph, 3),
                                GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 6) }, [&] {
        return hypergraph.pins(281474976710656);
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksPinsOfHyperedge4) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  verifyIterator<HypernodeID>({ GLOBAL_ID(hypergraph, 2),
                                GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) }, [&] {
        return hypergraph.pins(281474976710657);
      });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, VerifiesInitialNodeWeights) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  ASSERT_EQ(1, hypergraph.nodeWeight(GLOBAL_ID(hypergraph, 0)));
  ASSERT_EQ(1, hypergraph.nodeWeight(GLOBAL_ID(hypergraph, 1)));
  ASSERT_EQ(1, hypergraph.nodeWeight(GLOBAL_ID(hypergraph, 2)));
  ASSERT_EQ(1, hypergraph.nodeWeight(GLOBAL_ID(hypergraph, 3)));
  ASSERT_EQ(1, hypergraph.nodeWeight(GLOBAL_ID(hypergraph, 4)));
  ASSERT_EQ(1, hypergraph.nodeWeight(GLOBAL_ID(hypergraph, 5)));
  ASSERT_EQ(1, hypergraph.nodeWeight(GLOBAL_ID(hypergraph, 6)));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, VerifiesModifiedNodeWeights) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.setNodeWeight(GLOBAL_ID(hypergraph, 0), 2);
  hypergraph.setNodeWeight(GLOBAL_ID(hypergraph, 1), 3);
  hypergraph.setNodeWeight(GLOBAL_ID(hypergraph, 2), 4);
  hypergraph.setNodeWeight(GLOBAL_ID(hypergraph, 3), 5);
  hypergraph.setNodeWeight(GLOBAL_ID(hypergraph, 4), 6);
  hypergraph.setNodeWeight(GLOBAL_ID(hypergraph, 5), 7);
  hypergraph.setNodeWeight(GLOBAL_ID(hypergraph, 6), 8);

  ASSERT_EQ(2, hypergraph.nodeWeight(GLOBAL_ID(hypergraph, 0)));
  ASSERT_EQ(3, hypergraph.nodeWeight(GLOBAL_ID(hypergraph, 1)));
  ASSERT_EQ(4, hypergraph.nodeWeight(GLOBAL_ID(hypergraph, 2)));
  ASSERT_EQ(5, hypergraph.nodeWeight(GLOBAL_ID(hypergraph, 3)));
  ASSERT_EQ(6, hypergraph.nodeWeight(GLOBAL_ID(hypergraph, 4)));
  ASSERT_EQ(7, hypergraph.nodeWeight(GLOBAL_ID(hypergraph, 5)));
  ASSERT_EQ(8, hypergraph.nodeWeight(GLOBAL_ID(hypergraph, 6)));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, VerifiesInitialEdgeWeights) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  ASSERT_EQ(1, hypergraph.edgeWeight(0));
  ASSERT_EQ(1, hypergraph.edgeWeight(1));
  ASSERT_EQ(1, hypergraph.edgeWeight(281474976710656));
  ASSERT_EQ(1, hypergraph.edgeWeight(281474976710657));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, VerifiesModifiedEdgeWeights) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.setEdgeWeight(0, 2);
  hypergraph.setEdgeWeight(1, 3);
  hypergraph.setEdgeWeight(281474976710656, 4);
  hypergraph.setEdgeWeight(281474976710657, 5);

  ASSERT_EQ(2, hypergraph.edgeWeight(0));
  ASSERT_EQ(3, hypergraph.edgeWeight(1));
  ASSERT_EQ(4, hypergraph.edgeWeight(281474976710656));
  ASSERT_EQ(5, hypergraph.edgeWeight(281474976710657));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ContractsTwoHypernodes1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  HypernodeID u = id[0];
  HypernodeID v = id[2];
  hypergraph.contract(u, v);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(2, hypergraph.nodeWeight(u));

  verifyIterator<HyperedgeID>({ 0, 1, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(u);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[0], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ContractsTwoHypernodes2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  HypernodeID u = id[3];
  HypernodeID v = id[4];
  hypergraph.contract(u, v);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(2, hypergraph.nodeWeight(u));

  verifyIterator<HyperedgeID>({ 1, 281474976710656 }, [&] {
        return hypergraph.incidentEdges(u);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[1], id[3] }, { id[3], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ContractsTwoHypernodes3) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  HypernodeID u = id[6];
  HypernodeID v = id[3];
  hypergraph.contract(u, v);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(2, hypergraph.nodeWeight(u));

  verifyIterator<HyperedgeID>({ 1, 281474976710656, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(u);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[1], id[4], id[6] }, { id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ContractsTwoHypernodes4) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  HypernodeID u = id[5];
  HypernodeID v = id[0];
  hypergraph.contract(u, v);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(2, hypergraph.nodeWeight(u));

  verifyIterator<HyperedgeID>({ 0, 1, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(u);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[2], id[5] }, { id[1], id[3], id[4], id[5] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ContractsTwoHypernodes5) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  HypernodeID u = id[4];
  HypernodeID v = id[1];
  hypergraph.contract(u, v);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(2, hypergraph.nodeWeight(u));

  verifyIterator<HyperedgeID>({ 1, 281474976710656 }, [&] {
        return hypergraph.incidentEdges(u);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[3], id[4] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ContractsTwoHypernodes6) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  HypernodeID u = id[0];
  HypernodeID v = id[6];
  hypergraph.contract(u, v);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(2, hypergraph.nodeWeight(u));

  verifyIterator<HyperedgeID>({ 0, 1, 281474976710656, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(u);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[0], id[3], id[4] }, { id[0], id[2], id[5] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, HasEqualHashIfTwoHyperedgesBecomeParallel1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };

  ASSERT_NE(hypergraph.edgeHash(1), hypergraph.edgeHash(281474976710656));

  hypergraph.contract(id[0], id[6]);
  hypergraph.contract(id[0], id[1]);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(id[0]));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(id[1]));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(id[6]));
  ASSERT_EQ(3, hypergraph.nodeWeight(id[0]));

  verifyIterator<HyperedgeID>({ 0, 1, 281474976710656, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(id[0]);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[3], id[4] }, { id[0], id[3], id[4] }, { id[0], id[2], id[5] } });

  ASSERT_EQ(hypergraph.edgeHash(1), hypergraph.edgeHash(281474976710656));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, HasEqualHashIfTwoHyperedgesBecomeParallel2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };

  ASSERT_NE(hypergraph.edgeHash(0), hypergraph.edgeHash(281474976710657));

  hypergraph.contract(id[0], id[6]);
  hypergraph.contract(id[2], id[5]);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(id[0]));
  ASSERT_TRUE(hypergraph.nodeIsEnabled(id[2]));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(id[5]));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(id[6]));
  ASSERT_EQ(2, hypergraph.nodeWeight(id[0]));
  ASSERT_EQ(2, hypergraph.nodeWeight(id[2]));

  verifyIterator<HyperedgeID>({ 0, 1, 281474976710656, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(id[0]);
      });

  verifyIterator<HyperedgeID>({ 0, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(id[2]);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[0], id[3], id[4] }, { id[0], id[2] } });

  ASSERT_EQ(hypergraph.edgeHash(0), hypergraph.edgeHash(281474976710657));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, HasEqualHashIfTwoHyperedgesBecomeParallel3) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };

  ASSERT_NE(hypergraph.edgeHash(0), hypergraph.edgeHash(281474976710657));
  ASSERT_NE(hypergraph.edgeHash(1), hypergraph.edgeHash(281474976710656));

  hypergraph.contract(id[0], id[6]);
  hypergraph.contract(id[0], id[1]);
  hypergraph.contract(id[2], id[5]);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(id[0]));
  ASSERT_TRUE(hypergraph.nodeIsEnabled(id[2]));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(id[1]));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(id[5]));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(id[6]));
  ASSERT_EQ(3, hypergraph.nodeWeight(id[0]));
  ASSERT_EQ(2, hypergraph.nodeWeight(id[2]));

  verifyIterator<HyperedgeID>({ 0, 1, 281474976710656, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(id[0]);
      });

  verifyIterator<HyperedgeID>({ 0, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(id[2]);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[3], id[4] }, { id[0], id[3], id[4] }, { id[0], id[2] } });

  ASSERT_EQ(hypergraph.edgeHash(0), hypergraph.edgeHash(281474976710657));
  ASSERT_EQ(hypergraph.edgeHash(1), hypergraph.edgeHash(281474976710656));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, HasEqualHashIfTwoHyperedgesBecomeParallel4) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };

  ASSERT_NE(hypergraph.edgeHash(1), hypergraph.edgeHash(281474976710656));

  hypergraph.contract(id[4], id[6]);
  hypergraph.contract(id[0], id[1]);
  hypergraph.contract(id[0], id[3]);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(id[0]));
  ASSERT_TRUE(hypergraph.nodeIsEnabled(id[4]));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(id[1]));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(id[3]));
  ASSERT_FALSE(hypergraph.nodeIsEnabled(id[6]));
  ASSERT_EQ(3, hypergraph.nodeWeight(id[0]));
  ASSERT_EQ(2, hypergraph.nodeWeight(id[4]));

  verifyIterator<HyperedgeID>({ 0, 1, 281474976710656 }, [&] {
        return hypergraph.incidentEdges(id[0]);
      });

  verifyIterator<HyperedgeID>({ 1, 281474976710656, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(id[4]);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[4] }, { id[0], id[4] }, { id[2], id[4], id[5] } });

  ASSERT_EQ(hypergraph.edgeHash(1), hypergraph.edgeHash(281474976710656));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, UncontractsTwoHypernodes1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());

  HypernodeID u = id[0];
  HypernodeID v = id[2];
  auto memento = hypergraph.contract(u, v);
  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(memento, parallel_he_representative);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_TRUE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(1, hypergraph.nodeWeight(u));
  ASSERT_EQ(1, hypergraph.nodeWeight(v));

  verifyIterator<HyperedgeID>({ 0, 1 }, [&] {
        return hypergraph.incidentEdges(u);
      });

  verifyIterator<HyperedgeID>({ 0, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(v);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, UncontractsTwoHypernodes2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());

  HypernodeID u = id[3];
  HypernodeID v = id[4];
  auto memento = hypergraph.contract(u, v);
  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(memento, parallel_he_representative);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_TRUE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(1, hypergraph.nodeWeight(u));
  ASSERT_EQ(1, hypergraph.nodeWeight(v));

  verifyIterator<HyperedgeID>({ 1, 281474976710656 }, [&] {
        return hypergraph.incidentEdges(u);
      });

  verifyIterator<HyperedgeID>({ 1, 281474976710656 }, [&] {
        return hypergraph.incidentEdges(v);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, UncontractsTwoHypernodes3) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());

  HypernodeID u = id[6];
  HypernodeID v = id[3];
  auto memento = hypergraph.contract(u, v);
  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(memento, parallel_he_representative);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_TRUE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(1, hypergraph.nodeWeight(u));
  ASSERT_EQ(1, hypergraph.nodeWeight(v));

  verifyIterator<HyperedgeID>({ 281474976710656, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(u);
      });

  verifyIterator<HyperedgeID>({ 1, 281474976710656 }, [&] {
        return hypergraph.incidentEdges(v);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, UncontractsTwoHypernodes4) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());

  HypernodeID u = id[5];
  HypernodeID v = id[0];
  auto memento = hypergraph.contract(u, v);
  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(memento, parallel_he_representative);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_TRUE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(1, hypergraph.nodeWeight(u));
  ASSERT_EQ(1, hypergraph.nodeWeight(v));

  verifyIterator<HyperedgeID>({ 281474976710657 }, [&] {
        return hypergraph.incidentEdges(u);
      });

  verifyIterator<HyperedgeID>({ 0, 1 }, [&] {
        return hypergraph.incidentEdges(v);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, UncontractsTwoHypernodes5) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());

  HypernodeID u = id[4];
  HypernodeID v = id[1];
  auto memento = hypergraph.contract(u, v);
  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(memento, parallel_he_representative);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_TRUE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(1, hypergraph.nodeWeight(u));
  ASSERT_EQ(1, hypergraph.nodeWeight(v));

  verifyIterator<HyperedgeID>({ 1, 281474976710656 }, [&] {
        return hypergraph.incidentEdges(u);
      });

  verifyIterator<HyperedgeID>({ 1 }, [&] {
        return hypergraph.incidentEdges(v);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, UncontractsTwoHypernodes6) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());

  HypernodeID u = id[0];
  HypernodeID v = id[6];
  auto memento = hypergraph.contract(u, v);
  assignPartitionIDs(hypergraph);
  hypergraph.uncontract(memento, parallel_he_representative);

  ASSERT_TRUE(hypergraph.nodeIsEnabled(u));
  ASSERT_TRUE(hypergraph.nodeIsEnabled(v));
  ASSERT_EQ(1, hypergraph.nodeWeight(u));
  ASSERT_EQ(1, hypergraph.nodeWeight(v));

  verifyIterator<HyperedgeID>({ 0, 1 }, [&] {
        return hypergraph.incidentEdges(u);
      });

  verifyIterator<HyperedgeID>({ 281474976710656, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(v);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, FullyContractsHypergraphAndThenUncontract) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());

  auto memento_1 = hypergraph.contract(id[0], id[2]);
  auto memento_2 = hypergraph.contract(id[5], id[6]);
  auto memento_3 = hypergraph.contract(id[3], id[4]);
  auto memento_4 = hypergraph.contract(id[0], id[1]);
  auto memento_5 = hypergraph.contract(id[0], id[3]);
  auto memento_6 = hypergraph.contract(id[0], id[5]);

  // Initial State
  {
    ASSERT_TRUE(hypergraph.nodeIsEnabled(id[0]));
    ASSERT_EQ(7, hypergraph.nodeWeight(id[0]));

    verifyIterator<HyperedgeID>({ 0, 1, 281474976710656, 281474976710657 }, [&] {
          return hypergraph.incidentEdges(id[0]);
        });

    verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                       { { id[0] }, { id[0] }, { id[0] }, { id[0] } });
  }

  assignPartitionIDs(hypergraph);

  // First Uncontraction
  {
    hypergraph.uncontract(memento_6, parallel_he_representative);

    ASSERT_TRUE(hypergraph.nodeIsEnabled(id[0]));
    ASSERT_TRUE(hypergraph.nodeIsEnabled(id[5]));
    ASSERT_EQ(5, hypergraph.nodeWeight(id[0]));
    ASSERT_EQ(2, hypergraph.nodeWeight(id[5]));

    verifyIterator<HyperedgeID>({ 0, 1, 281474976710656, 281474976710657 }, [&] {
          return hypergraph.incidentEdges(id[0]);
        });

    verifyIterator<HyperedgeID>({ 281474976710656, 281474976710657 }, [&] {
          return hypergraph.incidentEdges(id[5]);
        });

    verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                       { { id[0] }, { id[0] }, { id[0], id[5] }, { id[0], id[5] } });
  }

  // Second Uncontraction
  {
    hypergraph.uncontract(memento_5, parallel_he_representative);

    ASSERT_TRUE(hypergraph.nodeIsEnabled(id[0]));
    ASSERT_TRUE(hypergraph.nodeIsEnabled(id[3]));
    ASSERT_EQ(3, hypergraph.nodeWeight(id[0]));
    ASSERT_EQ(2, hypergraph.nodeWeight(id[3]));

    verifyIterator<HyperedgeID>({ 0, 1, 281474976710657 }, [&] {
          return hypergraph.incidentEdges(id[0]);
        });

    verifyIterator<HyperedgeID>({ 1, 281474976710656 }, [&] {
          return hypergraph.incidentEdges(id[3]);
        });

    verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                       { { id[0] }, { id[0], id[3] }, { id[3], id[5] }, { id[0], id[5] } });
  }

  // Third Uncontraction
  {
    hypergraph.uncontract(memento_4, parallel_he_representative);

    ASSERT_TRUE(hypergraph.nodeIsEnabled(id[0]));
    ASSERT_TRUE(hypergraph.nodeIsEnabled(id[1]));
    ASSERT_EQ(2, hypergraph.nodeWeight(id[0]));
    ASSERT_EQ(1, hypergraph.nodeWeight(id[1]));

    verifyIterator<HyperedgeID>({ 0, 1, 281474976710657 }, [&] {
          return hypergraph.incidentEdges(id[0]);
        });

    verifyIterator<HyperedgeID>({ 1 }, [&] {
          return hypergraph.incidentEdges(id[1]);
        });

    verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                       { { id[0] }, { id[0], id[1], id[3] }, { id[3], id[5] }, { id[0], id[5] } });
  }

  // Fourth Uncontraction
  {
    hypergraph.uncontract(memento_3, parallel_he_representative);

    ASSERT_TRUE(hypergraph.nodeIsEnabled(id[3]));
    ASSERT_TRUE(hypergraph.nodeIsEnabled(id[4]));
    ASSERT_EQ(1, hypergraph.nodeWeight(id[3]));
    ASSERT_EQ(1, hypergraph.nodeWeight(id[4]));

    verifyIterator<HyperedgeID>({ 1, 281474976710656 }, [&] {
          return hypergraph.incidentEdges(id[3]);
        });

    verifyIterator<HyperedgeID>({ 1, 281474976710656 }, [&] {
          return hypergraph.incidentEdges(id[4]);
        });

    verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                       { { id[0] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[5] }, { id[0], id[5] } });
  }

  // Fifth Uncontraction
  {
    hypergraph.uncontract(memento_2, parallel_he_representative);

    ASSERT_TRUE(hypergraph.nodeIsEnabled(id[5]));
    ASSERT_TRUE(hypergraph.nodeIsEnabled(id[6]));
    ASSERT_EQ(1, hypergraph.nodeWeight(id[5]));
    ASSERT_EQ(1, hypergraph.nodeWeight(id[6]));

    verifyIterator<HyperedgeID>({ 281474976710657 }, [&] {
          return hypergraph.incidentEdges(id[5]);
        });

    verifyIterator<HyperedgeID>({ 281474976710656, 281474976710657 }, [&] {
          return hypergraph.incidentEdges(id[6]);
        });

    verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                       { { id[0] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[0], id[5], id[6] } });
  }

  // Sixth Uncontraction
  {
    hypergraph.uncontract(memento_1, parallel_he_representative);

    ASSERT_TRUE(hypergraph.nodeIsEnabled(id[0]));
    ASSERT_TRUE(hypergraph.nodeIsEnabled(id[2]));
    ASSERT_EQ(1, hypergraph.nodeWeight(id[0]));
    ASSERT_EQ(1, hypergraph.nodeWeight(id[2]));

    verifyIterator<HyperedgeID>({ 0, 1 }, [&] {
          return hypergraph.incidentEdges(id[0]);
        });

    verifyIterator<HyperedgeID>({ 0, 281474976710657 }, [&] {
          return hypergraph.incidentEdges(id[2]);
        });

    verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                       { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
  }
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, RemovesAnEdgeFromHypergraph1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };

  hypergraph.removeEdge(0);

  ASSERT_FALSE(hypergraph.edgeIsEnabled(0));

  verifyIterator<HyperedgeID>({ 1 }, [&] {
        return hypergraph.incidentEdges(id[0]);
      });

  verifyIterator<HyperedgeID>({ 281474976710657 }, [&] {
        return hypergraph.incidentEdges(id[2]);
      });

  verifyPinIterators(hypergraph, { 1, 281474976710656, 281474976710657 },
                     { { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, RemovesAnEdgeFromHypergraph2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };

  hypergraph.removeEdge(1);

  ASSERT_FALSE(hypergraph.edgeIsEnabled(1));

  verifyIterator<HyperedgeID>({ 0 }, [&] {
        return hypergraph.incidentEdges(id[0]);
      });

  verifyIterator<HyperedgeID>({ }, [&] {
        return hypergraph.incidentEdges(id[1]);
      });

  verifyIterator<HyperedgeID>({ 281474976710656 }, [&] {
        return hypergraph.incidentEdges(id[3]);
      });

  verifyIterator<HyperedgeID>({ 281474976710656 }, [&] {
        return hypergraph.incidentEdges(id[4]);
      });

  verifyPinIterators(hypergraph, { 0, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, RemovesAnEdgeFromHypergraph3) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };

  hypergraph.removeEdge(281474976710656);

  ASSERT_FALSE(hypergraph.edgeIsEnabled(281474976710656));

  verifyIterator<HyperedgeID>({ 1 }, [&] {
        return hypergraph.incidentEdges(id[3]);
      });

  verifyIterator<HyperedgeID>({ 1 }, [&] {
        return hypergraph.incidentEdges(id[4]);
      });

  verifyIterator<HyperedgeID>({ 281474976710657 }, [&] {
        return hypergraph.incidentEdges(id[6]);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, RemovesAnEdgeFromHypergraph4) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };

  hypergraph.removeEdge(281474976710657);

  ASSERT_FALSE(hypergraph.edgeIsEnabled(281474976710657));

  verifyIterator<HyperedgeID>({ 0 }, [&] {
        return hypergraph.incidentEdges(id[2]);
      });

  verifyIterator<HyperedgeID>({ }, [&] {
        return hypergraph.incidentEdges(id[5]);
      });

  verifyIterator<HyperedgeID>({ 281474976710656 }, [&] {
        return hypergraph.incidentEdges(id[6]);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656 },
                     { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, RemovesTwoEdgesFromHypergraph1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };

  hypergraph.removeEdge(0);
  hypergraph.removeEdge(1);

  ASSERT_FALSE(hypergraph.edgeIsEnabled(0));
  ASSERT_FALSE(hypergraph.edgeIsEnabled(1));

  verifyIterator<HyperedgeID>({ }, [&] {
        return hypergraph.incidentEdges(id[0]);
      });

  verifyIterator<HyperedgeID>({ }, [&] {
        return hypergraph.incidentEdges(id[1]);
      });

  verifyIterator<HyperedgeID>({ 281474976710657 }, [&] {
        return hypergraph.incidentEdges(id[2]);
      });

  verifyIterator<HyperedgeID>({ 281474976710656 }, [&] {
        return hypergraph.incidentEdges(id[3]);
      });

  verifyIterator<HyperedgeID>({ 281474976710656 }, [&] {
        return hypergraph.incidentEdges(id[4]);
      });

  verifyIterator<HyperedgeID>({ 281474976710657 }, [&] {
        return hypergraph.incidentEdges(id[5]);
      });

  verifyIterator<HyperedgeID>({ 281474976710656, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(id[6]);
      });

  verifyPinIterators(hypergraph, { 281474976710656, 281474976710657 },
                     { { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, RemovesTwoEdgesFromHypergraph2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };

  hypergraph.removeEdge(281474976710656);
  hypergraph.removeEdge(281474976710657);

  ASSERT_FALSE(hypergraph.edgeIsEnabled(281474976710656));
  ASSERT_FALSE(hypergraph.edgeIsEnabled(281474976710657));

  verifyIterator<HyperedgeID>({ 0, 1 }, [&] {
        return hypergraph.incidentEdges(id[0]);
      });

  verifyIterator<HyperedgeID>({ 1 }, [&] {
        return hypergraph.incidentEdges(id[1]);
      });

  verifyIterator<HyperedgeID>({ 0 }, [&] {
        return hypergraph.incidentEdges(id[2]);
      });

  verifyIterator<HyperedgeID>({ 1 }, [&] {
        return hypergraph.incidentEdges(id[3]);
      });

  verifyIterator<HyperedgeID>({ 1 }, [&] {
        return hypergraph.incidentEdges(id[4]);
      });

  verifyIterator<HyperedgeID>({ }, [&] {
        return hypergraph.incidentEdges(id[5]);
      });

  verifyIterator<HyperedgeID>({ }, [&] {
        return hypergraph.incidentEdges(id[6]);
      });

  verifyPinIterators(hypergraph, { 0, 1 },
                     { { id[0], id[2] }, { id[0], id[1], id[3], id[4] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, RestoreAnEdgeOfHypergraph1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };

  hypergraph.removeEdge(0);
  assignPartitionIDs(hypergraph);
  hypergraph.restoreEdge(0, 2);

  ASSERT_TRUE(hypergraph.edgeIsEnabled(0));

  verifyIterator<HyperedgeID>({ 0, 1 }, [&] {
        return hypergraph.incidentEdges(id[0]);
      });

  verifyIterator<HyperedgeID>({ 0, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(id[2]);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, RestoreAnEdgeOfHypergraph2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };

  hypergraph.removeEdge(1);
  assignPartitionIDs(hypergraph);
  hypergraph.restoreEdge(1, 4);

  ASSERT_TRUE(hypergraph.edgeIsEnabled(1));

  verifyIterator<HyperedgeID>({ 0, 1 }, [&] {
        return hypergraph.incidentEdges(id[0]);
      });

  verifyIterator<HyperedgeID>({ 1 }, [&] {
        return hypergraph.incidentEdges(id[1]);
      });

  verifyIterator<HyperedgeID>({ 1, 281474976710656 }, [&] {
        return hypergraph.incidentEdges(id[3]);
      });

  verifyIterator<HyperedgeID>({ 1, 281474976710656 }, [&] {
        return hypergraph.incidentEdges(id[4]);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, RestoreAnEdgeOfHypergraph3) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };

  hypergraph.removeEdge(281474976710656);
  assignPartitionIDs(hypergraph);
  hypergraph.restoreEdge(281474976710656, 3);

  ASSERT_TRUE(hypergraph.edgeIsEnabled(281474976710656));

  verifyIterator<HyperedgeID>({ 1, 281474976710656 }, [&] {
        return hypergraph.incidentEdges(id[3]);
      });

  verifyIterator<HyperedgeID>({ 1, 281474976710656 }, [&] {
        return hypergraph.incidentEdges(id[4]);
      });

  verifyIterator<HyperedgeID>({ 281474976710656, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(id[6]);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, RestoreAnEdgeOfHypergraph4) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };

  hypergraph.removeEdge(281474976710657);
  assignPartitionIDs(hypergraph);
  hypergraph.restoreEdge(281474976710657, 3);

  ASSERT_TRUE(hypergraph.edgeIsEnabled(281474976710657));

  verifyIterator<HyperedgeID>({ 0, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(id[2]);
      });

  verifyIterator<HyperedgeID>({ 281474976710657 }, [&] {
        return hypergraph.incidentEdges(id[5]);
      });

  verifyIterator<HyperedgeID>({ 281474976710656, 281474976710657 }, [&] {
        return hypergraph.incidentEdges(id[6]);
      });

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ContractionsIntermixedWithEdgeRemovals) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());

  auto memento_1 = hypergraph.contract(id[0], id[2]);
  hypergraph.removeEdge(0);
  auto memento_2 = hypergraph.contract(id[4], id[1]);
  auto memento_3 = hypergraph.contract(id[3], id[6]);
  hypergraph.removeEdge(281474976710656);

  verifyPinIterators(hypergraph, { 1, 281474976710657 },
                     { { id[0], id[3], id[4] }, { id[0], id[5], id[3] } });

  assignPartitionIDs(hypergraph);

  hypergraph.restoreEdge(281474976710656, 2);
  hypergraph.uncontract(memento_3, parallel_he_representative);
  hypergraph.uncontract(memento_2, parallel_he_representative);
  hypergraph.restoreSinglePinHyperedge(0);
  hypergraph.uncontract(memento_1, parallel_he_representative);

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { id[0], id[2] }, { id[0], id[1], id[3], id[4] }, { id[3], id[4], id[6] }, { id[2], id[5], id[6] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, StreamsCommunityIDsInParallelIntoHypergraph) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<PartitionID> communities = { 0, 0, 0, 1, 1, 2, 2 };

  tbb::parallel_for(0UL, hypergraph.initialNumNodes(), [&](const HypernodeID& hn) {
        hypergraph.setCommunityID(hypergraph.globalNodeID(hn), communities[hn]);
      });
  hypergraph.initializeCommunities();

  ASSERT_EQ(3, hypergraph.numCommunities());
  ASSERT_EQ(0, hypergraph.communityID(GLOBAL_ID(hypergraph, 0)));
  ASSERT_EQ(0, hypergraph.communityID(GLOBAL_ID(hypergraph, 1)));
  ASSERT_EQ(0, hypergraph.communityID(GLOBAL_ID(hypergraph, 2)));
  ASSERT_EQ(1, hypergraph.communityID(GLOBAL_ID(hypergraph, 3)));
  ASSERT_EQ(1, hypergraph.communityID(GLOBAL_ID(hypergraph, 4)));
  ASSERT_EQ(2, hypergraph.communityID(GLOBAL_ID(hypergraph, 5)));
  ASSERT_EQ(2, hypergraph.communityID(GLOBAL_ID(hypergraph, 6)));

  ASSERT_EQ(0, hypergraph.communityNodeId(GLOBAL_ID(hypergraph, 0)));
  ASSERT_EQ(1, hypergraph.communityNodeId(GLOBAL_ID(hypergraph, 1)));
  ASSERT_EQ(2, hypergraph.communityNodeId(GLOBAL_ID(hypergraph, 2)));
  ASSERT_EQ(0, hypergraph.communityNodeId(GLOBAL_ID(hypergraph, 3)));
  ASSERT_EQ(1, hypergraph.communityNodeId(GLOBAL_ID(hypergraph, 4)));
  ASSERT_EQ(0, hypergraph.communityNodeId(GLOBAL_ID(hypergraph, 5)));
  ASSERT_EQ(1, hypergraph.communityNodeId(GLOBAL_ID(hypergraph, 6)));

  ASSERT_EQ(3, hypergraph.numCommunityHypernodes(0));
  ASSERT_EQ(2, hypergraph.numCommunityHypernodes(1));
  ASSERT_EQ(2, hypergraph.numCommunityHypernodes(2));

  ASSERT_EQ(5, hypergraph.numCommunityPins(0));
  ASSERT_EQ(4, hypergraph.numCommunityPins(1));
  ASSERT_EQ(3, hypergraph.numCommunityPins(2));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ResetsPinsToOriginalIds) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  hypergraph.resetPinsToOriginalNodeIds(TBBArena::GLOBAL_TASK_GROUP);

  verifyPinIterators(hypergraph, { 0, 1, 281474976710656, 281474976710657 },
                     { { 0, 2 }, { 0, 1, 3, 4 }, { 3, 4, 6 }, { 2, 5, 6 } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, SetPartIdsOfVertex1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  ASSERT_TRUE(hypergraph.setNodePart(0, 0));
  ASSERT_EQ(0, hypergraph.partID(0));
  ASSERT_EQ(1, hypergraph.localPartWeight(0));
  ASSERT_EQ(1, hypergraph.localPartSize(0));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, SetPartIdsOfVertex2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  ASSERT_TRUE(hypergraph.setNodePart(281474976710656, 1));
  ASSERT_EQ(1, hypergraph.partID(281474976710656));
  ASSERT_EQ(1, hypergraph.localPartWeight(1));
  ASSERT_EQ(1, hypergraph.localPartSize(1));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, SetPartIdsOfAllVertices) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  assignPartitionIDs(hypergraph);

  ASSERT_EQ(0, hypergraph.partID(0));
  ASSERT_EQ(0, hypergraph.partID(1));
  ASSERT_EQ(0, hypergraph.partID(2));
  ASSERT_EQ(1, hypergraph.partID(281474976710656));
  ASSERT_EQ(1, hypergraph.partID(281474976710657));
  ASSERT_EQ(1, hypergraph.partID(281474976710658));
  ASSERT_EQ(1, hypergraph.partID(281474976710659));

  ASSERT_EQ(3, hypergraph.partWeight(0));
  ASSERT_EQ(3, hypergraph.partSize(0));
  ASSERT_EQ(4, hypergraph.partWeight(1));
  ASSERT_EQ(4, hypergraph.partSize(1));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, UpdatePartIdsOfOneVertex) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  assignPartitionIDs(hypergraph);
  ASSERT_TRUE(hypergraph.changeNodePart(2, 0, 1));

  ASSERT_EQ(1, hypergraph.partID(2));
  ASSERT_EQ(2, hypergraph.localPartWeight(0));
  ASSERT_EQ(2, hypergraph.localPartSize(0));
  ASSERT_EQ(5, hypergraph.localPartWeight(1));
  ASSERT_EQ(5, hypergraph.localPartSize(1));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, UpdatePartIdsOfTwoVertices) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  assignPartitionIDs(hypergraph);
  ASSERT_TRUE(hypergraph.changeNodePart(2, 0, 1));
  ASSERT_TRUE(hypergraph.changeNodePart(281474976710659, 1, 0));

  ASSERT_EQ(0, hypergraph.partID(281474976710659));
  ASSERT_EQ(3, hypergraph.localPartWeight(0));
  ASSERT_EQ(3, hypergraph.localPartSize(0));
  ASSERT_EQ(4, hypergraph.localPartWeight(1));
  ASSERT_EQ(4, hypergraph.localPartSize(1));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, UpdatePartIdsOfThreeVertices) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);

  assignPartitionIDs(hypergraph);
  ASSERT_TRUE(hypergraph.changeNodePart(2, 0, 1));
  ASSERT_TRUE(hypergraph.changeNodePart(281474976710659, 1, 0));
  ASSERT_TRUE(hypergraph.changeNodePart(281474976710658, 1, 0));

  ASSERT_EQ(0, hypergraph.partID(281474976710658));
  ASSERT_EQ(4, hypergraph.localPartWeight(0));
  ASSERT_EQ(4, hypergraph.localPartSize(0));
  ASSERT_EQ(3, hypergraph.localPartWeight(1));
  ASSERT_EQ(3, hypergraph.localPartSize(1));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksPartIdAssignmentAfterUncontraction1) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());

  HypernodeID u = id[0];
  HypernodeID v = id[2];
  auto memento = hypergraph.contract(u, v);

  assignPartitionIDs(hypergraph);
  ASSERT_EQ(3, hypergraph.partWeight(0));
  ASSERT_EQ(2, hypergraph.partSize(0));
  ASSERT_EQ(4, hypergraph.partWeight(1));
  ASSERT_EQ(4, hypergraph.partSize(1));

  hypergraph.uncontract(memento, parallel_he_representative);
  hypergraph.updateGlobalPartInfos();
  ASSERT_EQ(0, hypergraph.partID(v));
  ASSERT_EQ(3, hypergraph.partWeight(0));
  ASSERT_EQ(3, hypergraph.partSize(0));
  ASSERT_EQ(4, hypergraph.partWeight(1));
  ASSERT_EQ(4, hypergraph.partSize(1));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksPartIdAssignmentAfterUncontraction2) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());

  HypernodeID u = id[3];
  HypernodeID v = id[4];
  auto memento = hypergraph.contract(u, v);

  assignPartitionIDs(hypergraph);
  ASSERT_EQ(3, hypergraph.partWeight(0));
  ASSERT_EQ(3, hypergraph.partSize(0));
  ASSERT_EQ(4, hypergraph.partWeight(1));
  ASSERT_EQ(3, hypergraph.partSize(1));

  hypergraph.uncontract(memento, parallel_he_representative);
  hypergraph.updateGlobalPartInfos();
  ASSERT_EQ(1, hypergraph.partID(v));
  ASSERT_EQ(3, hypergraph.partWeight(0));
  ASSERT_EQ(3, hypergraph.partSize(0));
  ASSERT_EQ(4, hypergraph.partWeight(1));
  ASSERT_EQ(4, hypergraph.partSize(1));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksPartIdAssignmentAfterUncontraction3) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());

  HypernodeID u = id[0];
  HypernodeID v = id[6];
  auto memento = hypergraph.contract(u, v);

  assignPartitionIDs(hypergraph);
  ASSERT_EQ(4, hypergraph.partWeight(0));
  ASSERT_EQ(3, hypergraph.partSize(0));
  ASSERT_EQ(3, hypergraph.partWeight(1));
  ASSERT_EQ(3, hypergraph.partSize(1));

  hypergraph.uncontract(memento, parallel_he_representative);
  hypergraph.updateGlobalPartInfos();
  ASSERT_EQ(0, hypergraph.partID(v));
  ASSERT_EQ(4, hypergraph.partWeight(0));
  ASSERT_EQ(4, hypergraph.partSize(0));
  ASSERT_EQ(3, hypergraph.partWeight(1));
  ASSERT_EQ(3, hypergraph.partSize(1));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ChecksPartIdAssignmentAfterSeveralUncontractions) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2),
                                  GLOBAL_ID(hypergraph, 3), GLOBAL_ID(hypergraph, 4), GLOBAL_ID(hypergraph, 5), GLOBAL_ID(hypergraph, 6) };
  parallel::scalable_vector<HyperedgeID> parallel_he_representative(hypergraph.initialNumEdges());

  auto memento_1 = hypergraph.contract(id[0], id[1]);
  auto memento_2 = hypergraph.contract(id[3], id[4]);
  auto memento_3 = hypergraph.contract(id[5], id[2]);

  assignPartitionIDs(hypergraph);
  ASSERT_EQ(2, hypergraph.partWeight(0));
  ASSERT_EQ(1, hypergraph.partSize(0));
  ASSERT_EQ(5, hypergraph.partWeight(1));
  ASSERT_EQ(3, hypergraph.partSize(1));

  hypergraph.uncontract(memento_3, parallel_he_representative);
  hypergraph.updateGlobalPartInfos();
  ASSERT_EQ(1, hypergraph.partID(id[2]));
  ASSERT_EQ(2, hypergraph.partWeight(0));
  ASSERT_EQ(1, hypergraph.partSize(0));
  ASSERT_EQ(5, hypergraph.partWeight(1));
  ASSERT_EQ(4, hypergraph.partSize(1));

  hypergraph.uncontract(memento_2, parallel_he_representative);
  hypergraph.updateGlobalPartInfos();
  ASSERT_EQ(1, hypergraph.partID(id[4]));
  ASSERT_EQ(2, hypergraph.partWeight(0));
  ASSERT_EQ(1, hypergraph.partSize(0));
  ASSERT_EQ(5, hypergraph.partWeight(1));
  ASSERT_EQ(5, hypergraph.partSize(1));

  hypergraph.uncontract(memento_1, parallel_he_representative);
  hypergraph.updateGlobalPartInfos();
  ASSERT_EQ(1, hypergraph.partID(id[4]));
  ASSERT_EQ(2, hypergraph.partWeight(0));
  ASSERT_EQ(2, hypergraph.partSize(0));
  ASSERT_EQ(5, hypergraph.partWeight(1));
  ASSERT_EQ(5, hypergraph.partSize(1));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, RestoresParallelHyperedgeDuringUncontraction1) {
  TestHypergraph hypergraph = construct_hypergraph(3, { { 0, 1, 2 }, { 1, 2 } },
                                                   { 0, 0, 0 }, { 0, 0 }, { 0, 0, 0 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2) };
  // Hyperedge 1 becomes parallel to hyperedge 0
  auto memento = hypergraph.contract(id[0], id[1]);
  hypergraph.disableHyperedge(1);
  hypergraph.invalidateDisabledHyperedgesFromIncidentNets(TBBArena::GLOBAL_TASK_GROUP);
  assignPartitionIDs(hypergraph);

  parallel::scalable_vector<HyperedgeID> parallel_he_representative(
    hypergraph.initialNumEdges(), std::numeric_limits<HyperedgeID>::max());
  parallel_he_representative[1] = 0;  // Indicating that hyperedge 1 is parallel to 0

  verifyPinIterators(hypergraph, { 0 },
                     { { id[0], id[2] } });

  // Should restore hyperedge 1
  hypergraph.uncontract(memento, parallel_he_representative);

  verifyPinIterators(hypergraph, { 0, 1 },
                     { { id[0], id[1], id[2] }, { id[1], id[2] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, RestoresParallelHyperedgeDuringUncontraction2) {
  TestHypergraph hypergraph = construct_hypergraph(3, { { 0, 1, 2 }, { 1, 2 } },
                                                   { 0, 0, 0 }, { 0, 0 }, { 0, 0, 0 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2) };
  // Hyperedge 0 becomes parallel to hyperedge 1
  auto memento = hypergraph.contract(id[0], id[1]);
  hypergraph.disableHyperedge(0);
  hypergraph.invalidateDisabledHyperedgesFromIncidentNets(TBBArena::GLOBAL_TASK_GROUP);
  assignPartitionIDs(hypergraph);

  parallel::scalable_vector<HyperedgeID> parallel_he_representative(
    hypergraph.initialNumEdges(), std::numeric_limits<HyperedgeID>::max());
  parallel_he_representative[0] = 1;  // Indicating that hyperedge 0 is parallel to 1

  verifyPinIterators(hypergraph, { 1 },
                     { { id[0], id[2] } });

  // Should restore hyperedge 0
  hypergraph.uncontract(memento, parallel_he_representative);

  verifyPinIterators(hypergraph, { 0, 1 },
                     { { id[0], id[1], id[2] }, { id[1], id[2] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, RestoresParallelHyperedgeDuringUncontraction3) {
  TestHypergraph hypergraph = construct_hypergraph(3, { { 0, 2 }, { 1, 2 } },
                                                   { 0, 0, 0 }, { 0, 0 }, { 0, 0, 0 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2) };
  // Hyperedge 1 becomes parallel to hyperedge 0
  auto memento = hypergraph.contract(id[0], id[1]);
  hypergraph.disableHyperedge(1);
  hypergraph.invalidateDisabledHyperedgesFromIncidentNets(TBBArena::GLOBAL_TASK_GROUP);
  assignPartitionIDs(hypergraph);

  parallel::scalable_vector<HyperedgeID> parallel_he_representative(
    hypergraph.initialNumEdges(), std::numeric_limits<HyperedgeID>::max());
  parallel_he_representative[1] = 0;  // Indicating that hyperedge 1 is parallel to 0

  verifyPinIterators(hypergraph, { 0 },
                     { { id[0], id[2] } });

  // Should restore hyperedge 1
  hypergraph.uncontract(memento, parallel_he_representative);

  verifyPinIterators(hypergraph, { 0, 1 },
                     { { id[0], id[2] }, { id[1], id[2] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, RestoresParallelHyperedgeDuringUncontraction4) {
  TestHypergraph hypergraph = construct_hypergraph(3, { { 0, 2 }, { 1, 2 } },
                                                   { 0, 0, 0 }, { 0, 0 }, { 0, 0, 0 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2) };
  // Hyperedge 0 becomes parallel to hyperedge 1
  auto memento = hypergraph.contract(id[0], id[1]);
  hypergraph.disableHyperedge(0);
  hypergraph.invalidateDisabledHyperedgesFromIncidentNets(TBBArena::GLOBAL_TASK_GROUP);
  assignPartitionIDs(hypergraph);

  parallel::scalable_vector<HyperedgeID> parallel_he_representative(
    hypergraph.initialNumEdges(), std::numeric_limits<HyperedgeID>::max());
  parallel_he_representative[0] = 1;  // Indicating that hyperedge 0 is parallel to 1

  verifyPinIterators(hypergraph, { 1 },
                     { { id[0], id[2] } });

  // Should restore hyperedge 0
  hypergraph.uncontract(memento, parallel_he_representative);

  verifyPinIterators(hypergraph, { 0, 1 },
                     { { id[0], id[2] }, { id[1], id[2] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, RestoresParallelHyperedgeDuringUncontraction5) {
  TestHypergraph hypergraph = construct_hypergraph(4, { { 1, 2 }, { 0, 1, 2, 3 }, { 2, 3 } },
                                                   { 0, 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0, 0 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2), GLOBAL_ID(hypergraph, 3) };
  // Hyperedge 2 becomes parallel to hyperedge 0
  auto memento_1 = hypergraph.contract(id[1], id[3]);
  // Hyperede 1 becomes parallel to hyperedge 0
  auto memento_2 = hypergraph.contract(id[0], id[2]);
  hypergraph.disableHyperedge(1);
  hypergraph.disableHyperedge(2);
  hypergraph.invalidateDisabledHyperedgesFromIncidentNets(TBBArena::GLOBAL_TASK_GROUP);
  assignPartitionIDs(hypergraph);

  parallel::scalable_vector<HyperedgeID> parallel_he_representative(
    hypergraph.initialNumEdges(), std::numeric_limits<HyperedgeID>::max());
  parallel_he_representative[1] = 0;  // Indicating that hyperedge 1 is parallel to 0
  parallel_he_representative[2] = 1;  // Indicating that hyperedge 2 is parallel to 1

  verifyPinIterators(hypergraph, { 0 },
                     { { id[0], id[1] } });

  // Should restore hyperedge 1 and 2
  hypergraph.uncontract(memento_2, parallel_he_representative);

  verifyPinIterators(hypergraph, { 0, 1, 2 },
                     { { id[1], id[2] }, { id[0], id[1], id[2] }, { id[1], id[2] } });

  hypergraph.uncontract(memento_1, parallel_he_representative);

  verifyPinIterators(hypergraph, { 0, 1, 2 },
                     { { id[1], id[2] }, { id[0], id[1], id[2], id[3] }, { id[2], id[3] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, CopiesItselfWhereAllVerticesAOnOneNumaNode) {
  TestHypergraph hypergraph = construct_hypergraph(4, { { 1, 2 }, { 0, 1, 2, 3 }, { 2, 3 } },
                                                   { 0, 0, 0, 0 }, { 0, 0, 0 }, { 0, 0, 0, 0 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2), GLOBAL_ID(hypergraph, 3) };
  auto copy = hypergraph.copy(2, TBBArena::GLOBAL_TASK_GROUP);
  TestHypergraph& copy_hg = copy.first;
  auto& mapping = copy.second;
  std::vector<HypernodeID> copy_id = { GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[0])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[1])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[2])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[3])]) };

  verifyPinIterators(copy_hg, { copy_hg.globalEdgeID(0), copy_hg.globalEdgeID(1), copy_hg.globalEdgeID(2) },
                     { { copy_id[1], copy_id[2] }, { copy_id[0], copy_id[1], copy_id[2], copy_id[3] },
                       { copy_id[2], copy_id[3] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, CopiesItselfWhereVerticesAreOnDifferentNumaNodes) {
  TestHypergraph hypergraph = construct_hypergraph(4, { { 1, 2 }, { 0, 1, 2, 3 }, { 2, 3 } },
                                                   { 0, 0, 1, 1 }, { 0, 1, 0 }, { 0, 0, 0, 0 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2), GLOBAL_ID(hypergraph, 3) };
  auto copy = hypergraph.copy(2, TBBArena::GLOBAL_TASK_GROUP);
  TestHypergraph& copy_hg = copy.first;
  auto& mapping = copy.second;
  std::vector<HypernodeID> copy_id = { GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[0])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[1])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[2])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[3])]) };

  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(copy_id[0]));
  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(copy_id[1]));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(copy_id[2]));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(copy_id[3]));

  verifyPinIterators(copy_hg, { copy_hg.globalEdgeID(0), copy_hg.globalEdgeID(1), copy_hg.globalEdgeID(2) },
                     { { copy_id[1], copy_id[2] }, { copy_id[2], copy_id[3] },
                       { copy_id[0], copy_id[1], copy_id[2], copy_id[3] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, CopiesItselfAfterOneContraction) {
  TestHypergraph hypergraph = construct_hypergraph(4, { { 1, 2 }, { 0, 1, 2, 3 }, { 2, 3 } },
                                                   { 0, 0, 1, 1 }, { 0, 1, 0 }, { 0, 0, 0, 0 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2), GLOBAL_ID(hypergraph, 3) };
  hypergraph.contract(id[1], id[2]);

  auto copy = hypergraph.copy(2, TBBArena::GLOBAL_TASK_GROUP);
  TestHypergraph& copy_hg = copy.first;
  auto& mapping = copy.second;
  std::vector<HypernodeID> copy_id = { GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[0])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[1])]),
                                       0,
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[3])]) };

  ASSERT_EQ(1, copy_hg.nodeWeight(copy_id[0]));
  ASSERT_EQ(2, copy_hg.nodeWeight(copy_id[1]));
  ASSERT_EQ(1, copy_hg.nodeWeight(copy_id[3]));

  verifyPinIterators(copy_hg, { copy_hg.globalEdgeID(0), copy_hg.globalEdgeID(1), copy_hg.globalEdgeID(2) },
                     { { copy_id[1] }, { copy_id[1], copy_id[3] }, { copy_id[0], copy_id[1], copy_id[3] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, CopiesItselfAfterTwoContractions) {
  TestHypergraph hypergraph = construct_hypergraph(4, { { 1, 2 }, { 0, 1, 2, 3 }, { 2, 3 } },
                                                   { 0, 0, 1, 1 }, { 0, 1, 0 }, { 0, 0, 0, 0 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2), GLOBAL_ID(hypergraph, 3) };
  hypergraph.contract(id[1], id[2]);
  hypergraph.contract(id[1], id[3]);

  auto copy = hypergraph.copy(2, TBBArena::GLOBAL_TASK_GROUP);
  TestHypergraph& copy_hg = copy.first;
  auto& mapping = copy.second;
  std::vector<HypernodeID> copy_id = { GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[0])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[1])]),
                                       0, 0 };

  ASSERT_EQ(1, copy_hg.nodeWeight(copy_id[0]));
  ASSERT_EQ(3, copy_hg.nodeWeight(copy_id[1]));

  verifyPinIterators(copy_hg, { copy_hg.globalEdgeID(0), copy_hg.globalEdgeID(1), copy_hg.globalEdgeID(2) },
                     { { copy_id[1] }, { copy_id[1] }, { copy_id[0], copy_id[1] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, CopiesItselfAfterTwoContractionsWithDisabledHyperedge) {
  TestHypergraph hypergraph = construct_hypergraph(4, { { 1, 2 }, { 0, 1, 2, 3 }, { 2, 3 } },
                                                   { 0, 0, 1, 1 }, { 0, 1, 0 }, { 0, 0, 0, 0 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2), GLOBAL_ID(hypergraph, 3) };
  hypergraph.contract(id[1], id[2]);
  hypergraph.contract(id[1], id[3]);
  hypergraph.removeEdge(hypergraph.globalEdgeID(2));
  hypergraph.setEdgeWeight(hypergraph.globalEdgeID(0), 2);

  auto copy = hypergraph.copy(2, TBBArena::GLOBAL_TASK_GROUP);
  TestHypergraph& copy_hg = copy.first;
  auto& mapping = copy.second;
  std::vector<HypernodeID> copy_id = { GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[0])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[1])]),
                                       0, 0 };

  ASSERT_EQ(1, copy_hg.nodeWeight(copy_id[0]));
  ASSERT_EQ(3, copy_hg.nodeWeight(copy_id[1]));
  ASSERT_EQ(2, copy_hg.edgeWeight(copy_hg.globalEdgeID(0)));

  verifyPinIterators(copy_hg, { copy_hg.globalEdgeID(0), copy_hg.globalEdgeID(1) },
                     { { copy_id[1] }, { copy_id[0], copy_id[1] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, CopiesItselfWithCommunities) {
  TestHypergraph hypergraph = construct_hypergraph(4, { { 1, 2 }, { 0, 1, 2, 3 }, { 2, 3 } },
                                                   { 0, 0, 1, 1 }, { 0, 1, 0 }, { 0, 1, 2, 2 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2), GLOBAL_ID(hypergraph, 3) };
  auto copy = hypergraph.copy(2, TBBArena::GLOBAL_TASK_GROUP);
  TestHypergraph& copy_hg = copy.first;
  auto& mapping = copy.second;
  std::vector<HypernodeID> copy_id = { GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[0])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[1])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[2])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[3])]) };

  ASSERT_EQ(0, copy_hg.communityID(copy_id[0]));
  ASSERT_EQ(1, copy_hg.communityID(copy_id[1]));
  ASSERT_EQ(2, copy_hg.communityID(copy_id[2]));
  ASSERT_EQ(2, copy_hg.communityID(copy_id[3]));

  verifyPinIterators(copy_hg, { copy_hg.globalEdgeID(0), copy_hg.globalEdgeID(1), copy_hg.globalEdgeID(2) },
                     { { copy_id[1], copy_id[2] }, { copy_id[2], copy_id[3] },
                       { copy_id[0], copy_id[1], copy_id[2], copy_id[3] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, CopiesItselfWithCommunitiesAfterOneContraction) {
  TestHypergraph hypergraph = construct_hypergraph(4, { { 1, 2 }, { 0, 1, 2, 3 }, { 2, 3 } },
                                                   { 0, 0, 1, 1 }, { 0, 1, 0 }, { 0, 1, 2, 2 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2), GLOBAL_ID(hypergraph, 3) };
  hypergraph.contract(id[1], id[2]);

  auto copy = hypergraph.copy(2, TBBArena::GLOBAL_TASK_GROUP);
  TestHypergraph& copy_hg = copy.first;
  auto& mapping = copy.second;
  std::vector<HypernodeID> copy_id = { GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[0])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[1])]),
                                       0,
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[3])]) };

  ASSERT_EQ(0, copy_hg.communityID(copy_id[0]));
  ASSERT_EQ(1, copy_hg.communityID(copy_id[1]));
  ASSERT_EQ(2, copy_hg.communityID(copy_id[3]));

  verifyPinIterators(copy_hg, { copy_hg.globalEdgeID(0), copy_hg.globalEdgeID(1), copy_hg.globalEdgeID(2) },
                     { { copy_id[1] }, { copy_id[1], copy_id[3] }, { copy_id[0], copy_id[1], copy_id[3] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, CopiesItselfWithCommunitiesAfterTwoContractions) {
  TestHypergraph hypergraph = construct_hypergraph(4, { { 1, 2 }, { 0, 1, 2, 3 }, { 2, 3 } },
                                                   { 0, 0, 1, 1 }, { 0, 1, 0 }, { 0, 1, 2, 2 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2), GLOBAL_ID(hypergraph, 3) };
  hypergraph.contract(id[1], id[2]);
  hypergraph.contract(id[1], id[3]);

  auto copy = hypergraph.copy(2, TBBArena::GLOBAL_TASK_GROUP);
  TestHypergraph& copy_hg = copy.first;
  auto& mapping = copy.second;
  std::vector<HypernodeID> copy_id = { GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[0])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[1])]),
                                       0, 0 };

  ASSERT_EQ(0, copy_hg.communityID(copy_id[0]));
  ASSERT_EQ(1, copy_hg.communityID(copy_id[1]));

  verifyPinIterators(copy_hg, { copy_hg.globalEdgeID(0), copy_hg.globalEdgeID(1), copy_hg.globalEdgeID(2) },
                     { { copy_id[1] }, { copy_id[1] }, { copy_id[0], copy_id[1] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ExtractPartZeroOfPartitionAsHypergraphWithCutNetSplitting) {
  TestHypergraph hypergraph = construct_hypergraph(4, { { 1, 2 }, { 0, 1, 2, 3 }, { 2, 3 } },
                                                   { 0, 0, 1, 1 }, { 0, 1, 0 }, { 0, 0, 1, 1 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2), GLOBAL_ID(hypergraph, 3) };
  hypergraph.setNodePart(id[0], 0);
  hypergraph.setNodePart(id[1], 1);
  hypergraph.setNodePart(id[2], 0);
  hypergraph.setNodePart(id[3], 1);
  hypergraph.initializeNumCutHyperedges();
  hypergraph.updateGlobalPartInfos();

  auto copy = hypergraph.copy(2, TBBArena::GLOBAL_TASK_GROUP, 0, true /* cut net splitting */);
  TestHypergraph& copy_hg = copy.first;
  auto& mapping = copy.second;
  std::vector<HypernodeID> copy_id = { GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[0])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[2])]) };

  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(copy_id[0]));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(copy_id[1]));

  ASSERT_EQ(0, copy_hg.communityID(copy_id[0]));
  ASSERT_EQ(1, copy_hg.communityID(copy_id[1]));

  verifyPinIterators(copy_hg, { copy_hg.globalEdgeID(0), copy_hg.globalEdgeID(1), copy_hg.globalEdgeID(2) },
                     { { copy_id[1] }, { copy_id[1] }, { copy_id[0], copy_id[1] }});
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ExtractPartOneOfPartitionAsHypergraphWithCutNetSplitting) {
  TestHypergraph hypergraph = construct_hypergraph(4, { { 1, 2 }, { 0, 1, 2, 3 }, { 2, 3 } },
                                                   { 0, 0, 1, 1 }, { 0, 1, 0 }, { 0, 0, 1, 1 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2), GLOBAL_ID(hypergraph, 3) };
  hypergraph.setNodePart(id[0], 0);
  hypergraph.setNodePart(id[1], 1);
  hypergraph.setNodePart(id[2], 0);
  hypergraph.setNodePart(id[3], 1);
  hypergraph.initializeNumCutHyperedges();
  hypergraph.updateGlobalPartInfos();

  auto copy = hypergraph.copy(2, TBBArena::GLOBAL_TASK_GROUP, 1, true /* cut net splitting */);
  TestHypergraph& copy_hg = copy.first;
  auto& mapping = copy.second;
  std::vector<HypernodeID> copy_id = { GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[1])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[3])]) };

  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(copy_id[0]));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(copy_id[1]));

  ASSERT_EQ(0, copy_hg.communityID(copy_id[0]));
  ASSERT_EQ(1, copy_hg.communityID(copy_id[1]));

  verifyPinIterators(copy_hg, { copy_hg.globalEdgeID(0), copy_hg.globalEdgeID(1), copy_hg.globalEdgeID(2) },
                     { { copy_id[0] }, { copy_id[1] }, { copy_id[0], copy_id[1] }});
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ExtractPartOneOfPartitionAsHypergraphWithCutNetSplittingAfterContraction) {
  TestHypergraph hypergraph = construct_hypergraph(4, { { 1, 2 }, { 0, 1, 2, 3 }, { 2, 3 } },
                                                   { 0, 0, 1, 1 }, { 0, 1, 0 }, { 0, 0, 1, 1 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2), GLOBAL_ID(hypergraph, 3) };
  hypergraph.contract(id[1], id[3]);

  hypergraph.setNodePart(id[0], 0);
  hypergraph.setNodePart(id[1], 1);
  hypergraph.setNodePart(id[2], 0);
  hypergraph.initializeNumCutHyperedges();
  hypergraph.updateGlobalPartInfos();

  auto copy = hypergraph.copy(2, TBBArena::GLOBAL_TASK_GROUP, 1, true /* cut net splitting */);
  TestHypergraph& copy_hg = copy.first;
  auto& mapping = copy.second;
  std::vector<HypernodeID> copy_id = { GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[1])]) };

  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(copy_id[0]));

  ASSERT_EQ(0, copy_hg.communityID(copy_id[0]));

  verifyPinIterators(copy_hg, { copy_hg.globalEdgeID(0), copy_hg.globalEdgeID(1) },
                     { { copy_id[0] }, { copy_id[0] } });
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ExtractPartZeroOfPartitionAsHypergraphWithCutNetRemoval) {
  TestHypergraph hypergraph = construct_hypergraph(4, { { 1, 2 }, { 0, 1, 2, 3 }, { 2, 3 } },
                                                   { 0, 0, 1, 1 }, { 0, 1, 0 }, { 0, 0, 1, 1 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2), GLOBAL_ID(hypergraph, 3) };
  hypergraph.setNodePart(id[0], 1);
  hypergraph.setNodePart(id[1], 1);
  hypergraph.setNodePart(id[2], 0);
  hypergraph.setNodePart(id[3], 0);
  hypergraph.initializeNumCutHyperedges();
  hypergraph.updateGlobalPartInfos();

  auto copy = hypergraph.copy(2, TBBArena::GLOBAL_TASK_GROUP, 0, false /* cut net removal */);
  TestHypergraph& copy_hg = copy.first;
  auto& mapping = copy.second;
  std::vector<HypernodeID> copy_id = { GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[2])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[3])]) };

  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(copy_id[0]));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(copy_id[1]));

  ASSERT_EQ(1, copy_hg.communityID(copy_id[0]));
  ASSERT_EQ(1, copy_hg.communityID(copy_id[1]));

  verifyPinIterators(copy_hg, { copy_hg.globalEdgeID(0) }, { { copy_id[0], copy_id[1] }});
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ExtractPartOneOfPartitionAsHypergraphWithCutNetRemoval) {
  TestHypergraph hypergraph = construct_hypergraph(4, { { 1, 2 }, { 0, 1, 2, 3 }, { 2, 3 } },
                                                   { 0, 0, 1, 1 }, { 0, 1, 0 }, { 0, 0, 1, 1 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2), GLOBAL_ID(hypergraph, 3) };
  hypergraph.setNodePart(id[0], 0);
  hypergraph.setNodePart(id[1], 1);
  hypergraph.setNodePart(id[2], 1);
  hypergraph.setNodePart(id[3], 1);
  hypergraph.initializeNumCutHyperedges();
  hypergraph.updateGlobalPartInfos();

  auto copy = hypergraph.copy(2, TBBArena::GLOBAL_TASK_GROUP, 1, false /* cut net removal */);
  TestHypergraph& copy_hg = copy.first;
  auto& mapping = copy.second;
  std::vector<HypernodeID> copy_id = { GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[1])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[2])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[3])]) };

  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(copy_id[0]));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(copy_id[1]));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(copy_id[2]));

  ASSERT_EQ(0, copy_hg.communityID(copy_id[0]));
  ASSERT_EQ(1, copy_hg.communityID(copy_id[1]));
  ASSERT_EQ(1, copy_hg.communityID(copy_id[2]));

  verifyPinIterators(copy_hg, { copy_hg.globalEdgeID(0), copy_hg.globalEdgeID(1) },
                     { { copy_id[0], copy_id[1] }, { copy_id[1], copy_id[2] }});
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ExtractPartOneOfPartitionAsHypergraphWithCutNetRemovalAfterContraction) {
  TestHypergraph hypergraph = construct_hypergraph(4, { { 1, 2 }, { 0, 1, 2, 3 }, { 2, 3 } },
                                                   { 0, 0, 1, 1 }, { 0, 1, 0 }, { 0, 0, 1, 1 });
  std::vector<HypernodeID> id = { GLOBAL_ID(hypergraph, 0), GLOBAL_ID(hypergraph, 1), GLOBAL_ID(hypergraph, 2), GLOBAL_ID(hypergraph, 3) };
  hypergraph.contract(id[2], id[3]);

  hypergraph.setNodePart(id[0], 0);
  hypergraph.setNodePart(id[1], 1);
  hypergraph.setNodePart(id[2], 1);
  hypergraph.initializeNumCutHyperedges();
  hypergraph.updateGlobalPartInfos();

  auto copy = hypergraph.copy(2, TBBArena::GLOBAL_TASK_GROUP, 1, false /* cut net removal */);
  TestHypergraph& copy_hg = copy.first;
  auto& mapping = copy.second;
  std::vector<HypernodeID> copy_id = { GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[1])]),
                                       GLOBAL_ID(copy_hg, mapping[hypergraph.originalNodeID(id[2])]) };

  ASSERT_EQ(0, TestStreamingHypergraph::get_numa_node_of_vertex(copy_id[0]));
  ASSERT_EQ(1, TestStreamingHypergraph::get_numa_node_of_vertex(copy_id[1]));

  ASSERT_EQ(0, copy_hg.communityID(copy_id[0]));
  ASSERT_EQ(1, copy_hg.communityID(copy_id[1]));

  verifyPinIterators(copy_hg, { copy_hg.globalEdgeID(0), copy_hg.globalEdgeID(1) },
                     { { copy_id[0], copy_id[1] }, { copy_id[1] }});
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ResetPartWeightAndSizeAfterResetPartitioning) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  assignPartitionIDs(hypergraph);
  hypergraph.resetPartition();

  ASSERT_EQ(0, hypergraph.partWeight(0));
  ASSERT_EQ(0, hypergraph.partSize(0));
  ASSERT_EQ(0, hypergraph.partWeight(1));
  ASSERT_EQ(0, hypergraph.partSize(1));
  ASSERT_EQ(0, hypergraph.localPartWeight(0));
  ASSERT_EQ(0, hypergraph.localPartSize(0));
  ASSERT_EQ(0, hypergraph.localPartWeight(1));
  ASSERT_EQ(0, hypergraph.localPartSize(1));
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ResetPartitionIdsAfterResetPartitioning) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  assignPartitionIDs(hypergraph);
  hypergraph.resetPartition();

  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    ASSERT_EQ(-1, hypergraph.partID(hn));
    ASSERT_EQ(0, hypergraph.numIncidentCutHyperedges(hn));
  }
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ResetConnectivitySetsAfterResetPartitioning) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  assignPartitionIDs(hypergraph);
  hypergraph.resetPartition();

  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(0, hypergraph.connectivity(he));
  }
}

TEST_F(AHypergraphWithTwoStreamingHypergraphs, ResetPinCountInPartAfterResetPartitioning) {
  TestHypergraph hypergraph = construct_test_hypergraph(*this);
  assignPartitionIDs(hypergraph);
  hypergraph.resetPartition();

  for ( const HyperedgeID& he : hypergraph.edges() ) {
    ASSERT_EQ(0, hypergraph.pinCountInPart(he, 0));
    ASSERT_EQ(0, hypergraph.pinCountInPart(he, 1));
  }
}


}  // namespace ds
}  // namespace mt_kahypar
