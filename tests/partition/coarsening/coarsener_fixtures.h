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
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/coarsening/multilevel_coarsener.h"

using ::testing::Test;
using ::testing::Eq;
using ::testing::Le;

namespace mt_kahypar {

using TestTypeTraits = ds::TestTypeTraits<2>;
using HyperGraph = typename TestTypeTraits::HyperGraph;
using HyperGraphFactory = typename TestTypeTraits::HyperGraphFactory;
using PartitionedHyperGraph = typename TestTypeTraits::template PartitionedHyperGraph<>;
using HwTopology = typename TestTypeTraits::HwTopology;
using TBB = typename TestTypeTraits::TBB;

class ACoarsener : public Test {
 private:

 public:
  ACoarsener() :
    hypergraph(HyperGraphFactory::construct(TBB::GLOBAL_TASK_GROUP,
      16, 18, { { 0, 1 }, { 0, 1, 3 }, { 1, 2, 3 }, { 2, 3, 4 }, { 2, 4 },
                { 4, 5 }, { 4, 5, 7 }, { 5, 6, 7 }, { 6, 7, 8 }, { 6, 8 },
                { 8, 9 }, { 8, 9, 11 }, { 9, 10, 11 }, { 10, 11, 12 }, { 10, 12 },
                { 12, 13 }, { 12, 13, 15 }, { 13, 14, 15 } })),
    context(),
    nullptr_refiner(nullptr) {
    for ( const HypernodeID& hn : hypergraph.nodes() ) {
      hypergraph.setCommunityID(hn, hypergraph.originalNodeID(hn) / 4);
    }
    hypergraph.initializeCommunities(TBB::GLOBAL_TASK_GROUP);
    hypergraph.setCommunityNodeMapping({ 0, 0, 1, 1 });

    context.partition.k = 2;
    context.partition.objective = kahypar::Objective::km1;
    context.coarsening.max_allowed_node_weight = std::numeric_limits<HypernodeWeight>::max();
    context.coarsening.contraction_limit = 8;
    context.coarsening.use_adaptive_max_allowed_node_weight = false;
    context.coarsening.minimum_shrink_factor = 1.0;
    context.coarsening.maximum_shrink_factor = 4.0;
    context.setupPartWeights(hypergraph.totalWeight());
  }

  static void SetUpTestSuite() {
    TBB::instance(HwTopology::instance().num_cpus());
  }

  HyperGraph hypergraph;
  Context context;
  std::unique_ptr<IRefinerT<TestTypeTraits>> nullptr_refiner;
};

template <class HyperGraph>
void assignPartitionIDs(HyperGraph& hypergraph) {
  for (const HypernodeID& hn : hypergraph.nodes()) {
    PartitionID part_id = common::get_numa_node_of_vertex(hn);
    hypergraph.setNodePart(hn, part_id);
  }
  hypergraph.initializeNumCutHyperedges(TBB::GLOBAL_TASK_GROUP);
}

template <class HyperGraph>
HypernodeID currentNumNodes(HyperGraph& hypergraph) {
  HypernodeID num_nodes = 0;
  for (const HypernodeID& hn : hypergraph.nodes()) {
    unused(hn);
    ++num_nodes;
  }
  return num_nodes;
}

template <class HyperGraph>
HyperedgeID currentNumEdges(HyperGraph& hypergraph) {
  HyperedgeID num_edges = 0;
  for (const HyperedgeID& he : hypergraph.edges()) {
    unused(he);
    ++num_edges;
  }
  return num_edges;
}

template <class HyperGraph>
HypernodeID currentNumPins(HyperGraph& hypergraph) {
  HypernodeID num_pins = 0;
  for (const HypernodeID& he : hypergraph.edges()) {
    num_pins += hypergraph.edgeSize(he);
  }
  return num_pins;
}

template <class Coarsener>
void doCoarsening(Coarsener& coarsener) {
  coarsener.disableRandomization();
  coarsener.coarsen();
}

template <class Coarsener>
void decreasesNumberOfPins(Coarsener& coarsener,
                           const size_t number_of_pins) {
  doCoarsening(coarsener);
  ASSERT_THAT(currentNumPins(coarsener.coarsestHypergraph()), Eq(number_of_pins));
}

template <class Coarsener>
void decreasesNumberOfHyperedges(Coarsener& coarsener,
                                 const HyperedgeID num_hyperedges) {
  doCoarsening(coarsener);
  ASSERT_THAT(currentNumEdges(coarsener.coarsestHypergraph()), Eq(num_hyperedges));
}

template <class Coarsener, class HyperGraph>
void removesHyperedgesOfSizeOneDuringCoarsening(Coarsener& coarsener,
                                                HyperGraph& hypergraph,
                                                const std::vector<HyperedgeID>& single_node_hes) {
  doCoarsening(coarsener);
  for (const HyperedgeID& he : single_node_hes) {
    ASSERT_THAT(hypergraph.edgeIsEnabled(he), Eq(false)) << V(he);
  }
}

template <class Coarsener, class HyperGraph, class TypeTraits>
void reAddsHyperedgesOfSizeOneDuringUncoarsening(Coarsener& coarsener,
                                                 std::unique_ptr<IRefinerT<TypeTraits>>& refiner,
                                                 HyperGraph& hypergraph,
                                                 const std::vector<HyperedgeID>& single_node_hes) {
  doCoarsening(coarsener);
  for (const HyperedgeID& he : single_node_hes) {
    ASSERT_THAT(hypergraph.edgeIsEnabled(he), Eq(false)) << V(he);
  }
  PartitionedHyperGraph& partitioned_hypergraph = coarsener.coarsestPartitionedHypergraph();
  assignPartitionIDs(partitioned_hypergraph);
  coarsener.uncoarsen(refiner, refiner);
  for (const HyperedgeID& he : single_node_hes) {
    ASSERT_THAT(hypergraph.edgeIsEnabled(he), Eq(true)) << V(he);
  }
}

template <class Coarsener, class HyperGraph>
void removesParallelHyperedgesDuringCoarsening(Coarsener& coarsener,
                                               HyperGraph& hypergraph,
                                               const std::vector<HyperedgeID>& parallel_hes) {
  doCoarsening(coarsener);
  for (const HyperedgeID& he : parallel_hes) {
    ASSERT_THAT(hypergraph.edgeIsEnabled(he), Eq(false)) << V(he);
  }
}

template <class Coarsener, class HyperGraph>
void updatesEdgeWeightOfRepresentativeHyperedgeOnParallelHyperedgeRemoval(Coarsener& coarsener,
                                                                          HyperGraph& hypergraph,
                                                                          const std::vector<std::pair<HyperedgeID, HyperedgeWeight> >& he_weights) {
  doCoarsening(coarsener);
  for (const auto& he_weight : he_weights) {
    HyperedgeID he = he_weight.first;
    HyperedgeWeight weight = he_weight.second;
    ASSERT_THAT(hypergraph.edgeIsEnabled(he), Eq(true)) << V(he);
    ASSERT_THAT(hypergraph.edgeWeight(he), Eq(weight));
  }
}

template <class Coarsener, class HyperGraph, class TypeTraits>
void restoresParallelHyperedgesDuringUncoarsening(Coarsener& coarsener,
                                                  std::unique_ptr<IRefinerT<TypeTraits>>& refiner,
                                                  HyperGraph& hypergraph,
                                                  const std::vector<HyperedgeID>& parallel_hes) {
  doCoarsening(coarsener);
  for (const HyperedgeID& he : parallel_hes) {
    ASSERT_THAT(hypergraph.edgeIsEnabled(he), Eq(false)) << V(he);
  }
  PartitionedHyperGraph& partitioned_hypergraph = coarsener.coarsestPartitionedHypergraph();
  assignPartitionIDs(partitioned_hypergraph);
  coarsener.uncoarsen(refiner, refiner);
  for (const HyperedgeID& he : parallel_hes) {
    ASSERT_THAT(hypergraph.edgeIsEnabled(he), Eq(true)) << V(he);
  }
}

template <class Coarsener, class HyperGraph>
void doesNotCoarsenUntilCoarseningLimit(Coarsener& coarsener,
                                        HyperGraph& hypergraph,
                                        Context& context,
                                        const HypernodeID contraction_limit,
                                        const HypernodeWeight max_allowed_node_weight,
                                        const size_t expected_num_nodes) {
  context.coarsening.contraction_limit = contraction_limit;
  context.coarsening.max_allowed_node_weight = max_allowed_node_weight;
  doCoarsening(coarsener);
  for (const HypernodeID& hn : hypergraph.nodes()) {
    ASSERT_THAT(hypergraph.nodeWeight(hn), Le(context.coarsening.max_allowed_node_weight));
  }
  ASSERT_THAT(currentNumNodes(hypergraph), Eq(expected_num_nodes));
}
}  // namespace mt_kahypar
