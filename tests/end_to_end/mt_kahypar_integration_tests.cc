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

#include "mt-kahypar/application/command_line_options.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/io/hypergraph_io.h"
#include "mt-kahypar/mt_kahypar.h"
#include "mt-kahypar/partition/context.h"
#include "mt-kahypar/partition/partitioner.h"

using ::testing::Test;

namespace mt_kahypar {

template< PartitionID k,
          kahypar::Objective objective,
          CommunityLoadBalancingStrategy balancing_strategy,
          InitialPartitioningMode initial_partitioning_mode,
          LabelPropagationAlgorithm lp_algorithm,
          bool use_batch_uncontractions >
struct TestConfig {
  static constexpr PartitionID K = k;
  static constexpr kahypar::Objective OBJECTIVE = objective;
  static constexpr CommunityLoadBalancingStrategy BALANCING_STRATEGY = balancing_strategy;
  static constexpr InitialPartitioningMode INITIAL_PARTITIONING_MODE = initial_partitioning_mode;
  static constexpr LabelPropagationAlgorithm LP_ALGORITHM = lp_algorithm;
  static constexpr bool USE_BATCH_UNCONTRACTIONS = use_batch_uncontractions;
};

template< typename Config >
class MtKaHyPar : public Test {
  static size_t num_threads;

 public:
  MtKaHyPar() :
    context() {
    parseIniToContext(context, "../../../config/shared_memory_context.ini");
    context.partition.graph_filename = "test_instances/ibm01.hgr";
    context.partition.k = Config::K;
    context.partition.objective = Config::OBJECTIVE;
    context.partition.epsilon = 0.03;
    context.partition.seed = 42;
    context.partition.verbose_output = true;
    context.partition.detailed_timings = true;
    context.preprocessing.community_detection.load_balancing_strategy = Config::BALANCING_STRATEGY;
    context.initial_partitioning.mode = Config::INITIAL_PARTITIONING_MODE;
    context.initial_partitioning.context_file = "../../../config/fast_initial_partitioning.ini";
    context.refinement.use_batch_uncontractions = Config::USE_BATCH_UNCONTRACTIONS;
    context.refinement.label_propagation.algorithm = Config::LP_ALGORITHM;
    context.shared_memory.num_threads = num_threads;
  }

  static void SetUpTestSuite() {
    TBBNumaArena::instance(num_threads);
  }

  Context context;
};

template< typename Config >
size_t MtKaHyPar<Config>::num_threads = HardwareTopology::instance().num_cpus();

void verifyThatHypergraphsAreEquivalent(const Hypergraph& hypergraph,
                                        const Hypergraph& reference) {
  // Verify equivallence of hypernodes and incident nets
  for (const HypernodeID& hn : reference.nodes()) {
    const HypernodeID original_id = reference.originalNodeID(hn);
    const HypernodeID u = hypergraph.globalNodeID(original_id);
    ASSERT_TRUE(hypergraph.nodeIsEnabled(u));

    std::set<HyperedgeID> incident_nets;
    for (const HyperedgeID& he : reference.incidentEdges(hn)) {
      const HyperedgeID original_edge_id = reference.originalEdgeID(he);
      incident_nets.insert(original_edge_id);
    }

    size_t num_incident_nets = 0;
    for (const HyperedgeID& he : hypergraph.incidentEdges(u)) {
      const HyperedgeID original_edge_id = hypergraph.originalEdgeID(he);
      ASSERT_TRUE(incident_nets.find(original_edge_id) != incident_nets.end()) << V(u) << V(original_edge_id);
      ++num_incident_nets;
    }
    ASSERT_EQ(num_incident_nets, incident_nets.size());
  }

  // Verify equivallence of hyperedges and pins
  for (const HyperedgeID& he : reference.edges()) {
    const HyperedgeID original_id = reference.originalEdgeID(he);
    const HyperedgeID e = hypergraph.globalEdgeID(original_id);
    ASSERT_TRUE(hypergraph.edgeIsEnabled(e));

    std::set<HypernodeID> pins;
    for (const HypernodeID& pin : reference.pins(he)) {
      const HypernodeID original_pin_id = reference.originalNodeID(pin);
      pins.insert(original_pin_id);
    }

    size_t num_pins = 0;
    for (const HypernodeID& pin : hypergraph.pins(e)) {
      const HypernodeID original_pin_id = hypergraph.originalNodeID(pin);
      ASSERT_TRUE(pins.find(original_pin_id) != pins.end()) << V(e) << V(original_pin_id);
      ++num_pins;
    }
    ASSERT_EQ(num_pins, pins.size());
  }
}

template< PartitionID k,
          CommunityLoadBalancingStrategy balancing_strategy,
          InitialPartitioningMode initial_partitioning_mode,
          bool use_batch_uncontractions>
using Km1Config = TestConfig<k,
                             kahypar::Objective::km1,
                             balancing_strategy,
                             initial_partitioning_mode,
                             LabelPropagationAlgorithm::label_propagation_km1,
                             use_batch_uncontractions>;

template< PartitionID k,
          CommunityLoadBalancingStrategy balancing_strategy,
          InitialPartitioningMode initial_partitioning_mode,
          bool use_batch_uncontractions>
using CutConfig = TestConfig<k,
                             kahypar::Objective::cut,
                             balancing_strategy,
                             initial_partitioning_mode,
                             LabelPropagationAlgorithm::label_propagation_cut,
                             use_batch_uncontractions>;

typedef ::testing::Types<CutConfig<2,
                                   CommunityLoadBalancingStrategy::none,
                                   InitialPartitioningMode::direct,
                                   false>,
                         CutConfig<4,
                                   CommunityLoadBalancingStrategy::none,
                                   InitialPartitioningMode::direct,
                                   false>,
                         CutConfig<8,
                                   CommunityLoadBalancingStrategy::none,
                                   InitialPartitioningMode::direct,
                                   false>,
                         CutConfig<2,
                                   CommunityLoadBalancingStrategy::none,
                                   InitialPartitioningMode::recursive,
                                   false>,
                         CutConfig<4,
                                   CommunityLoadBalancingStrategy::none,
                                   InitialPartitioningMode::recursive,
                                   false>,
                         CutConfig<8,
                                   CommunityLoadBalancingStrategy::none,
                                   InitialPartitioningMode::recursive,
                                   false>,
                         Km1Config<2,
                                   CommunityLoadBalancingStrategy::none,
                                   InitialPartitioningMode::direct,
                                   false>,
                         Km1Config<4,
                                   CommunityLoadBalancingStrategy::none,
                                   InitialPartitioningMode::direct,
                                   false>,
                         Km1Config<8,
                                   CommunityLoadBalancingStrategy::none,
                                   InitialPartitioningMode::direct,
                                   false>,
                         Km1Config<2,
                                   CommunityLoadBalancingStrategy::none,
                                   InitialPartitioningMode::recursive,
                                   false>,
                         Km1Config<4,
                                   CommunityLoadBalancingStrategy::none,
                                   InitialPartitioningMode::recursive,
                                   false>,
                         Km1Config<8,
                                   CommunityLoadBalancingStrategy::none,
                                   InitialPartitioningMode::recursive,
                                   false>,
                         Km1Config<2,
                                   CommunityLoadBalancingStrategy::none,
                                   InitialPartitioningMode::recursive,
                                   true>,
                         Km1Config<4,
                                   CommunityLoadBalancingStrategy::none,
                                   InitialPartitioningMode::recursive,
                                   true>,
                         Km1Config<8,
                                   CommunityLoadBalancingStrategy::none,
                                   InitialPartitioningMode::recursive,
                                   true>,
                         Km1Config<2,
                                   CommunityLoadBalancingStrategy::size_constraint,
                                   InitialPartitioningMode::direct,
                                   false> > TestConfigs;

TYPED_TEST_CASE(MtKaHyPar, TestConfigs);

TYPED_TEST(MtKaHyPar, PartitionsAHypergraph) {
  // Partition Hypergraph
  Hypergraph hypergraph = io::readHypergraphFile(
    this->context.partition.graph_filename, this->context.partition.k,
    this->context.shared_memory.initial_distribution);

  HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
  partition::Partitioner().partition(hypergraph, this->context);
  HighResClockTimepoint end = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed_seconds(end - start);
  mt_kahypar::io::printPartitioningResults(hypergraph, this->context, elapsed_seconds);

  // Verify that partitioned hypergraph is
  // equivalent with input hypergraph
  Hypergraph reference = io::readHypergraphFile(
    this->context.partition.graph_filename, this->context.partition.k,
    this->context.shared_memory.initial_distribution);
  verifyThatHypergraphsAreEquivalent(hypergraph, reference);

  HypernodeID num_hypernodes = 0;
  std::vector<HypernodeWeight> weights(hypergraph.k(), 0);
  std::vector<size_t> sizes(hypergraph.k(), 0);
  for ( const HypernodeID& hn : hypergraph.nodes() ) {
    // Make sure that each hypernode has a block id between [0,k)
    PartitionID part_id = hypergraph.partID(hn);
    ASSERT_NE(-1, part_id);
    ASSERT_LE(part_id, this->context.partition.k);
    weights[part_id] += hypergraph.nodeWeight(hn);
    ++sizes[part_id];
    ++num_hypernodes;

    // Verify that border nodes are stored correctly
    bool is_border_node = false;
    for ( const HyperedgeID& he : hypergraph.incidentEdges(hn) ) {
      if ( hypergraph.connectivity(he) > 1 ) {
        is_border_node = true;
        break;
      }
    }
    ASSERT_EQ(is_border_node, hypergraph.isBorderNode(hn));
  }
  // Verify that hypergraph is fully uncontracted
  ASSERT_EQ(hypergraph.initialNumNodes(), num_hypernodes);

  // Verify block weights and sizes
  for ( PartitionID k = 0; k < hypergraph.k(); ++k ) {
    ASSERT_EQ(weights[k], hypergraph.partWeight(k));
    ASSERT_EQ(sizes[k], hypergraph.partSize(k));
  }

  // Verify connectivity (sets) and pin count in part
  size_t num_hyperedges = 0;
  for ( const HyperedgeID& he : hypergraph.edges() ) {
    std::vector<HypernodeID> pin_count_in_part(hypergraph.k(), 0);
    PartitionID connectivity = 0;
    for ( const HypernodeID& pin : hypergraph.pins(he) ) {
      PartitionID part_id = hypergraph.partID(pin);
      if ( pin_count_in_part[part_id] == 0 ) {
        ++connectivity;
      }
      ++pin_count_in_part[part_id];
    }
    ASSERT_EQ(connectivity, hypergraph.connectivity(he));

    for ( PartitionID k = 0; k < hypergraph.k(); ++k ) {
      ASSERT_EQ(pin_count_in_part[k], hypergraph.pinCountInPart(he, k));
    }

    PartitionID connectivity_2 = 0;
    for ( const PartitionID& k : hypergraph.connectivitySet(he) ) {
      ++connectivity_2;
      ASSERT_GT(pin_count_in_part[k], 0);
    }
    ASSERT_EQ(connectivity, connectivity_2);

    ++num_hyperedges;
  }
  ASSERT_EQ(hypergraph.initialNumEdges(), num_hyperedges);
}


}  // namespace mt_kahypar
