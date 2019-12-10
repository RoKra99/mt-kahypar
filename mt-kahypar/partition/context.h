/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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
#pragma once

#include "kahypar/definitions.h"
#include "kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/context_enum_classes.h"

namespace mt_kahypar {
struct PartitioningParameters {
  kahypar::Mode mode = kahypar::Mode::UNDEFINED;
  kahypar::Objective objective = kahypar::Objective::UNDEFINED;
  double epsilon = std::numeric_limits<double>::max();
  PartitionID k = std::numeric_limits<PartitionID>::max();
  int seed = 0;

  int time_limit = 0;
  std::vector<HypernodeWeight> perfect_balance_part_weights;
  std::vector<HypernodeWeight> max_part_weights;
  HyperedgeID hyperedge_size_threshold = 1000;

  bool verbose_output = false;
  bool quiet_mode = false;
  bool detailed_timings = false;
  bool sp_process_output = false;
  bool write_partition_file = false;

  std::string graph_filename { };
  std::string graph_partition_filename { };
  std::string graph_community_filename { };
};

inline std::ostream & operator<< (std::ostream& str, const PartitioningParameters& params) {
  str << "Partitioning Parameters:" << std::endl;
  str << "  Hypergraph:                         " << params.graph_filename << std::endl;
  str << "  Partition File:                     " << params.graph_partition_filename << std::endl;
  str << "  Community File:                     " << params.graph_community_filename << std::endl;
  str << "  Mode:                               " << params.mode << std::endl;
  str << "  Objective:                          " << params.objective << std::endl;
  str << "  k:                                  " << params.k << std::endl;
  str << "  epsilon:                            " << params.epsilon << std::endl;
  str << "  seed:                               " << params.seed << std::endl;
  str << "  time limit:                         " << params.time_limit << "s" << std::endl;
  str << "  hyperedge size threshold:           " << params.hyperedge_size_threshold << std::endl;
  return str;
}

struct CommunityDetectionParameters {
  CommunityLoadBalancingStrategy load_balancing_strategy = CommunityLoadBalancingStrategy::none;
  size_t size_constraint_factor = 0;
  LouvainEdgeWeight edge_weight_function = LouvainEdgeWeight::UNDEFINED;
  uint32_t max_pass_iterations = std::numeric_limits<uint32_t>::max();
  long double min_eps_improvement = std::numeric_limits<long double>::max();
};

inline std::ostream & operator<< (std::ostream& str, const CommunityDetectionParameters& params) {
  str << "  Community Detection Parameters:" << std::endl;
  str << "    Load Balancing Strategy:          " << params.load_balancing_strategy << std::endl;
  if (params.load_balancing_strategy == CommunityLoadBalancingStrategy::size_constraint) {
    str << "    Size Constraint Factor:           " << params.size_constraint_factor << std::endl;
  }
  str << "    Edge Weight Function:             " << params.edge_weight_function << std::endl;
  str << "    Maximum Louvain-Pass Iterations:  " << params.max_pass_iterations << std::endl;
  str << "    Minimum Quality Improvement:      " << params.min_eps_improvement << std::endl;
  return str;
}

struct CommunityRedistributionParameters {
  bool use_community_redistribution = false;
  CommunityAssignmentObjective assignment_objective = CommunityAssignmentObjective::UNDEFINED;
  CommunityAssignmentStrategy assignment_strategy = CommunityAssignmentStrategy::UNDEFINED;
};

inline std::ostream & operator<< (std::ostream& str, const CommunityRedistributionParameters& params) {
  str << "  Community Detection Parameters:" << std::endl;
  str << "    Use Community Redistribution:     " << std::boolalpha << params.use_community_redistribution << std::endl;
  str << "    Community Assignment Objective:   " << params.assignment_objective << std::endl;
  str << "    Community Assignment Strategy:    " << params.assignment_strategy << std::endl;
  return str;
}

struct PreprocessingParameters {
  bool use_community_structure_from_file = false;
  CommunityDetectionParameters community_detection = { };
  CommunityRedistributionParameters community_redistribution = { };
};

inline std::ostream & operator<< (std::ostream& str, const PreprocessingParameters& params) {
  str << "Preprocessing Parameters:" << std::endl;
  str << "  Use Community Structure from File:  " << std::boolalpha << params.use_community_structure_from_file << std::endl;
  if (!params.use_community_structure_from_file) {
    str << std::endl << params.community_detection;
  }
  str << std::endl << params.community_redistribution;
  return str;
}

struct RatingParameters {
  RatingFunction rating_function = RatingFunction::UNDEFINED;
  HeavyNodePenaltyPolicy heavy_node_penalty_policy = HeavyNodePenaltyPolicy::UNDEFINED;
  AcceptancePolicy acceptance_policy = AcceptancePolicy::UNDEFINED;
};

inline std::ostream & operator<< (std::ostream& str, const RatingParameters& params) {
  str << "  Rating Parameters:" << std::endl;
  str << "    Rating Function:                  " << params.rating_function << std::endl;
  str << "    Heavy Node Penalty:               " << params.heavy_node_penalty_policy << std::endl;
  str << "    Acceptance Policy:                " << params.acceptance_policy << std::endl;
  return str;
}

struct CoarseningParameters {
  CoarseningAlgorithm algorithm = CoarseningAlgorithm::UNDEFINED;
  RatingParameters rating = { };
  HypernodeID contraction_limit_multiplier = std::numeric_limits<HypernodeID>::max();
  double max_allowed_weight_multiplier = std::numeric_limits<double>::max();
  bool use_high_degree_vertex_threshold = false;

  // Those will be determined dynamically
  HypernodeWeight max_allowed_node_weight = 0;
  HypernodeID contraction_limit = 0;
  double hypernode_weight_fraction = 0.0;
  HyperedgeID high_degree_vertex_threshold = std::numeric_limits<HyperedgeID>::max();
};

inline std::ostream & operator<< (std::ostream& str, const CoarseningParameters& params) {
  str << "Coarsening Parameters:" << std::endl;
  str << "  Algorithm:                          " << params.algorithm << std::endl;
  str << "  max allowed weight multiplier:      " << params.max_allowed_weight_multiplier << std::endl;
  str << "  maximum allowed hypernode weight:   " << params.max_allowed_node_weight << std::endl;
  str << "  contraction limit multiplier:       " << params.contraction_limit_multiplier << std::endl;
  str << "  contraction limit:                  " << params.contraction_limit << std::endl;
  if ( params.use_high_degree_vertex_threshold ) {
    str << "  high degree degree threshold:       " << params.high_degree_vertex_threshold << std::endl;
  }
  str << std::endl << params.rating;
  return str;
}

struct InitialPartitioningParameters {
  std::string context_file = "";
  InitialPartitioningMode mode = InitialPartitioningMode::UNDEFINED;
  bool call_kahypar_multiple_times = false;
  size_t runs = 1;
};

inline std::ostream & operator<< (std::ostream& str, const InitialPartitioningParameters& params) {
  str << "Initial Partitioning Parameters:" << std::endl;
  str << "  Initial Partitioning Context:       " << params.context_file << std::endl;
  str << "  Initial Partitioning Mode:          " << params.mode << std::endl;
  str << "  Call KaHyPar multiple times:        " << std::boolalpha << params.call_kahypar_multiple_times << std::endl;
  str << "  Number of Runs:                     " << params.runs << std::endl;
  return str;
}

struct RefinementParameters {
  struct FlowParameter{
    double alpha = 16.0;
    bool use_improvement_history = true;
    bool use_most_balanced_minimum_cut = true;
  };
  RefinementAlgorithm algorithm = RefinementAlgorithm::do_nothing;
  size_t maximum_iterations = 1;
  size_t part_weight_update_frequency = 100;
  bool numa_aware = false;
  bool rebalancing = true;
  bool use_batch_uncontractions = true;
  size_t batch_size = 1;
  ExecutionType execution_policy = ExecutionType::UNDEFINED;
  double execution_policy_alpha = 2.0;
  FlowParameter flow{};
};

inline std::ostream& operator<< (std::ostream& str, const RefinementParameters& params) {
  str << "Refinement Parameters:" << std::endl;
  if(params.algorithm == RefinementAlgorithm::flow){
    str << "  Flow Parameters:" << std::endl;
    str << "    alpha:                            " << params.flow.alpha << std::endl;
    str << "    use improvement history:          "
        << std::boolalpha << params.flow.use_improvement_history << std::endl;
        str << "    use most balanced minimum cut:"
        << std::boolalpha << params.flow.use_most_balanced_minimum_cut << std::endl;

  }else{
    str << "  Label Propagation Parameters:" << std::endl;
  }
  str << "    Algorithm:                        " << params.algorithm << std::endl;
  str << "    Maximum Iterations:               " << params.maximum_iterations << std::endl;
  str << "    Part Weight Update Frequency:     " << params.part_weight_update_frequency << std::endl;
  str << "    Numa Aware:                       " << std::boolalpha << params.numa_aware << std::endl;
  str << "    Rebalancing:                      " << std::boolalpha << params.rebalancing << std::endl;
  str << "    Execution Policy:                 " << params.execution_policy << std::endl;
  str << "    Execution Policy Alpha:           " << params.execution_policy_alpha << std::endl;
  return str;
}


struct SharedMemoryParameters {
  size_t num_threads = 1;
  InitialHyperedgeDistribution initial_distribution = InitialHyperedgeDistribution::UNDEFINED;
};

inline std::ostream & operator<< (std::ostream& str, const SharedMemoryParameters& params) {
  str << "Shared Memory Parameters:             " << std::endl;
  str << "  Number of Threads:                  " << params.num_threads << std::endl;
  str << "  Initial Hyperedge Distribution:     " << params.initial_distribution << std::endl;
  return str;
}

class Context {
 public:
  PartitioningParameters partition { };
  PreprocessingParameters preprocessing { };
  CoarseningParameters coarsening { };
  InitialPartitioningParameters initial_partitioning { };
  RefinementParameters refinement { };
  SharedMemoryParameters shared_memory { };
  kahypar::ContextType type = kahypar::ContextType::main;

  Context() { }

  bool isMainRecursiveBisection() const {
    return partition.mode == kahypar::Mode::recursive_bisection &&
           type == kahypar::ContextType::main;
  }

  void setupPartWeights(const HypernodeWeight total_hypergraph_weight) {
    partition.perfect_balance_part_weights.clear();
    partition.perfect_balance_part_weights.push_back(ceil(
                                                       total_hypergraph_weight
                                                       / static_cast<double>(partition.k)));
    for (PartitionID part = 1; part != partition.k; ++part) {
      partition.perfect_balance_part_weights.push_back(
        partition.perfect_balance_part_weights[0]);
    }
    partition.max_part_weights.clear();
    partition.max_part_weights.push_back((1 + partition.epsilon)
                                         * partition.perfect_balance_part_weights[0]);
    for (PartitionID part = 1; part != partition.k; ++part) {
      partition.max_part_weights.push_back(partition.max_part_weights[0]);
    }
  }

  void sanityCheck() {
    if (partition.objective == kahypar::Objective::cut &&
        refinement.algorithm == RefinementAlgorithm::label_propagation_km1) {
      ALGO_SWITCH("Refinement algorithm" << refinement.algorithm << "only works for km1 metric."
                                         << "Do you want to use the cut version of the label propagation refiner (Y/N)?",
                  "Partitioning with" << refinement.algorithm
                                         << "refiner in combination with cut metric is not possible!",
                  refinement.algorithm,
                  RefinementAlgorithm::label_propagation_cut);
    } else if (partition.objective == kahypar::Objective::km1 &&
               refinement.algorithm == RefinementAlgorithm::label_propagation_cut) {
      ALGO_SWITCH("Refinement algorithm" << refinement.algorithm << "only works for cut metric."
                                         << "Do you want to use the km1 version of the label propagation refiner (Y/N)?",
                  "Partitioning with" << refinement.algorithm
                                         << "refiner in combination with km1 metric is not possible!",
                  refinement.algorithm,
                  RefinementAlgorithm::label_propagation_km1);
    }
  }
};

inline std::ostream & operator<< (std::ostream& str, const Context& context) {
  str << "*******************************************************************************\n"
      << "*                            Partitioning Context                             *\n"
      << "*******************************************************************************\n"
      << context.partition
      << "-------------------------------------------------------------------------------\n"
      << context.preprocessing
      << "-------------------------------------------------------------------------------\n"
      << context.coarsening
      << "-------------------------------------------------------------------------------\n"
      << context.initial_partitioning
      << "-------------------------------------------------------------------------------\n"
      << context.refinement
      << "-------------------------------------------------------------------------------\n"
      << context.shared_memory
      << "-------------------------------------------------------------------------------";
  return str;
}
}  // namespace mt_kahypar
