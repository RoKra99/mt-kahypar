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
  Paradigm paradigm = Paradigm::multilevel;
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
  bool show_detailed_timings = false;
  bool show_detailed_clustering_timings = false;
  bool show_memory_consumption = false;
  bool enable_progress_bar = false;
  bool sp_process_output = false;
  bool write_partition_file = false;

  bool enable_profiler = false;
  int snapshot_interval = std::numeric_limits<int>::max();

  std::string graph_filename { };
  std::string graph_partition_filename { };
  std::string graph_community_filename { };
};

inline std::ostream & operator<< (std::ostream& str, const PartitioningParameters& params) {
  str << "Partitioning Parameters:" << std::endl;
  str << "  Hypergraph:                         " << params.graph_filename << std::endl;
  str << "  Partition File:                     " << params.graph_partition_filename << std::endl;
  str << "  Community File:                     " << params.graph_community_filename << std::endl;
  str << "  Paradigm:                           " << params.paradigm << std::endl;
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
  LouvainEdgeWeight edge_weight_function = LouvainEdgeWeight::UNDEFINED;
  uint32_t max_pass_iterations = std::numeric_limits<uint32_t>::max();
  long double min_eps_improvement = std::numeric_limits<long double>::max();
  size_t vertex_degree_sampling_threshold = std::numeric_limits<size_t>::max();
};

inline std::ostream & operator<< (std::ostream& str, const CommunityDetectionParameters& params) {
  str << "  Community Detection Parameters:" << std::endl;
  str << "    Edge Weight Function:             " << params.edge_weight_function << std::endl;
  str << "    Maximum Louvain-Pass Iterations:  " << params.max_pass_iterations << std::endl;
  str << "    Minimum Quality Improvement:      " << params.min_eps_improvement << std::endl;
  str << "    Vertex Degree Sampling Threshold: " << params.vertex_degree_sampling_threshold << std::endl;
  return str;
}

struct CommunityRedistributionParameters {
  CommunityAssignmentObjective assignment_objective = CommunityAssignmentObjective::UNDEFINED;
  CommunityAssignmentStrategy assignment_strategy = CommunityAssignmentStrategy::UNDEFINED;
};

inline std::ostream & operator<< (std::ostream& str, const CommunityRedistributionParameters& params) {
  str << "  Community Detection Parameters:" << std::endl;
  str << "    Community Assignment Objective:   " << params.assignment_objective << std::endl;
  str << "    Community Assignment Strategy:    " << params.assignment_strategy << std::endl;
  return str;
}

struct PreprocessingParameters {
  bool use_community_detection = false;
  bool use_community_redistribution = false;
  CommunityDetectionParameters community_detection = { };
  CommunityRedistributionParameters community_redistribution = { };
};

inline std::ostream & operator<< (std::ostream& str, const PreprocessingParameters& params) {
  str << "Preprocessing Parameters:" << std::endl;
  str << "  Use Community Detection:            " << std::boolalpha << params.use_community_detection << std::endl;
  str << "  Use Community Redistribution:       " << std::boolalpha << params.use_community_redistribution << std::endl;
  if (params.use_community_detection) {
    str << std::endl << params.community_detection;
  }
  if ( params.use_community_redistribution ) {
    str << std::endl << params.community_redistribution;
  }
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
  bool use_adaptive_max_allowed_node_weight = false;
  double max_allowed_weight_fraction = std::numeric_limits<double>::max();
  double adaptive_node_weight_shrink_factor_threshold = std::numeric_limits<double>::max();
  double max_allowed_weight_multiplier = std::numeric_limits<double>::max();
  double minimum_shrink_factor = std::numeric_limits<double>::max();
  double maximum_shrink_factor = std::numeric_limits<double>::max();
  size_t vertex_degree_sampling_threshold = std::numeric_limits<size_t>::max();

  // Those will be determined dynamically
  HypernodeWeight max_allowed_node_weight = 0;
  HypernodeID contraction_limit = 0;
};

inline std::ostream & operator<< (std::ostream& str, const CoarseningParameters& params) {
  str << "Coarsening Parameters:" << std::endl;
  str << "  Algorithm:                          " << params.algorithm << std::endl;
  str << "  use adaptive max node weight:       " << std::boolalpha << params.use_adaptive_max_allowed_node_weight << std::endl;
  if ( params.use_adaptive_max_allowed_node_weight ) {
    str << "  max allowed weight fraction:        " << params.max_allowed_weight_fraction << std::endl;
    str << "  adaptive node weight threshold:     " << params.adaptive_node_weight_shrink_factor_threshold << std::endl;
    str << "  initial max hypernode weight:       " << params.max_allowed_node_weight << std::endl;
  } else {
    str << "  max allowed weight multiplier:      " << params.max_allowed_weight_multiplier << std::endl;
    str << "  maximum allowed hypernode weight:   " << params.max_allowed_node_weight << std::endl;
  }
  str << "  contraction limit multiplier:       " << params.contraction_limit_multiplier << std::endl;
  str << "  contraction limit:                  " << params.contraction_limit << std::endl;
  if ( params.algorithm == CoarseningAlgorithm::multilevel_coarsener ) {
    str << "  minimum shrink factor:              " << params.minimum_shrink_factor << std::endl;
    str << "  maximum shrink factor:              " << params.maximum_shrink_factor << std::endl;
  }
  str << "  vertex degree sampling threshold:   " << params.vertex_degree_sampling_threshold << std::endl;
  str << std::endl << params.rating;
  return str;
}

struct LabelPropagationParameters {
  LabelPropagationAlgorithm algorithm = LabelPropagationAlgorithm::do_nothing;
  size_t maximum_iterations = 1;
  bool numa_aware = false;
  bool rebalancing = true;
  bool execute_sequential = false;
  size_t hyperedge_size_activation_threshold = std::numeric_limits<size_t>::max();
};

inline std::ostream & operator<< (std::ostream& str, const LabelPropagationParameters& params) {
  str << "  Label Propagation Parameters:" << std::endl;
  str << "    Algorithm:                        " << params.algorithm << std::endl;
  str << "    Maximum Iterations:               " << params.maximum_iterations << std::endl;
  str << "    Numa Aware:                       " << std::boolalpha << params.numa_aware << std::endl;
  str << "    Rebalancing:                      " << std::boolalpha << params.rebalancing << std::endl;
  str << "    HE Size Activation Threshold:     " << std::boolalpha << params.hyperedge_size_activation_threshold << std::endl;
  return str;
}

struct FlowParameters {
  FlowAlgorithm algorithm = FlowAlgorithm::do_nothing;
  double alpha = 16.0;
  bool use_most_balanced_minimum_cut = true;
  bool use_improvement_history = true;
  double scheduling_thresh = 3.0;
};

inline std::ostream & operator<< (std::ostream& str, const FlowParameters& params) {
  str << "  Flow Parameters:" << std::endl;
  str << "    Algorithm:                        " << params.algorithm << std::endl;
  str << "    Alpha:                            " << params.alpha << std::endl;
  str << "    Use Most Balanced Minimum Cut:    " << std::boolalpha << params.use_most_balanced_minimum_cut << std::endl;
  str << "    Use Improvement History:          " << std::boolalpha << params.use_improvement_history << std::endl;
  if(params.algorithm == FlowAlgorithm::flow_opt){
    str << "    Scheduling Threshold:             " << params.scheduling_thresh << std::endl;
  }
  return str;
}

struct RefinementParameters {
  LabelPropagationParameters label_propagation {};
  FlowParameters flow {};
};

inline std::ostream& operator<< (std::ostream& str, const RefinementParameters& params) {
  str << "Refinement Parameters:" << std::endl;
  if ( params.label_propagation.algorithm != LabelPropagationAlgorithm::do_nothing ) {
    str << std::endl << params.label_propagation;
  }
  if ( params.flow.algorithm != FlowAlgorithm::do_nothing ) {
    str << std::endl << params.flow;
  }
  return str;
}

struct SparsificationParameters {
  bool use_degree_zero_contractions = false;
  bool use_heavy_net_removal = false;
  bool use_similiar_net_removal = false;
  double hyperedge_pin_weight_fraction = 0.0;
  size_t min_hash_footprint_size = 0;
  double jaccard_threshold = 1.0;
  SimiliarNetCombinerStrategy similiar_net_combiner_strategy = SimiliarNetCombinerStrategy::UNDEFINED;
  // Those will be determined dynamically
  HypernodeWeight max_hyperedge_pin_weight = std::numeric_limits<HypernodeWeight>::max();
};

inline std::ostream & operator<< (std::ostream& str, const SparsificationParameters& params) {
  str << "Sparsification Parameters:" << std::endl;
  str << "  use degree-zero HN contractions:    " << std::boolalpha << params.use_degree_zero_contractions << std::endl;
  str << "  use heavy net removal:              " << std::boolalpha << params.use_heavy_net_removal << std::endl;
  str << "  use similiar net removal:           " << std::boolalpha << params.use_similiar_net_removal << std::endl;
  if ( params.use_heavy_net_removal ) {
    str << "  hyperedge pin weight fraction:      " << params.hyperedge_pin_weight_fraction << std::endl;
    str << "  maximum hyperedge pin weight:       " << params.max_hyperedge_pin_weight << std::endl;
  }
  if ( params.use_similiar_net_removal ) {
    str << "  min-hash footprint size:            " << params.min_hash_footprint_size << std::endl;
    str << "  jaccard threshold:                  " << params.jaccard_threshold << std::endl;
    str << "  similiar net combiner strategy:     " << params.similiar_net_combiner_strategy << std::endl;
  }
  return str;
}

struct InitialPartitioningParameters {
  InitialPartitioningMode mode = InitialPartitioningMode::UNDEFINED;
  RefinementParameters refinement = { };
  size_t runs = 1;
  bool use_adaptive_epsilon = false;
  size_t lp_maximum_iterations = 1;
  size_t lp_initial_block_size = 1;
};

inline std::ostream & operator<< (std::ostream& str, const InitialPartitioningParameters& params) {
  str << "Initial Partitioning Parameters:" << std::endl;
  str << "  Initial Partitioning Mode:          " << params.mode << std::endl;
  str << "  Number of Runs:                     " << params.runs << std::endl;
  str << "  Use Adaptive Epsilon:               " << std::boolalpha << params.use_adaptive_epsilon << std::endl;
  str << "  Maximum Iterations of LP IP:        " << params.lp_maximum_iterations << std::endl;
  str << "  Initial Block Size of LP IP:        " << params.lp_initial_block_size << std::endl;
  str << "\nInitial Partitioning ";
  str << params.refinement << std::endl;
  return str;
}

struct SharedMemoryParameters {
  size_t num_threads = 1;
  bool use_localized_random_shuffle = false;
  size_t shuffle_block_size = 2;
};

inline std::ostream & operator<< (std::ostream& str, const SharedMemoryParameters& params) {
  str << "Shared Memory Parameters:             " << std::endl;
  str << "  Number of Threads:                  " << params.num_threads << std::endl;
  str << "  Number of used NUMA nodes:          " << TBBNumaArena::instance().num_used_numa_nodes() << std::endl;
  str << "  Use Localized Random Shuffle:       " << std::boolalpha << params.use_localized_random_shuffle << std::endl;
  str << "  Random Shuffle Block Size:          " << params.shuffle_block_size << std::endl;
  return str;
}

class Context {
 public:
  PartitioningParameters partition { };
  PreprocessingParameters preprocessing { };
  CoarseningParameters coarsening { };
  InitialPartitioningParameters initial_partitioning { };
  RefinementParameters refinement { };
  SparsificationParameters sparsification { };
  SharedMemoryParameters shared_memory { };
  kahypar::ContextType type = kahypar::ContextType::main;

  Context() { }

  bool useSparsification() const {
    return sparsification.use_degree_zero_contractions ||
           sparsification.use_heavy_net_removal ||
           sparsification.use_similiar_net_removal;
  }

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

    setupSparsificationParameters();
  }

  void setupContractionLimit(const HypernodeWeight total_hypergraph_weight) {
    // Setup contraction limit
    if (initial_partitioning.mode == InitialPartitioningMode::recursive) {
      coarsening.contraction_limit =
        2 * std::max(shared_memory.num_threads, static_cast<size_t>(partition.k)) *
        coarsening.contraction_limit_multiplier;
    } else {
      coarsening.contraction_limit =
        coarsening.contraction_limit_multiplier * partition.k;
    }

    // Setup maximum allowed vertex and high-degree vertex weight
    setupMaximumAllowedNodeWeight(total_hypergraph_weight);
  }

  void setupMaximumAllowedNodeWeight(const HypernodeWeight total_hypergraph_weight) {
    HypernodeWeight min_block_weight = std::numeric_limits<HypernodeWeight>::max();
    for ( PartitionID part_id = 0; part_id < partition.k; ++part_id ) {
      min_block_weight = std::min(min_block_weight, partition.max_part_weights[part_id]);
    }

    double hypernode_weight_fraction =
      coarsening.max_allowed_weight_multiplier
      / coarsening.contraction_limit;
    coarsening.max_allowed_node_weight =
      std::ceil(hypernode_weight_fraction * total_hypergraph_weight);
    coarsening.max_allowed_node_weight =
      std::min(coarsening.max_allowed_node_weight, min_block_weight);
  }

  void setupSparsificationParameters() {
    if ( sparsification.use_heavy_net_removal ) {
      HypernodeWeight max_block_weight = 0;
      for ( PartitionID block = 0; block < partition.k; ++block ) {
        max_block_weight = std::max(max_block_weight, partition.max_part_weights[block]);
      }

      sparsification.max_hyperedge_pin_weight = max_block_weight /
        sparsification.hyperedge_pin_weight_fraction;
    }
  }

  void sanityCheck() {
    if (partition.objective == kahypar::Objective::cut &&
        refinement.label_propagation.algorithm == LabelPropagationAlgorithm::label_propagation_km1) {
      ALGO_SWITCH("Refinement algorithm" << refinement.label_propagation.algorithm << "only works for km1 metric."
                                         << "Do you want to use the cut version of the label propagation refiner (Y/N)?",
                  "Partitioning with" << refinement.label_propagation.algorithm
                                         << "refiner in combination with cut metric is not possible!",
                  refinement.label_propagation.algorithm,
                  LabelPropagationAlgorithm::label_propagation_cut);
    } else if (partition.objective == kahypar::Objective::km1 &&
               refinement.label_propagation.algorithm == LabelPropagationAlgorithm::label_propagation_cut) {
      ALGO_SWITCH("Refinement algorithm" << refinement.label_propagation.algorithm << "only works for cut metric."
                                         << "Do you want to use the km1 version of the label propagation refiner (Y/N)?",
                  "Partitioning with" << refinement.label_propagation.algorithm
                                         << "refiner in combination with km1 metric is not possible!",
                  refinement.label_propagation.algorithm,
                  LabelPropagationAlgorithm::label_propagation_km1);
    }

    if (partition.objective == kahypar::Objective::cut &&
        initial_partitioning.refinement.label_propagation.algorithm ==
        LabelPropagationAlgorithm::label_propagation_km1) {
      ALGO_SWITCH("Initial Partitioning Refinement algorithm"
                    << initial_partitioning.refinement.label_propagation.algorithm
                    << "only works for km1 metric."
                    << "Do you want to use the cut version of the label propagation refiner (Y/N)?",
                  "Partitioning with" << initial_partitioning.refinement.label_propagation.algorithm
                    << "refiner in combination with cut metric is not possible!",
                  initial_partitioning.refinement.label_propagation.algorithm,
                  LabelPropagationAlgorithm::label_propagation_cut);
    } else if (partition.objective == kahypar::Objective::km1 &&
               initial_partitioning.refinement.label_propagation.algorithm ==
               LabelPropagationAlgorithm::label_propagation_cut) {
      ALGO_SWITCH("Initial Partitioning Refinement algorithm"
                                         << initial_partitioning.refinement.label_propagation.algorithm
                                         << "only works for cut metric."
                                         << "Do you want to use the km1 version of the label propagation refiner (Y/N)?",
                  "Partitioning with" << initial_partitioning.refinement.label_propagation.algorithm
                                         << "refiner in combination with km1 metric is not possible!",
                  initial_partitioning.refinement.label_propagation.algorithm,
                  LabelPropagationAlgorithm::label_propagation_km1);
    }

    if ( !preprocessing.use_community_detection ) {
      if ( preprocessing.use_community_redistribution ) {
        ALGO_SWITCH("Community redistribution only works if community detection is enabled."
                    << "Do you want to enable community detection (Y/N)?",
                    "Community redistribution without community detection is not possible!",
                    preprocessing.use_community_detection,
                    true);
      }
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
      << context.sparsification
      << "-------------------------------------------------------------------------------\n"
      << context.shared_memory
      << "-------------------------------------------------------------------------------";
  return str;
}
}  // namespace mt_kahypar
