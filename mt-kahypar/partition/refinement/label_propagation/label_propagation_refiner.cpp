/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
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

#include "label_propagation_refiner.h"

namespace mt_kahypar {

  template <template <typename> class GainPolicy>
  bool LabelPropagationRefiner<GainPolicy>::refineImpl(PartitionedHypergraph& hypergraph,
                  kahypar::Metrics& best_metrics,
                  const double)  {
          _gain.reset();
          _next_active.reset();

          // Perform Label Propagation
          labelPropagation(hypergraph);

          // Update global part weight and sizes
          best_metrics.imbalance = metrics::imbalance(hypergraph, _context);

          // Update metrics statistics
          HyperedgeWeight current_metric = best_metrics.getMetric(
          kahypar::Mode::direct_kway, _context.partition.objective);
          Gain delta = _gain.delta();
          ASSERT(delta <= 0, "LP refiner worsen solution quality");

          HEAVY_REFINEMENT_ASSERT(current_metric + delta ==
                                  metrics::objective(hypergraph, _context.partition.objective,
                                                     !_context.refinement.label_propagation.execute_sequential),
                                  V(current_metric) << V(delta) <<
                                                    V(metrics::objective(hypergraph, _context.partition.objective,
                                                                         _context.refinement.label_propagation.execute_sequential)));

          best_metrics.updateMetric(current_metric + delta, kahypar::Mode::direct_kway, _context.partition.objective);
          utils::Stats::instance().update_stat("lp_improvement", std::abs(delta));
          return delta < 0;
  }


  template <template <typename> class GainPolicy>
  void LabelPropagationRefiner<GainPolicy>::labelPropagation(PartitionedHypergraph& hypergraph) {
    NextActiveNodes next_active_nodes;
    for (size_t i = 0; i < _context.refinement.label_propagation.maximum_iterations; ++i) {
      DBG << "Starting Label Propagation Round" << i;

      utils::Timer::instance().start_timer(
              "lp_round_" + std::to_string(i), "Label Propagation Round " + std::to_string(i), true);

      if ( _active_nodes.size() > 0 ) {
        labelPropagationRound(hypergraph, next_active_nodes);
      }

      if ( _context.refinement.label_propagation.execute_sequential ) {
        _active_nodes = next_active_nodes.copy_sequential();
        next_active_nodes.clear_sequential();
      } else {
        _active_nodes = next_active_nodes.copy_parallel();
        next_active_nodes.clear_parallel();
      }
      utils::Timer::instance().stop_timer("lp_round_" + std::to_string(i));

      if ( _active_nodes.size() == 0 ) {
        break;
      }
    }
  }

  template <template <typename> class GainPolicy>
  bool LabelPropagationRefiner<GainPolicy>::labelPropagationRound(PartitionedHypergraph& hypergraph,
                             NextActiveNodes& next_active_nodes) {
    _visited_he.reset();
    _next_active.reset();
    // This function is passed as lambda to the changeNodePart function and used
    // to calculate the "real" delta of a move (in terms of the used objective function).
    auto objective_delta = [&](const HyperedgeID he,
                               const HyperedgeWeight edge_weight,
                               const HypernodeID edge_size,
                               const HypernodeID pin_count_in_from_part_after,
                               const HypernodeID pin_count_in_to_part_after) {
      _gain.computeDeltaForHyperedge(he, edge_weight, edge_size,
                                     pin_count_in_from_part_after, pin_count_in_to_part_after);
    };

    // Shuffle Vector
    bool converged = true;
    if ( _context.refinement.label_propagation.execute_sequential ) {
      utils::Randomize::instance().shuffleVector(
              _active_nodes, 0UL, _active_nodes.size(), sched_getcpu());

      for ( size_t j = 0; j < _active_nodes.size(); ++j ) {
        const HypernodeID hn = _active_nodes[j];
        converged &= !moveVertex(hypergraph, hn,
                                 next_active_nodes, objective_delta);
      }
    } else {
      utils::Randomize::instance().parallelShuffleVector(
              _active_nodes, 0UL, _active_nodes.size());

      tbb::parallel_for(0UL, _active_nodes.size(), [&](const size_t& j) {
        const HypernodeID hn = _active_nodes[j];
        converged &= !moveVertex(hypergraph, hn,
                                 next_active_nodes, objective_delta);
      });
    }
    return converged;
  }

  template <template <typename> class GainPolicy>
  void LabelPropagationRefiner<GainPolicy>::initializeImpl(PartitionedHypergraph& hypergraph) {
    ActiveNodes tmp_active_nodes;
    _active_nodes = std::move(tmp_active_nodes);

    if ( _context.refinement.label_propagation.execute_sequential ) {
      // Setup active nodes sequential
      for ( const HypernodeID hn : hypergraph.nodes() ) {
        if ( _context.refinement.label_propagation.rebalancing || hypergraph.isBorderNode(hn) ) {
          _active_nodes.push_back(hn);
        }
      }
    } else {
      // Setup active nodes in parallel
      // A node is active, if it is a border vertex.
      NextActiveNodes tmp_active_nodes;

      hypergraph.doParallelForAllNodes([&](const HypernodeID& hn) {
        if ( _context.refinement.label_propagation.rebalancing || hypergraph.isBorderNode(hn) ) {
          tmp_active_nodes.stream(hn);
        }
      });

      _active_nodes = tmp_active_nodes.copy_parallel();
    }
  }


  // explicitly instantiate so the compiler can generate them when compiling this cpp file
  template class LabelPropagationRefiner<Km1Policy>;
  template class LabelPropagationRefiner<CutPolicy>;
}