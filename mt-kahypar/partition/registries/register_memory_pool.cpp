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


#include "register_memory_pool.h"

#include "mt-kahypar/datastructures/pin_count_in_part.h"
#include "mt-kahypar/datastructures/connectivity_set.h"
#include "mt-kahypar/parallel/memory_pool.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/memory_tree.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar {

  void register_memory_pool(const Hypergraph& hypergraph,
                            const Context& context) {

    // ########## Preprocessing Memory ##########

    const HypernodeID num_hypernodes = hypergraph.initialNumNodes();
    const HyperedgeID num_hyperedges = hypergraph.initialNumEdges();
    const HypernodeID num_pins = hypergraph.initialNumPins();
    const HypernodeID max_edge_size = hypergraph.maxEdgeSize();

    auto& pool = parallel::MemoryPool::instance();

    if ( context.preprocessing.use_community_detection ) {
      const bool is_graph = hypergraph.maxEdgeSize() == 2;
      const size_t num_star_expansion_nodes = num_hypernodes + (is_graph ? 0 : num_hyperedges);
      const size_t num_star_expansion_edges = is_graph ? num_pins : (2UL * num_pins);

      pool.register_memory_group("Preprocessing", 1);
      // Hyperedge Modularity
      pool.register_memory_chunk("Preprocessing", "node_volumes", num_hypernodes, sizeof(HyperedgeWeight));
      pool.register_memory_chunk("Preprocessing", "d_edge_weights", max_edge_size + 1, sizeof(HyperedgeWeight));
      pool.register_memory_chunk("Preprocessing", "clearlist_edge_contribution", num_hypernodes, sizeof(HyperedgeWeight));
      pool.register_memory_chunk("Preprocessing", "powers_of_source_community", max_edge_size + 1, sizeof(Volume));
      pool.register_memory_chunk("Preprocessing", "tmp_community_volumes", num_hypernodes, sizeof(parallel::AtomicWrapper<HyperedgeWeight>));

      // Old Modualrity
      // pool.register_memory_chunk("Preprocessing", "indices", num_star_expansion_nodes + 1, sizeof(size_t));
      // pool.register_memory_chunk("Preprocessing", "arcs", num_star_expansion_edges, sizeof(Arc));
      // pool.register_memory_chunk("Preprocessing", "node_volumes", num_star_expansion_nodes, sizeof(ArcWeight));
      // pool.register_memory_chunk("Preprocessing", "tmp_indices",
      //                            num_star_expansion_nodes + 1, sizeof(parallel::IntegralAtomicWrapper<size_t>));
      // pool.register_memory_chunk("Preprocessing", "tmp_pos",
      //                            num_star_expansion_nodes, sizeof(parallel::IntegralAtomicWrapper<size_t>));
      // pool.register_memory_chunk("Preprocessing", "tmp_arcs", num_star_expansion_edges, sizeof(Arc));
      // pool.register_memory_chunk("Preprocessing", "valid_arcs", num_star_expansion_edges, sizeof(size_t));
      // pool.register_memory_chunk("Preprocessing", "tmp_node_volumes",
      //                            num_star_expansion_nodes, sizeof(parallel::AtomicWrapper<ArcWeight>));
    }

    // ########## Coarsening Memory ##########

    pool.register_memory_group("Coarsening", 2);
    pool.register_memory_chunk("Coarsening", "mapping", num_hypernodes, sizeof(size_t));
    pool.register_memory_chunk("Coarsening", "tmp_hypernodes", num_hypernodes, Hypergraph::SIZE_OF_HYPERNODE);
    pool.register_memory_chunk("Coarsening", "tmp_incident_nets", num_pins, sizeof(HyperedgeID));
    pool.register_memory_chunk("Coarsening", "tmp_num_incident_nets",
                               num_hypernodes, sizeof(parallel::IntegralAtomicWrapper<size_t>));
    pool.register_memory_chunk("Coarsening", "hn_weights",
                               num_hypernodes, sizeof(parallel::IntegralAtomicWrapper<HypernodeWeight>));
    pool.register_memory_chunk("Coarsening", "tmp_hyperedges", num_hyperedges, Hypergraph::SIZE_OF_HYPEREDGE);
    pool.register_memory_chunk("Coarsening", "tmp_incidence_array", num_pins, sizeof(HypernodeID));
    pool.register_memory_chunk("Coarsening", "he_sizes", num_hyperedges, sizeof(size_t));
    pool.register_memory_chunk("Coarsening", "valid_hyperedges", num_hyperedges, sizeof(size_t));

    // ########## Refinement Memory ##########

    pool.register_memory_group("Refinement", 3);
    const HypernodeID max_he_size = hypergraph.maxEdgeSize();
    pool.register_memory_chunk("Refinement", "part_ids", num_hypernodes, sizeof(PartitionID));
    pool.register_memory_chunk("Refinement", "pin_count_in_part",
                               ds::PinCountInPart::num_elements(num_hyperedges, context.partition.k, max_he_size),
                               sizeof(ds::PinCountInPart::Value));
    pool.register_memory_chunk("Refinement", "connectivity_set",
                               ds::ConnectivitySets::num_elements(num_hyperedges, context.partition.k),
                               sizeof(ds::ConnectivitySets::UnsafeBlock));
    pool.register_memory_chunk("Refinement", "move_to_penalty",
                               static_cast<size_t>(num_hypernodes) * context.partition.k, sizeof(CAtomic<HyperedgeWeight>));
    pool.register_memory_chunk("Refinement", "move_from_penalty",
                               num_hypernodes, sizeof(CAtomic<HyperedgeWeight>));
    pool.register_memory_chunk("Refinement", "pin_count_update_ownership",
                               num_hyperedges, sizeof(parallel::IntegralAtomicWrapper<bool>));

    if ( context.refinement.fm.algorithm != FMAlgorithm::do_nothing && context.refinement.fm.revert_parallel ) {
      pool.register_memory_chunk("Refinement", "remaining_original_pins",
                                 static_cast<size_t>(hypergraph.numNonGraphEdges()) * context.partition.k, sizeof(CAtomic<HypernodeID>));
      pool.register_memory_chunk("Refinement", "first_move_in",
                                 static_cast<size_t>(hypergraph.numNonGraphEdges()) * context.partition.k, sizeof(CAtomic<MoveID>));
      pool.register_memory_chunk("Refinement", "last_move_out",
                                 static_cast<size_t>(hypergraph.numNonGraphEdges()) * context.partition.k, sizeof(CAtomic<MoveID>));
    }

    // Allocate Memory
    utils::Timer::instance().start_timer("memory_pool_allocation", "Memory Pool Allocation");
    pool.allocate_memory_chunks<TBBNumaArena>();
    utils::Timer::instance().stop_timer("memory_pool_allocation");
  }


} // namespace mt_kahypar