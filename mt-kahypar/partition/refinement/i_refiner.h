/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2014 Sebastian Schlag <sebastian.schlag@kit.edu>
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

#include <array>
#include <string>
#include <utility>
#include <vector>

#include "kahypar/partition/metrics.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
template<typename TypeTraits,
         bool track_border_vertices = TRACK_BORDER_VERTICES>
class IRefinerT {
  using HyperGraph = typename TypeTraits::template PartitionedHyperGraph<track_border_vertices>;

 public:
  IRefinerT(const IRefinerT&) = delete;
  IRefinerT(IRefinerT&&) = delete;
  IRefinerT & operator= (const IRefinerT &) = delete;
  IRefinerT & operator= (IRefinerT &&) = delete;

  virtual ~IRefinerT() = default;

  void initialize(HyperGraph& hypergraph) {
    initializeImpl(hypergraph);
  }

  bool refine(HyperGraph& hypergraph,
              const parallel::scalable_vector<HypernodeID>& refinement_nodes,
              kahypar::Metrics& best_metrics) {
    return refineImpl(hypergraph, refinement_nodes, best_metrics);
  }

 protected:
  IRefinerT() = default;

 private:
  virtual void initializeImpl(HyperGraph& hypergraph) = 0;

  virtual bool refineImpl(HyperGraph& hypergraph,
                          const parallel::scalable_vector<HypernodeID>& refinement_nodes,
                          kahypar::Metrics& best_metrics) = 0;
};

using IRefiner = IRefinerT<GlobalTypeTraits>;

}  // namespace mt_kahypar
