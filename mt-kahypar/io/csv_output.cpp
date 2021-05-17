/*******************************************************************************
 * This file is part of MT-KaHyPar.
 *
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

#include "csv_output.h"

#include <sstream>

#include "mt-kahypar/partition/metrics.h"
#include "mt-kahypar/utils/initial_partitioning_stats.h"
#include "mt-kahypar/utils/stats.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_modularity.h"

namespace mt_kahypar::io::csv {

  std::string header() {
    return "algorithm,threads,graph,k,seed,epsilon,imbalance,"
           "objective,km1,cut,partitionTime,fmTime,lpTime,coarseningTime,ipTime,preprocessingTime,lmRound0,kmRound1,lmRound2,lmRound3,kmRound4"
           "\n";
  }

  std::string serialize(const PartitionedHypergraph& phg, const Context& context,
                                      const std::chrono::duration<double>& elapsed_seconds) {
    const char sep = ',';
    std::stringstream s;

    s << context.algorithm_name;
    if (context.algorithm_name == "MT-KaHyPar") {
      if (context.partition.preset_file.find("fast") != std::string::npos) {
        s << "-Fast";
      } else if (context.partition.preset_file.find("quality") != std::string::npos) {
        s << "-Eco";
      }
    }
    s << sep;

    s << context.shared_memory.num_threads << sep;
    s << context.partition.graph_filename.substr(context.partition.graph_filename.find_last_of('/') + 1) << sep;
    s << context.partition.k << sep;
    s << context.partition.seed << sep;

    s << context.partition.epsilon << sep;
    s << metrics::imbalance(phg, context) << sep;

    s << context.partition.objective << sep;
    s << metrics::km1(phg) << sep;
    s << metrics::hyperedgeCut(phg) << sep;
    s << elapsed_seconds.count() << sep;

    utils::Timer& timer = utils::Timer::instance(context.partition.show_detailed_timings);
    s << (timer.get("fm") + timer.get("initialize_fm_refiner"))<< sep;
    s << (timer.get("label_propagation") + timer.get("initialize_lp_refiner")) << sep;
    s << timer.get("coarsening") << sep;
    s << timer.get("initial_partitioning") << sep;
    s << timer.get("preprocessing") << sep;
    s << timer.get("local_moving 0") << sep;
    s << timer.get("local_moving 1") << sep;
    s << timer.get("local_moving 2") << sep;
    s << timer.get("local_moving 3") << sep;
    s << timer.get("local_moving 4");
    //s << timer.get("edge_contribution") << sep;
    //s << timer.get("exp_edge_contribution");

    return s.str();
  }
}