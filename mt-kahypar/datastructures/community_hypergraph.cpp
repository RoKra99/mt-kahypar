#include "community_hypergraph.h"

#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/datastructures/concurrent_bucket_map.h"
#include "mt-kahypar/utils/timer.h"

#include <tbb/parallel_reduce.h>

namespace mt_kahypar::ds {

StaticHypergraph CommunityHypergraph::contract(parallel::scalable_vector<HypernodeID>& communities) {
    return _hg->contract(communities, 0, false);
}

CommunityHypergraph CommunityHypergraph::mapContractedVolumes(StaticHypergraph& hypergraph) {
    CommunityHypergraph chg;
    chg._hg = &hypergraph;
    chg._vol_v = _vol_v;
    chg._valid_edge_sizes = std::move(_valid_edge_sizes);
    chg._d_edge_weights = std::move(_d_edge_weights);
    chg._total_edge_weight = _total_edge_weight;
    // map according to prefixsum
    // chg._community_volumes
    // chg._node_volumes
    //############################## Keine Ahnung wie ########################
    // chg.community_counts
    return chg;
}

}