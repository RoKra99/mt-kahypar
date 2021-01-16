#include "community_hypergraph.h"

#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/datastructures/concurrent_bucket_map.h"
#include "mt-kahypar/utils/timer.h"

#include <tbb/parallel_reduce.h>

namespace mt_kahypar::ds {

//TODO: parallel contract community count
//TODO: remove "single pin edges" / save them for dbg purposes
//TODO: bug affecting edge contribution calculation (slight difference)
CommunityHypergraph CommunityHypergraph::contract(StaticHypergraph& hypergraph, parallel::scalable_vector<HypernodeID>& communities) {
    hypergraph = _hg->contract(communities, 0, false);
    CommunityHypergraph chg;
    chg._hg = &hypergraph;
    chg._vol_v = _vol_v;
    chg._valid_edge_sizes = std::move(_valid_edge_sizes);
    chg._d_edge_weights = std::move(_d_edge_weights);
    chg._total_edge_weight = _total_edge_weight;

    // reset buffer values
    Array<parallel::AtomicWrapper<HyperedgeWeight>>& tmp_node_volumes = _tmp_community_hypergraph_buffer->tmp_node_volumes;
    auto& tmp_community_counts = _tmp_community_hypergraph_buffer->tmp_community_counts;

    tbb::parallel_for(0U, initialNumNodes(), [&](const HypernodeID i) {
        tmp_node_volumes[i].store(0);
        });

    // map according to hypergraph mapping
    tbb::parallel_for(0U, initialNumNodes(), [&](const HypernodeID i) {
        const HypernodeID coarse_i = communities[i];
        tmp_node_volumes[coarse_i] += nodeVolume(i);
        });

    chg._node_volumes.resize(chg.initialNumNodes());

    // actually writing the volumes to the graph
    tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID i) {
        chg._node_volumes[i] = std::move(tmp_node_volumes[i]);
        });

    for (HyperedgeID i = 0; i < initialNumEdges(); ++i) {
        const HyperedgeID coarse_i = hypergraph._tmp_contraction_buffer->valid_hyperedges[i];
        if (coarse_i && !tmp_community_counts[coarse_i - 1]) {
            tmp_community_counts.insert(tmp_community_counts.begin() + coarse_i - 1, std::move(_community_counts[i]));
        }
    }

    tbb::parallel_for(0U, chg.initialNumEdges(), [&](const HyperedgeID i) {
            tmp_community_counts[i]->contract(communities);
        });

    chg._community_counts = std::move(tmp_community_counts);


    chg._tmp_community_hypergraph_buffer = _tmp_community_hypergraph_buffer;
    _tmp_community_hypergraph_buffer = nullptr;

    return chg;
}
}