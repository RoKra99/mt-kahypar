#include "community_hypergraph.h"

#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar::ds {

//TODO:  "single pin edges" save them for dbg purposes
CommunityHypergraph CommunityHypergraph::contract(StaticHypergraph& hypergraph, parallel::scalable_vector<HypernodeID>& communities) {
    //utils::Timer::instance().start_timer("community_hypergraph_contract", "CommunityHypergaph Contraction");
    hypergraph = _hg->contract(communities, 0);
    //utils::Timer::instance().start_timer("community_specific", "Community Specific");
    CommunityHypergraph chg(_context);
    chg._hg = &hypergraph;
    chg._vol_v = _vol_v;
    chg._valid_edge_sizes = std::move(_valid_edge_sizes);
    chg._d_edge_weights = std::move(_d_edge_weights);
    chg._total_edge_weight = _total_edge_weight;
    chg._community_count_locks = std::move(_community_count_locks);

    Array<parallel::AtomicWrapper<HyperedgeWeight>>& tmp_node_volumes = _tmp_community_hypergraph_buffer->tmp_node_volumes;

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

    chg._community_counts.resize(chg.initialNumEdges());
    tbb::parallel_for(0U, chg.initialNumEdges(), [&](const HyperedgeID he) {
        chg._community_counts[he] = chg.edgeSize(he) > _context.preprocessing.community_detection.hyperedge_size_caching_threshold
            ? std::make_unique<CommunityCount<Map>>(chg.edgeSize(he), chg.pins(he), true) : std::unique_ptr<CommunityCount<Map>>(nullptr);
        });

    chg._tmp_community_hypergraph_buffer = _tmp_community_hypergraph_buffer;
    _tmp_community_hypergraph_buffer = nullptr;
    //utils::Timer::instance().stop_timer("community_specific");
    //utils::Timer::instance().stop_timer("community_hypergraph_contract");
    return chg;
}
}