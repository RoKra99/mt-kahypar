#include "hypergraph_local_moving_modularity.h"
#include <iostream>

namespace mt_kahypar::metrics {
Volume hyp_modularity(const ds::CommunityHypergraph& chg, const parallel::scalable_vector<HypernodeID>& communities) {
    static constexpr bool debug = false;
    const HyperedgeWeight vol_total = chg.totalVolume();
    const size_t number_of_communities = static_cast<size_t>(*std::max_element(communities.begin(), communities.end())) + 1;
    parallel::scalable_vector<parallel::AtomicWrapper<HyperedgeWeight>> community_volumes(number_of_communities);
    tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID hn) {
        community_volumes[communities[hn]] += chg.nodeVolume(hn);
    });

    Volume edge_contribution = 0.0L;
    // zero indicates, that this community is not a neighbour
    std::vector<HyperedgeWeight> weight_to_community(chg.initialNumNodes(), 0);
    std::vector<PartitionID> neigh_communities;
    for (const HyperedgeID& he : chg.edges()) {
        for (const HypernodeID& n : chg.pins(he)) {
            const PartitionID community_n = communities[n];
            if (!weight_to_community[community_n]) { // first time finding this community
                weight_to_community[community_n] = chg.edgeWeight(he);
                neigh_communities.emplace_back(community_n);
            }
        }
        for (const PartitionID& p : neigh_communities) {
            edge_contribution += weight_to_community[p];
            weight_to_community[p] = 0;
        }
        neigh_communities.clear();
    }
    DBG << "edge contribution" << edge_contribution;
    Volume exp_edge_contribution = 0.0;
    for (const size_t d : chg.edgeSizes()) {
        Volume d_chance = 0.0L;
        for (const HyperedgeWeight vol_c : community_volumes) {
            if (vol_c != vol_total) {
                d_chance += 1.0L - math::fast_power(vol_total - vol_c, d) / math::fast_power(vol_total, d);
            }
        }
        exp_edge_contribution += chg.edgeWeightBySize(d) * d_chance;
    }
    DBG << "exp edge contribution" << exp_edge_contribution;
    return  (edge_contribution - exp_edge_contribution);
}
}