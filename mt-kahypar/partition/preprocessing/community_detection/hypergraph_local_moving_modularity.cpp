#include "hypergraph_local_moving_modularity.h"
#include <iostream>

namespace mt_kahypar::metrics {
Volume hyp_modularity(const ds::CommunityHypergraph& hypergraph) {
    static constexpr bool debug = false;
    const HyperedgeWeight vol_V = hypergraph.totalVolume();
    Volume edge_contribution = 0.0L;
    // zero indicates, that this community is not a neighbour
    std::vector<HyperedgeWeight> weight_to_community(hypergraph.initialNumNodes(), 0);
    std::vector<PartitionID> neigh_communities;
    for (const HyperedgeID& he : hypergraph.edges()) {
        for (const HypernodeID& n : hypergraph.pins(he)) {
            const PartitionID community_n = hypergraph.communityID(n);
            if (!weight_to_community[community_n]) { // first time finding this community
                weight_to_community[community_n] = hypergraph.edgeWeight(he);
                neigh_communities.emplace_back(community_n);
            }
        }
        for (const PartitionID& p : neigh_communities) {
            edge_contribution += weight_to_community[p];
            weight_to_community[p] = 0;
        }
        neigh_communities.clear();
    }
    Volume exp_edge_contribution = 0.0;
    for (const size_t d : hypergraph.edgeSizes()) {
        Volume d_chance = 0.0L;
        for (const HyperedgeWeight& vol_c : hypergraph.communityVolumes()) {
            DBG << "Volumes: " << vol_c;
            d_chance += 1.0L - static_cast<Volume>(math::fast_power(vol_V - vol_c , d)) / math::fast_power(vol_V,d);
        }
        exp_edge_contribution += hypergraph.edgeWeightBySize(d) * d_chance;
        DBG << "Edgesize: " << d << ", d_chance " << d_chance << " what is added: " << hypergraph.edgeWeightBySize(d) * d_chance;
    }
    return  (edge_contribution - exp_edge_contribution) / hypergraph.totalEdgeWeight();
}
}