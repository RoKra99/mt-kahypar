#include "hypergraph_local_moving_modularity.h"
#include <iostream>

namespace mt_kahypar::metrics {
Volume hyp_modularity(const ds::CommunityHypergraph& hypergraph, const parallel::scalable_vector<HypernodeID>& communities, const community_detection::HypergraphLocalMovingModularity& hlmm) {
    static constexpr bool debug = false;
    const HyperedgeWeight vol_V = hypergraph.totalVolume();
    Volume edge_contribution = 0.0L;
    // zero indicates, that this community is not a neighbour
    std::vector<HyperedgeWeight> weight_to_community(hypergraph.initialNumNodes(), 0);
    std::vector<PartitionID> neigh_communities;
    for (const HyperedgeID& he : hypergraph.edges()) {
        for (const HypernodeID& n : hypergraph.pins(he)) {
            const PartitionID community_n = communities[n];
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
    DBG << "edge contribution" << edge_contribution;
    Volume exp_edge_contribution = 0.0;
    for (const auto& d_pair : hypergraph.edgeSizes()) {
        Volume d_chance = 0.0L;
        const size_t d = d_pair.index;
        for (const HyperedgeWeight& vol_c : hlmm.communityVolumes(hypergraph)) {
            d_chance += 1.0L - static_cast<Volume>(math::fast_power(vol_V - vol_c, d)) / math::fast_power(vol_V, d);
        }
        exp_edge_contribution += d_pair.weight * d_chance;
    }
    DBG << "exp edge contribution" << exp_edge_contribution;
    return  (edge_contribution - exp_edge_contribution);
}
}