#include "hyp_local_moving_modularity.h"
#include <iostream>

namespace mt_kahypar::metrics {
double hyp_modularity(const ds::CommunityHypergraph& hypergraph) {
    const HyperedgeWeight vol_V = hypergraph.totalVolume();
    double edge_contribution = 0.0;
    double exp_edge_contribution = 0.0;
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
    for (HypernodeID d = 0; d < hypergraph.maxEdgeSize(); ++d) {
        double d_chance = 0.0;
        for (const HyperedgeWeight& vol_c : hypergraph.communityVolumes()) {
            d_chance += 1 - pow(1.0 - static_cast<double>(vol_c) / vol_V, d);
        }
        exp_edge_contribution += hypergraph.dEdgeWeight(d) * d_chance;
    }
    return  (edge_contribution - exp_edge_contribution) / hypergraph.totalEdgeWeight();
}
}