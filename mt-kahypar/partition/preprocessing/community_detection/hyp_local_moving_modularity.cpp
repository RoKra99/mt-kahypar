#include "hyp_local_moving_modularity.h"

namespace mt_kahypar::metrics {
double hyp_modularity(const ds::CommunityHypergraph& hypergraph, ds::Clustering& communities) {
    const HyperedgeWeight vol_V = hypergraph.totalVolume();
    double edge_contribution = 0.0;
    double exp_edge_contribution = 0.0;
    // zero indicates, that this community is not a neighbour
    std::vector<HyperedgeWeight> weight_to_community(communities.size(), 0);
    std::vector<PartitionID> neigh_communities;
    for (const HyperedgeID& he : hypergraph.edges()) {
        for (const HypernodeID& n : hypergraph.pins(he)) {
            const PartitionID community_n = communities[n];
            if (!weight_to_community[community_n]) { // first time finding this community
                weight_to_community[community_n] = hypergraph.edgeWeight(he);
                neigh_communities.emplace_back(community_n);
                exp_edge_contribution += hypergraph.edgeWeight(he) * (1 - pow(1 - hypergraph.communityVolume(community_n) / hypergraph.totalVolume(), hypergraph.edgeSize(he)));
            }
        }
        for (const PartitionID& p : neigh_communities) {
            edge_contribution += weight_to_community[p];
            weight_to_community[p] = 0;
        }
        neigh_communities.clear();
    }
    return  (edge_contribution - exp_edge_contribution) / hypergraph.totalEdgeWeight();
}
}