#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_louvain.h"

namespace mt_kahypar::community_detection {
ds::Clustering hypergraph_local_moving_contract_recurse(ds::CommunityHypergraph& chg, HypergraphLocalMovingModularity& hlmm) {
    ds::Clustering communities(chg.initialNumNodes());
    bool clustering_changed = hlmm.localMoving(communities);

    if (clustering_changed) {

        ds::CommunityHypergraph coars_chg = chg.contract(/*communities*/);

        ds::Clustering coarse_communities = hypergraph_local_moving_contract_recurse(chg, hlmm);
        for (size_t i = 0; i < chg.initialNumNodes(); ++i) {
            communities[i] = coarse_communities[communities[i]];
        }
    }
    return communities;
}

ds::Clustering hypergraph_louvain(ds::CommunityHypergraph& chg) {
    HypergraphLocalMovingModularity hlmm(chg);
    ds::Clustering communities = community_detection::hypergraph_local_moving_contract_recurse(chg, hlmm);
    return communities;
}
}