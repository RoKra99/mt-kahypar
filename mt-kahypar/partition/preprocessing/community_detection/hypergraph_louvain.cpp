#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_louvain.h"

namespace mt_kahypar::community_detection {
parallel::scalable_vector<HypernodeID> hypergraph_local_moving_contract_recurse(ds::CommunityHypergraph& chg, HypergraphLocalMovingModularity& hlmm) {
    parallel::scalable_vector<HypernodeID> communities(chg.initialNumNodes());
    bool clustering_changed = hlmm.localMoving(chg, communities);
    if (clustering_changed) {
        ds::StaticHypergraph coarse_hg;
        ds::CommunityHypergraph coarse_chg = chg.contract(coarse_hg, communities);

        parallel::scalable_vector<HypernodeID> coarse_communities = hypergraph_local_moving_contract_recurse(chg, hlmm);
        for (size_t i = 0; i < chg.initialNumNodes(); ++i) {
            communities[i] = coarse_communities[communities[i]];
        }
    }
    return communities;
}

parallel::scalable_vector<HypernodeID> hypergraph_louvain(ds::CommunityHypergraph& chg, const bool deactivate_random) {
    HypergraphLocalMovingModularity hlmm(chg, deactivate_random);
    parallel::scalable_vector<HypernodeID> communities = community_detection::hypergraph_local_moving_contract_recurse(chg, hlmm);
    return communities;
}
}