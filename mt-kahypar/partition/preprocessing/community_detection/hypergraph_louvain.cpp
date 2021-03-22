#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_louvain.h"

namespace mt_kahypar::community_detection {

static size_t local_moving_round = 0;
template<typename HypergraphLocalMoving>
parallel::scalable_vector<HypernodeID> hypergraph_local_moving_contract_recurse(ds::CommunityHypergraph& chg, HypergraphLocalMoving& hlmm) {
    static constexpr bool debug = false;
    parallel::scalable_vector<HypernodeID> communities(chg.initialNumNodes());
    utils::Timer::instance().start_timer("local_moving " + std::to_string(local_moving_round), "Local Moving" + std::to_string(local_moving_round));
    //utils::Timer::instance().start_timer("hyp_local_moving", "Hypergraph Local Moving");
    bool clustering_changed = hlmm.localMoving(chg, communities);
    //utils::Timer::instance().stop_timer("hyp_local_moving");
    utils::Timer::instance().stop_timer("local_moving " + std::to_string(local_moving_round));
    local_moving_round++;
    if (clustering_changed) {
        ds::StaticHypergraph coarse_hg;
        if (debug) {
            //Volume after = metrics::hyp_modularity(chg, communities, hlmm);
            Volume after = metrics::hyp_map_equation(chg, communities);
            LOG << "after" << after;
        }
        ds::CommunityHypergraph coarse_chg = chg.contract(coarse_hg, communities);
        if (debug) {
            // hlmm.initializeCommunityVolumes(chg, communities);
            parallel::scalable_vector<HypernodeID> comm(coarse_chg.initialNumNodes());
            for (HypernodeID i = 0; i < coarse_chg.initialNumNodes(); ++i) {
                comm[i] = i;
            }
            // Volume contraction = metrics::hyp_modularity(coarse_chg, comm, hlmm);
            Volume contraction = metrics::hyp_map_equation(coarse_chg, comm);
            LOG << "contraction" << contraction;
        }
        parallel::scalable_vector<HypernodeID> coarse_communities = hypergraph_local_moving_contract_recurse(coarse_chg, hlmm);
        tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID i) {
            ASSERT(communities[i] < static_cast<HypernodeID>(coarse_communities.size()));
            communities[i] = coarse_communities[communities[i]];
        });

    }
    return communities;
}

parallel::scalable_vector<HypernodeID> hypergraph_louvain(ds::CommunityHypergraph& chg, const Context& context, const bool deactivate_random) {
    //HypergraphLocalMovingModularity hlmm(chg, context, deactivate_random);
    HypergraphLocalMovingMapEquation hlmm(chg, context, deactivate_random);
    parallel::scalable_vector<HypernodeID> communities = community_detection::hypergraph_local_moving_contract_recurse(chg, hlmm);
    // for (size_t i = local_moving_round; i <= 5; ++i) {
    //     std::cout << 0 << ',' << 0 << ',';
    // }
    return communities;
}
}