#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_louvain.h"

namespace mt_kahypar::community_detection {

static size_t local_moving_round = 0;
parallel::scalable_vector<HypernodeID> hypergraph_local_moving_contract_recurse(ds::CommunityHypergraph& chg, HypergraphLocalMovingModularity& hlmm) {
    static constexpr bool debug = false;
    parallel::scalable_vector<HypernodeID> communities(chg.initialNumNodes());
    std::cout << chg.initialNumNodes() << ',';
    utils::Timer::instance().start_timer("local_moving " + std::to_string(local_moving_round), "Local Moving" + std::to_string(local_moving_round));
    //utils::Timer::instance().start_timer("hyp_local_moving", "Hypergraph Local Moving");
    bool clustering_changed = hlmm.localMoving(chg, communities);
    //utils::Timer::instance().stop_timer("hyp_local_moving");
    utils::Timer::instance().stop_timer("local_moving " + std::to_string(local_moving_round));
    local_moving_round++;
    if (clustering_changed) {
        ds::StaticHypergraph coarse_hg;
        if (debug) {
            Volume after = metrics::hyp_modularity(chg, communities, hlmm);
            LOG << "after" << after;
        }
        //utils::Timer::instance().start_timer("community_contraction", "Community Hypergraph Contraction");
        ds::CommunityHypergraph coarse_chg = chg.contract(coarse_hg, communities);
        //utils::Timer::instance().stop_timer("community_contraction");

        if (debug) {
            hlmm.initializeCommunityVolumes(chg, communities);
            parallel::scalable_vector<HypernodeID> comm(coarse_chg.initialNumNodes());
            for (HypernodeID i = 0; i < coarse_chg.initialNumNodes(); ++i) {
                comm[i] = i;
            }
            Volume contraction = metrics::hyp_modularity(coarse_chg, comm, hlmm);
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
    HypergraphLocalMovingModularity hlmm(chg, context, deactivate_random);
    parallel::scalable_vector<HypernodeID> communities = community_detection::hypergraph_local_moving_contract_recurse(chg, hlmm);
    // LOG << "Edge Contribution: " << hlmm.edge_contribution_time;
    // LOG << "Expected Edge Contribution: " << hlmm.exp_edge_contribution_time;
    std::cout << hlmm.overall_checks << ',' << hlmm.pruned_by_old << ',' << hlmm.edge_contribution_time << ',' << hlmm.exp_edge_contribution_time << ',';
    return communities;
}
}