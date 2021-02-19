#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_louvain.h"

namespace mt_kahypar::community_detection {

static size_t local_moving_round = 0;
parallel::scalable_vector<HypernodeID> hypergraph_local_moving_contract_recurse(ds::CommunityHypergraph& chg, HypergraphLocalMovingModularity& hlmm) {
    static constexpr bool debug = false;
    parallel::scalable_vector<HypernodeID> communities(chg.initialNumNodes());
    //LOG << local_moving_round << ": " << chg.initialNumNodes();
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
    //LOG << "Edge Contribution: " << hlmm.edge_contribution_time;
    //LOG << "Expected Edge Contribution: " << hlmm.exp_edge_contribution_time;
    //std::cout << hlmm.overall_checks << ',' << hlmm.pruned_by_old << ',' << hlmm.edge_contribution_time << ',' << hlmm.exp_edge_contribution_time << ',';
    auto sep = ',';
    std::cout /*<< "Max" */ << *std::max_element(hlmm.distance.begin(), hlmm.distance.end()) << sep << *std::max_element(hlmm.com_neighbours.begin(), hlmm.com_neighbours.end()) << sep;
    std::cout /*<< "Min"*/ << *std::min_element(hlmm.distance.begin(), hlmm.distance.end()) << sep << *std::min_element(hlmm.com_neighbours.begin(), hlmm.com_neighbours.end()) << sep;
    std::cout /*<< "Avg"*/ << static_cast<Volume>(std::accumulate(hlmm.distance.begin(), hlmm.distance.end(), 0, std::plus<int>())) / hlmm.distance.size()
  << sep << static_cast<Volume>(std::accumulate(hlmm.com_neighbours.begin(), hlmm.com_neighbours.end(), 0, std::plus<int>())) / hlmm.com_neighbours.size() << sep;
    const int middle = hlmm.distance.size() / 2;
    std::sort(hlmm.distance.begin(), hlmm.distance.end());
    std::sort(hlmm.com_neighbours.begin(), hlmm.com_neighbours.end());
    std::cout /*<< "Mean" */ << hlmm.distance[middle] << sep << hlmm.com_neighbours[middle] << sep;
    size_t count = 0;
    while (hlmm.distance[count] == 0) {
        ++count;
    }
    std::cout /*<< "Zeros"*/ << count  << sep << hlmm.distance.size() << sep;
    std::cout /*<< "positive edge count"*/ << hlmm.positive_edgeContribution_count << sep;
    std::cout << hlmm.edge_contribution_time << ',' << hlmm.exp_edge_contribution_time << ',';
    return communities;
}
}