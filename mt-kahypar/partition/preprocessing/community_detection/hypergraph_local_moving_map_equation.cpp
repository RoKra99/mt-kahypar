#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_map_equation.h" 

namespace mt_kahypar::metrics {
// currently only correct for graphs with one community per node
// TODO: implement for evry graph
double hyp_map_equation(
    const ds::CommunityHypergraph& chg,
    const parallel::scalable_vector<HypernodeID>& communities) {
    ASSERT(chg.initialNumNodes() == communities.size());
    const double reciprocal_vol_total = 1.0 / chg.totalVolume();
    const size_t community_count = *std::max_element(communities.begin(), communities.end()) + 1;
    parallel::scalable_vector<parallel::AtomicWrapper<HyperedgeWeight>> community_volumes(community_count);
    parallel::scalable_vector<parallel::AtomicWrapper<double>> exit_prob_vol_total(community_count);
    tbb::enumerable_thread_specific<double> sum_exit_prob_local(0.0);
    const auto plogp_rel = [reciprocal_vol_total](double p) -> double {
        if (p > 0.0) {
            p *= reciprocal_vol_total;
            return p * log2(p);
        }
        return 0.0;
    };

    tbb::enumerable_thread_specific<double> sum_plogp_prob_in_node_local(0.0);
    tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID hn) {
        community_volumes[communities[hn]] += chg.nodeVolume(hn);
        // this is not neccessary for optimization
        sum_plogp_prob_in_node_local.local() += plogp_rel(chg.nodeVolume(hn));
        });

    // this is not neccessary for optimization
    const double sum_plogp_prob_in_node = sum_plogp_prob_in_node_local.combine(std::plus<>());

    tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeWeight>> overlap_local(community_count, 0);
    tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeID>> neighbouring_communities;

    tbb::parallel_for(0U, chg.initialNumEdges(), [&](const HyperedgeID he) {
        for (const auto& mp : chg.multipins(he)) {
            const size_t community_mp = communities[mp.id];
            if (overlap_local.local()[community_mp] == 0) {
                neighbouring_communities.local().push_back(community_mp);
            }
            overlap_local.local()[community_mp] += mp.multiplicity;
        }

        for (const HypernodeID comm : neighbouring_communities.local()) {
            const HypernodeWeight pincount_in_edge = overlap_local.local()[comm];
            const HypernodeWeight edge_size = static_cast<HypernodeWeight>(chg.edgeSize(he));
            //TODO: Not sure if / (edge_size - 1) is better (since that results in an equal model to the original map equation)
            exit_prob_vol_total[comm] += static_cast<double>(pincount_in_edge * chg.edgeWeight(he) * (edge_size - pincount_in_edge)) / edge_size;
            overlap_local.local()[comm] = 0;
        }
        neighbouring_communities.local().clear();
        });

    tbb::parallel_for(0UL, community_count, [&](const HypernodeID com) {
        sum_exit_prob_local.local() += exit_prob_vol_total[com];
        });
    const double sum_exit_prob = sum_exit_prob_local.combine(std::plus<>());

    tbb::enumerable_thread_specific<double> sum_plogp_exit_prob_local(0.0);
    tbb::enumerable_thread_specific<double> sum_plogp_exit_prob_plus_com_vol_local(0.0);

    tbb::parallel_for(0UL, community_count, [&](const HypernodeID i) {
        sum_plogp_exit_prob_local.local() += plogp_rel(exit_prob_vol_total[i]);
        sum_plogp_exit_prob_plus_com_vol_local.local() += plogp_rel(exit_prob_vol_total[i] + community_volumes[i]);
        });

    const double plogp_sum_exit_prob = plogp_rel(sum_exit_prob);
    const double sum_plogp_exit_prob = sum_plogp_exit_prob_local.combine(std::plus<>());
    const double sum_plogp_exit_prob_plus_com_vol = sum_plogp_exit_prob_plus_com_vol_local.combine(std::plus<>());;
    // LOG << plogp_sum_exit_prob;
    // LOG << sum_plogp_exit_prob;
    // LOG << sum_plogp_exit_prob_plus_com_vol;
    // LOG << sum_plogp_prob_in_node;

    return plogp_sum_exit_prob - 2 * sum_plogp_exit_prob + sum_plogp_exit_prob_plus_com_vol - sum_plogp_prob_in_node;
}
}