#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_map_equation.h" 

namespace mt_kahypar::metrics {
double hyp_map_equation(
    const ds::CommunityHypergraph& chg,
    const parallel::scalable_vector<HypernodeID>& communities) {
    ASSERT(chg.initialNumNodes() == communities.size());
    const double reciprocal_vol_total = 1.0 / chg.totalVolume();
    parallel::scalable_vector<parallel::AtomicWrapper<HyperedgeWeight>> community_volumes(chg.initialNumNodes());
    parallel::scalable_vector<parallel::AtomicWrapper<uint32_t>> exit_prob_vol_total(chg.initialNumNodes());
    tbb::enumerable_thread_specific<uint32_t> sum_exit_prob_local(0);
    const auto plogp_rel = [reciprocal_vol_total](double p) -> double {
        if (p > 0.0) {
            p *= reciprocal_vol_total;
            return p * log2(p);
        }
        return 0.0;
    };
    tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID hn) {
        community_volumes[communities[hn]] += chg.nodeVolume(hn);
        for (const HyperedgeID he : chg.incidentEdges(hn)) {
            exit_prob_vol_total[communities[hn]] += chg.edgeSize(he) * chg.edgeWeight(he);
        }
        exit_prob_vol_total[hn] -= chg.nodeVolume(hn);
        });
    tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID hn) {
        sum_exit_prob_local.local() += exit_prob_vol_total[hn];
        });
    const uint32_t sum_exit_prob = sum_exit_prob_local.combine(std::plus<>());


    tbb::enumerable_thread_specific<double> sum_plogp_exit_prob_local(0.0);
    tbb::enumerable_thread_specific<double> sum_plogp_exit_prob_plus_com_vol_local(0.0);
    tbb::enumerable_thread_specific<double> sum_plogp_prob_in_node_local(0.0);

    tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID i) {
        sum_plogp_exit_prob_local.local() += plogp_rel(exit_prob_vol_total[i]);
        sum_plogp_exit_prob_plus_com_vol_local.local() += plogp_rel(exit_prob_vol_total[i] + community_volumes[i]);
        // this is not neccessary for optimization
        sum_plogp_prob_in_node_local.local() += plogp_rel(chg.nodeVolume(i));
        });
    const double plogp_sum_exit_prob = plogp_rel(sum_exit_prob);
    const double sum_plogp_exit_prob = sum_plogp_exit_prob_local.combine(std::plus<>());
    const double sum_plogp_exit_prob_plus_com_vol = sum_plogp_exit_prob_plus_com_vol_local.combine(std::plus<>());;
    // this is not neccessary for optimization
    const double sum_plogp_prob_in_node = sum_plogp_prob_in_node_local.combine(std::plus<>());

    return plogp_sum_exit_prob - 2 * sum_plogp_exit_prob + sum_plogp_exit_prob_plus_com_vol - sum_plogp_prob_in_node;
}
}