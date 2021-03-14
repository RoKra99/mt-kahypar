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
    parallel::scalable_vector<parallel::AtomicWrapper<size_t>> exit_prob_vol_total(community_count);
    tbb::enumerable_thread_specific<size_t> sum_exit_prob_local(0);
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

    tbb::enumerable_thread_specific<parallel::scalable_vector<size_t>> overlap_local(community_count, 0);
    tbb::enumerable_thread_specific<parallel::scalable_vector<HypernodeID>> neighbouring_communities;

    tbb::parallel_for(0U, chg.initialNumEdges(), [&](const HyperedgeID he) {
        const size_t edge_size = chg.edgeSize(he);
        const size_t edge_weight = chg.edgeWeight(he);

        for (const auto& mp : chg.multipins(he)) {
            const HypernodeID community = communities[mp.id];
            const size_t pincount = mp.multiplicity;
            if (overlap_local.local()[community] == 0) {
                neighbouring_communities.local().push_back(community);
            }
            overlap_local.local()[community] += pincount;
        }
        for (const HypernodeID community : neighbouring_communities.local()) {
            exit_prob_vol_total[community] += (edge_size - overlap_local.local()[community]) * edge_weight;
            overlap_local.local()[community] = 0;
        }
        neighbouring_communities.local().clear();
        });

    tbb::parallel_for(0UL, community_count, [&](const HypernodeID com) {
        sum_exit_prob_local.local() += exit_prob_vol_total[com];
        });
    const size_t sum_exit_prob = sum_exit_prob_local.combine(std::plus<>());

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