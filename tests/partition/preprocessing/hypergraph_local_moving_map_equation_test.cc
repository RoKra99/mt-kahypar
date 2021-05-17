#include "gmock/gmock.h"

#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_map_equation.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "tests/datastructures/hypergraph_fixtures.h"
#include <mt-kahypar/io/hypergraph_io.h>
#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_louvain.h"

using ::testing::Test;

namespace mt_kahypar {
namespace community_detection {

class AHypergraphLocalMovingMapEquation : public ds::HypergraphFixture<Hypergraph, HypergraphFactory> {
    using Base = ds::HypergraphFixture<Hypergraph, HypergraphFactory>;

public:
    AHypergraphLocalMovingMapEquation() :
        Base(),
        chg(nullptr),
        context() {
        context.preprocessing.community_detection.max_pass_iterations = 5;
        context.preprocessing.community_detection.min_vertex_move_fraction = 0.01;
        context.shared_memory.num_threads = 1;
        context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;

        chg = std::make_unique<ds::CommunityHypergraph>(hypergraph, context, true);
    }

    void verifyInitialProbabilities(
        const ds::CommunityHypergraph& chyper,
        const parallel::scalable_vector<parallel::AtomicWrapper<HyperedgeWeight>>& community_volumes,
        const parallel::scalable_vector<parallel::AtomicWrapper<double>> community_exit_probability_mul_vol_total,
        const bool log = false) {
        // testing if the probability of being in any community is equal to 1
        HyperedgeWeight sum = 0;
        for (const auto& v : community_volumes) {
            sum += v;
        }
        if (log) LOG << V(sum) << V(chyper.totalVolume());
        ASSERT_EQ(sum, chyper.totalVolume());

        // testing if the exit probability is always less than or equal 1
        double sum_of_exit_probs_mul_vol_total = 0.0;
        for (const auto& q : community_exit_probability_mul_vol_total) {
            if (log) LOG << V(q) << V(chyper.totalVolume());
            if (log) LOG << static_cast<double>(q) / chyper.totalVolume();
            sum_of_exit_probs_mul_vol_total += q;
            ASSERT_LE(q, chyper.totalVolume());
        }
        ASSERT_LE(sum_of_exit_probs_mul_vol_total, chyper.totalVolume());

    }

    void verifyInitialProbabilitiesExactly(
        const parallel::scalable_vector<parallel::AtomicWrapper<double>>& community_exit_probability_mul_vol_total,
        const parallel::scalable_vector<double>& expected,
        const bool log = false) {
        for (size_t i = 0; i < community_exit_probability_mul_vol_total.size(); ++i) {
            if (log) LOG << V(expected[i]) << V(community_exit_probability_mul_vol_total[i]);
            ASSERT_EQ(expected[i], community_exit_probability_mul_vol_total[i]);
        }
    }

    using Base::hypergraph;
    std::unique_ptr<ds::CommunityHypergraph> chg;
    Context context;
};

TEST_F(AHypergraphLocalMovingMapEquation, InitializesTheExitProbabilities) {
    HypergraphLocalMovingMapEquation hlmme(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmme.initializeCommunityVolumes(*chg, communities);
    verifyInitialProbabilities(*chg, hlmme._community_volumes, hlmme._community_exit_probability_mul_vol_total);
    // parallel::scalable_vector<double> expected = { 4,3,3,5,5,2,4 };
    // verifyInitialProbabilitiesExactly(hlmme._community_exit_probability_mul_vol_total, expected);
}

TEST_F(AHypergraphLocalMovingMapEquation, InitializesTheExitProbabilities0) {
    ds::StaticHypergraph hg = io::readHypergraphFile("../tests/instances/karate_club.graph.hgr", 0);
    ds::CommunityHypergraph chyper(hg, context, true);
    HypergraphLocalMovingMapEquation hlmme(chyper, context);
    parallel::scalable_vector<HypernodeID> communities(chyper.initialNumNodes());
    for (size_t i = 0; i < chyper.initialNumNodes(); ++i) {
        communities[i] = i;
    }
    hlmme.initializeCommunityVolumes(chyper, communities);
    verifyInitialProbabilities(chyper, hlmme._community_volumes, hlmme._community_exit_probability_mul_vol_total);
}

TEST_F(AHypergraphLocalMovingMapEquation, InitializesTheExitProbabilities1) {
    ds::StaticHypergraph hg = io::readHypergraphFile("../tests/instances/powersim.mtx.hgr", 0);
    ds::CommunityHypergraph chyper(hg, context, true);
    HypergraphLocalMovingMapEquation hlmme(chyper, context);
    parallel::scalable_vector<HypernodeID> communities(chyper.initialNumNodes());
    for (size_t i = 0; i < chyper.initialNumNodes(); ++i) {
        communities[i] = i;
    }
    hlmme.initializeCommunityVolumes(chyper, communities);
    verifyInitialProbabilities(chyper, hlmme._community_volumes, hlmme._community_exit_probability_mul_vol_total);
}

TEST_F(AHypergraphLocalMovingMapEquation, InitializesTheExitProbabilities2) {
    ds::StaticHypergraph hg = io::readHypergraphFile("../tests/instances/sat14_atco_enc1_opt2_10_16.cnf.primal.hgr", 0);
    ds::CommunityHypergraph chyper(hg, context, true);
    HypergraphLocalMovingMapEquation hlmme(chyper, context);
    parallel::scalable_vector<HypernodeID> communities(chyper.initialNumNodes());
    for (size_t i = 0; i < chyper.initialNumNodes(); ++i) {
        communities[i] = i;
    }
    hlmme.initializeCommunityVolumes(chyper, communities);
    verifyInitialProbabilities(chyper, hlmme._community_volumes, hlmme._community_exit_probability_mul_vol_total);
}

TEST_F(AHypergraphLocalMovingMapEquation, InitializesTheExitProbabilities3) {
    ds::StaticHypergraph hg = io::readHypergraphFile("../tests/instances/test_instance.hgr", 0);
    ds::CommunityHypergraph chyper(hg, context, true);
    HypergraphLocalMovingMapEquation hlmme(chyper, context);
    parallel::scalable_vector<HypernodeID> communities(chyper.initialNumNodes());
    for (size_t i = 0; i < chyper.initialNumNodes(); ++i) {
        communities[i] = i;
    }
    hlmme.initializeCommunityVolumes(chyper, communities);
    verifyInitialProbabilities(chyper, hlmme._community_volumes, hlmme._community_exit_probability_mul_vol_total);
}

TEST_F(AHypergraphLocalMovingMapEquation, InitializesTheExitProbabilities4) {
    ds::StaticHypergraph hg = io::readHypergraphFile("../tests/instances/twocenters.hgr", 0);
    ds::CommunityHypergraph chyper(hg, context, true);
    HypergraphLocalMovingMapEquation hlmme(chyper, context);
    parallel::scalable_vector<HypernodeID> communities(chyper.initialNumNodes());
    for (size_t i = 0; i < chyper.initialNumNodes(); ++i) {
        communities[i] = i;
    }
    hlmme.initializeCommunityVolumes(chyper, communities);
    verifyInitialProbabilities(chyper, hlmme._community_volumes, hlmme._community_exit_probability_mul_vol_total);
}

TEST_F(AHypergraphLocalMovingMapEquation, ComputesTheBestNeighboringCommunity0) {
    HypergraphLocalMovingMapEquation hlmme(*chg, context);
    parallel::scalable_vector<HypernodeID> communities(chg->initialNumNodes());
    for (size_t i = 0; i < chg->initialNumNodes(); ++i) {
        communities[i] = i;
    }
    hlmme.initializeCommunityVolumes(*chg, communities);
    hlmme.recalculateAllExitProbabilities(*chg, communities);
    Probability before = metrics::hyp_map_equation(*chg, communities);
    const auto move = hlmme.calculateBestMove(*chg, communities, 0, hlmme._large_overlap_map.local(), hlmme._large_deltas_mul_vol_total.local());
    communities[0] = move;
    Probability after = metrics::hyp_map_equation(*chg, communities);
    ASSERT_EQ(move, 2);
    ASSERT_GT(before - after, 0);
    //ASSERT_DOUBLE_EQ(move.delta, before - after);
}

TEST_F(AHypergraphLocalMovingMapEquation, ComputesTheBestNeighboringCommunity1) {
    HypergraphLocalMovingMapEquation hlmme(*chg, context);
    parallel::scalable_vector<HypernodeID> communities(chg->initialNumNodes());
    for (size_t i = 0; i < chg->initialNumNodes(); ++i) {
        communities[i] = i;
    }
    hlmme.initializeCommunityVolumes(*chg, communities);
    hlmme.recalculateAllExitProbabilities(*chg, communities);
    Probability before = metrics::hyp_map_equation(*chg, communities);
    const auto move = hlmme.calculateBestMove(*chg, communities, 1, hlmme._large_overlap_map.local(), hlmme._large_deltas_mul_vol_total.local());
    communities[1] = move;
    Probability after = metrics::hyp_map_equation(*chg, communities);
    ASSERT_TRUE(move == 0 || move == 3 || move == 4);
    ASSERT_GT(before - after, 0);
    //ASSERT_DOUBLE_EQ(move.delta, before - after);
}
TEST_F(AHypergraphLocalMovingMapEquation, ComputesTheBestNeighboringCommunity2) {
    HypergraphLocalMovingMapEquation hlmme(*chg, context);
    parallel::scalable_vector<HypernodeID> communities(chg->initialNumNodes());
    for (size_t i = 0; i < chg->initialNumNodes(); ++i) {
        communities[i] = i;
    }
    hlmme.initializeCommunityVolumes(*chg, communities);
    hlmme.recalculateAllExitProbabilities(*chg, communities);
    Probability before = metrics::hyp_map_equation(*chg, communities);
    const auto move = hlmme.calculateBestMove(*chg, communities, 2, hlmme._large_overlap_map.local(), hlmme._large_deltas_mul_vol_total.local());
    communities[2] = move;
    Probability after = metrics::hyp_map_equation(*chg, communities);
    ASSERT_EQ(move, 0);
    ASSERT_GT(before - after, 0);
    //ASSERT_DOUBLE_EQ(move.delta, before - after);
}
TEST_F(AHypergraphLocalMovingMapEquation, ComputesTheBestNeighboringCommunity3) {
    HypergraphLocalMovingMapEquation hlmme(*chg, context);
    parallel::scalable_vector<HypernodeID> communities(chg->initialNumNodes());
    for (size_t i = 0; i < chg->initialNumNodes(); ++i) {
        communities[i] = i;
    }
    hlmme.initializeCommunityVolumes(*chg, communities);
    hlmme.recalculateAllExitProbabilities(*chg, communities);
    Probability before = metrics::hyp_map_equation(*chg, communities);
    const auto move = hlmme.calculateBestMove(*chg, communities, 3, hlmme._large_overlap_map.local(), hlmme._large_deltas_mul_vol_total.local());
    communities[3] = move;
    Probability after = metrics::hyp_map_equation(*chg, communities);
    ASSERT_EQ(move, 4);
    ASSERT_GT(before - after, 0);
    //ASSERT_DOUBLE_EQ(move.delta, before - after);
}
TEST_F(AHypergraphLocalMovingMapEquation, ComputesTheBestNeighboringCommunity4) {
    HypergraphLocalMovingMapEquation hlmme(*chg, context);
    parallel::scalable_vector<HypernodeID> communities(chg->initialNumNodes());
    for (size_t i = 0; i < chg->initialNumNodes(); ++i) {
        communities[i] = i;
    }
    hlmme.initializeCommunityVolumes(*chg, communities);
    hlmme.recalculateAllExitProbabilities(*chg, communities);
    Probability before = metrics::hyp_map_equation(*chg, communities);
    const auto move = hlmme.calculateBestMove(*chg, communities, 4, hlmme._large_overlap_map.local(), hlmme._large_deltas_mul_vol_total.local());
    communities[4] = move;
    Probability after = metrics::hyp_map_equation(*chg, communities);
    ASSERT_EQ(move, 3);
    ASSERT_GT(before - after, 0);
    //ASSERT_DOUBLE_EQ(move.delta, before - after);
}
TEST_F(AHypergraphLocalMovingMapEquation, ComputesTheBestNeighboringCommunity5) {
    HypergraphLocalMovingMapEquation hlmme(*chg, context);
    parallel::scalable_vector<HypernodeID> communities(chg->initialNumNodes());
    for (size_t i = 0; i < chg->initialNumNodes(); ++i) {
        communities[i] = i;
    }
    hlmme.initializeCommunityVolumes(*chg, communities);
    hlmme.recalculateAllExitProbabilities(*chg, communities);
    Probability before = metrics::hyp_map_equation(*chg, communities);
    const auto move = hlmme.calculateBestMove(*chg, communities, 5, hlmme._large_overlap_map.local(), hlmme._large_deltas_mul_vol_total.local());
    communities[5] = move;
    Probability after = metrics::hyp_map_equation(*chg, communities);
    ASSERT_EQ(move, 2);
    ASSERT_GT(before - after, 0);
    //ASSERT_DOUBLE_EQ(move.delta, before - after);
}
TEST_F(AHypergraphLocalMovingMapEquation, ComputesTheBestNeighboringCommunity6) {
    HypergraphLocalMovingMapEquation hlmme(*chg, context);
    parallel::scalable_vector<HypernodeID> communities(chg->initialNumNodes());
    LOG << chg->initialNumNodes();
    for (size_t i = 0; i < chg->initialNumNodes(); ++i) {
        communities[i] = i;
    }
    hlmme.initializeCommunityVolumes(*chg, communities);
    hlmme.recalculateAllExitProbabilities(*chg, communities);
    Probability before = metrics::hyp_map_equation(*chg, communities);
    const auto move = hlmme.calculateBestMove(*chg, communities, 6, hlmme._large_overlap_map.local(), hlmme._large_deltas_mul_vol_total.local());
    communities[6] = move;
    Probability after = metrics::hyp_map_equation(*chg, communities);
    ASSERT_EQ(move, 5);
    ASSERT_GT(before - after, 0);
    //ASSERT_DOUBLE_EQ(move.delta, before - after);
}
}
}