#include "gmock/gmock.h"

#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_map_equation.h"
#include "mt-kahypar/datastructures/static_hypergraph_factory.h"
#include "tests/datastructures/hypergraph_fixtures.h"
#include <mt-kahypar/io/hypergraph_io.h>

using ::testing::Test;

namespace mt_kahypar {
namespace community_detection {

class AHyperGraphLocalMovingMapEquation : public ds::HypergraphFixture<Hypergraph, HypergraphFactory> {
    using Base = ds::HypergraphFixture<Hypergraph, HypergraphFactory>;

public:
    AHyperGraphLocalMovingMapEquation() :
        Base(),
        chg(nullptr),
        context() {
        context.preprocessing.community_detection.max_pass_iterations = 100;
        context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
        context.shared_memory.num_threads = 1;
        context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;

        chg = std::make_unique<ds::CommunityHypergraph>(hypergraph, context);
    }

    void verifyInitialProbabilities(
        const ds::CommunityHypergraph& chyper,
        const parallel::scalable_vector<parallel::AtomicWrapper<HyperedgeWeight>>& community_volumes,
        const parallel::scalable_vector<parallel::AtomicWrapper<uint32_t>> community_exit_probability_mul_vol_total,
        const bool log = false) {
        // testing if the probability of being in any community is equal to 1
        HyperedgeWeight sum = 0;
        for (const auto& v : community_volumes) {
            sum += v;
        }
        if (log) LOG << V(sum) << V(chyper.totalVolume());
        ASSERT_EQ(sum, chyper.totalVolume());

        // testing if the exit probability is always less than 1
        for (const auto& q : community_exit_probability_mul_vol_total) {
            if (log) LOG << V(q) << V(chyper.totalVolume());
            ASSERT_LT(q, chyper.totalVolume());
        }

        // testing if the probability to use the inter community book is less than 1
        double prob_to_use_inter_community_codebook = 0.0;
        for (size_t i = 0; i < chyper.initialNumNodes(); ++i) {
            const double prob_in_i_and_exiting = static_cast<double>(community_exit_probability_mul_vol_total[i] * community_volumes[i]) / (chyper.totalVolume() * chyper.totalVolume());
            if (log) LOG << V(prob_in_i_and_exiting);
            ASSERT_LT(prob_in_i_and_exiting, 1.0);
            prob_to_use_inter_community_codebook += prob_in_i_and_exiting;
        }
        if (log) LOG << V(prob_to_use_inter_community_codebook);
        ASSERT_LT(prob_to_use_inter_community_codebook, 1.0);
    }

    void verifyInitialProbabilitiesExactly(
        const parallel::scalable_vector<parallel::AtomicWrapper<uint32_t>>& community_exit_probability_mul_vol_total,
        const parallel::scalable_vector<uint32_t>& expected,
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

TEST_F(AHyperGraphLocalMovingMapEquation, InitializesTheExitProbabilities) {
    HypergraphLocalMovingMapEquation hlmme(*chg, context);
    parallel::scalable_vector<HypernodeID> communities = { 0,1,2,3,4,5,6 };
    hlmme.initializeCommunityVolumes(*chg, communities);
    verifyInitialProbabilities(*chg, hlmme._community_volumes, hlmme._community_exit_probability_mul_vol_total);
    parallel::scalable_vector<uint32_t> expected = { 4,3,3,5,5,2,4 };
    verifyInitialProbabilitiesExactly(hlmme._community_exit_probability_mul_vol_total, expected);
    LOG << metrics::hyp_map_equation(*chg, communities);
}

TEST_F(AHyperGraphLocalMovingMapEquation, InitializesTheExitProbabilities0) {
    ds::StaticHypergraph hg = io::readHypergraphFile("../tests/instances/karate_club.graph.hgr", 0);
    ds::CommunityHypergraph chyper(hg, context);
    HypergraphLocalMovingMapEquation hlmme(chyper, context);
    parallel::scalable_vector<HypernodeID> communities(chyper.initialNumNodes());
    for (size_t i = 0; i < chyper.initialNumNodes(); ++i) {
        communities[i] = i;
    }
    hlmme.initializeCommunityVolumes(chyper, communities);
    verifyInitialProbabilities(chyper, hlmme._community_volumes, hlmme._community_exit_probability_mul_vol_total);
    LOG << metrics::hyp_map_equation(chyper, communities);
}

TEST_F(AHyperGraphLocalMovingMapEquation, InitializesTheExitProbabilities1) {
    ds::StaticHypergraph hg = io::readHypergraphFile("../tests/instances/powersim.mtx.hgr", 0);
    ds::CommunityHypergraph chyper(hg, context);
    HypergraphLocalMovingMapEquation hlmme(chyper, context);
    parallel::scalable_vector<HypernodeID> communities(chyper.initialNumNodes());
    for (size_t i = 0; i < chyper.initialNumNodes(); ++i) {
        communities[i] = i;
    }
    hlmme.initializeCommunityVolumes(chyper, communities);
    verifyInitialProbabilities(chyper, hlmme._community_volumes, hlmme._community_exit_probability_mul_vol_total);
    LOG << metrics::hyp_map_equation(chyper, communities);
}

TEST_F(AHyperGraphLocalMovingMapEquation, InitializesTheExitProbabilities2) {
    ds::StaticHypergraph hg = io::readHypergraphFile("../tests/instances/sat14_atco_enc1_opt2_10_16.cnf.primal.hgr", 0);
    ds::CommunityHypergraph chyper(hg, context);
    HypergraphLocalMovingMapEquation hlmme(chyper, context);
    parallel::scalable_vector<HypernodeID> communities(chyper.initialNumNodes());
    for (size_t i = 0; i < chyper.initialNumNodes(); ++i) {
        communities[i] = i;
    }
    hlmme.initializeCommunityVolumes(chyper, communities);
    verifyInitialProbabilities(chyper, hlmme._community_volumes, hlmme._community_exit_probability_mul_vol_total);
    LOG << metrics::hyp_map_equation(chyper, communities);
}

TEST_F(AHyperGraphLocalMovingMapEquation, InitializesTheExitProbabilities3) {
    ds::StaticHypergraph hg = io::readHypergraphFile("../tests/instances/test_instance.hgr", 0);
    ds::CommunityHypergraph chyper(hg, context);
    HypergraphLocalMovingMapEquation hlmme(chyper, context);
    parallel::scalable_vector<HypernodeID> communities(chyper.initialNumNodes());
    for (size_t i = 0; i < chyper.initialNumNodes(); ++i) {
        communities[i] = i;
    }
    hlmme.initializeCommunityVolumes(chyper, communities);
    verifyInitialProbabilities(chyper, hlmme._community_volumes, hlmme._community_exit_probability_mul_vol_total);
    LOG << metrics::hyp_map_equation(chyper, communities);
}

TEST_F(AHyperGraphLocalMovingMapEquation, InitializesTheExitProbabilities4) {
    ds::StaticHypergraph hg = io::readHypergraphFile("../tests/instances/twocenters.hgr", 0);
    ds::CommunityHypergraph chyper(hg, context);
    HypergraphLocalMovingMapEquation hlmme(chyper, context);
    parallel::scalable_vector<HypernodeID> communities(chyper.initialNumNodes());
    for (size_t i = 0; i < chyper.initialNumNodes(); ++i) {
        communities[i] = i;
    }
    hlmme.initializeCommunityVolumes(chyper, communities);
    verifyInitialProbabilities(chyper, hlmme._community_volumes, hlmme._community_exit_probability_mul_vol_total);
    LOG << metrics::hyp_map_equation(chyper, communities);
}
}
}