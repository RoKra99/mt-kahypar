#include "gmock/gmock.h"

#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_modularity.h"
#include <mt-kahypar/io/hypergraph_io.h>


using ::testing::Test;

namespace mt_kahypar {
namespace ds {

TEST(AHypergraphLocalMovingSpeed, TestsTheSpeedOfTheDeltaCalculation) {
    StaticHypergraph hg = io::readHypergraphFile("../tests/instances/powersim.mtx.hgr", 0);
    utils::Timer::instance().start_timer("comm_init", "Initialize Community Hypergraph");
    CommunityHypergraph chg(hg);
    utils::Timer::instance().stop_timer("comm_init");
    community_detection::HypergraphLocalMovingModularity hlmm(chg);
    for (size_t k = 0; k < 1000; ++k) {
        for (size_t i = 1; i <= 15838; ++i) {
            hlmm.makeMove(hlmm.calculateBestMove(i));
        }
    }
    LOG << utils::Timer::instance(true);
    ASSERT_TRUE(true);
}

TEST(ReusableExponentiation, CalculatesThePowersCorrectly) {
    StaticHypergraph hg = io::readHypergraphFile("../tests/instances/powersim.mtx.hgr", 0);
    CommunityHypergraph chg(hg);
    const HyperedgeWeight maxEdgeSize = 40;
    std::vector<Volume> powers(maxEdgeSize + 1, 0.0L);
    powers[0] = 1.0L;
    PartitionID biggest_d_yet = 0;
    for (const PartitionID d : chg.edgeSizes()) {
        ASSERT_DOUBLE_EQ(math::fast_power(0.9L, d), math::reusable_power(powers, biggest_d_yet, 0.9L, d));
        biggest_d_yet = d;
    }
}
}
}