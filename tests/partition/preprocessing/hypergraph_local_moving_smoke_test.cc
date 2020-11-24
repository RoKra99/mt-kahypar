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
    for (size_t k = 0; k < 10; ++k) {
        for (HypernodeID hn = 0; hn < 15838; ++hn) {
            hlmm.makeMove(hlmm.calculateBestMove(hn));
        }
    }
    LOG << utils::Timer::instance(true);
    ASSERT_TRUE(true);
}
}
}