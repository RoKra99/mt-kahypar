#include "gmock/gmock.h"

#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_modularity.h"
#include <mt-kahypar/io/hypergraph_io.h>


using ::testing::Test;

namespace mt_kahypar {
namespace ds {

TEST(FindBestMove, Example1) {
    StaticHypergraph hg = io::readHypergraphFile("../../tests/instances/powersim.mtx.hgr", 0);
    utils::Timer::instance().start_timer("comm_init", "Initialize Community Hypergraph");
    CommunityHypergraph chg(hg);
    utils::Timer::instance().stop_timer("comm_init");
    community_detection::HypergraphLocalMovingModularity hlmm(chg);
    for (size_t k = 0; k < 10; ++k) {
        for (size_t i = 1; i <= 15838; ++i) {
            hlmm.makeMove(hlmm.calculateBestMove(i));
        }
    }
    math::fast_power(100,100);
    LOG << utils::Timer::instance(true);
    ASSERT_TRUE(true);
}
}
}