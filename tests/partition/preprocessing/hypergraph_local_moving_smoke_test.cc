#include "gmock/gmock.h"

#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_modularity.h"
#include <mt-kahypar/io/hypergraph_io.h>


using ::testing::Test;

namespace mt_kahypar {
namespace ds {

TEST(AHypergraphLocalMovingSpeed, TestsTheSpeedOfTheDeltaCalculation) {
    //StaticHypergraph hg = io::readHypergraphFile("../tests/instances/sat14_atco_enc1_opt2_10_16.cnf.primal.hgr", 0);
    StaticHypergraph hg = io::readHypergraphFile("../tests/instances/powersim.mtx.hgr", 0);
    //StaticHypergraph hg = io::readHypergraphFile("../tests/instances/karate_club.graph.hgr", 0);
    utils::Timer::instance().start_timer("comm_init", "Initialize Community Hypergraph");
    CommunityHypergraph chg(hg);
    utils::Timer::instance().stop_timer("comm_init");
    community_detection::HypergraphLocalMovingModularity hlmm(chg);
    //for (size_t k = 0; k < 1000; ++k) {
    bool moved = true;
    size_t counter = 0;
    while (moved) {
        moved = false;
        for (HypernodeID hn = 0; hn < hg.initialNumNodes(); ++hn) {
            moved |= hlmm.makeMove(hlmm.calculateBestMove(hn));
        }
        ++counter;
    }
    //Volume min = *std::min_element(hlmm.ratios.begin(), hlmm.ratios.end());
    //Volume max = *std::max_element(hlmm.ratios.begin(), hlmm.ratios.end());
    //Volume avg =  std::accumulate(hlmm.ratios.begin(), hlmm.ratios.end(), 0.0L) / hlmm.ratios.size();
    //int add = 1 + 1;
    //hlmm.ratios.clear();
//}
    LOG << utils::Timer::instance(true);
    //LOG << "Pruning count" << hlmm.pruned;
    LOG << "rounds until nothing moved: " << counter;
    ASSERT_TRUE(true);
}
}
}