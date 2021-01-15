#include "gmock/gmock.h"

#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_modularity.h"
#include <mt-kahypar/io/hypergraph_io.h>
#include "mt-kahypar/partition/preprocessing/community_detection/local_moving_modularity.h"
#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_louvain.h"


using ::testing::Test;

namespace mt_kahypar {
namespace ds {

TEST(AHypergraphLocalMovingSpeed, TestsTheSpeedOfTheDeltaCalculation) {
    //StaticHypergraph hg = io::readHypergraphFile("../tests/instances/sat14_atco_enc1_opt2_10_16.cnf.primal.hgr", 0);
    //StaticHypergraph hg = io::readHypergraphFile("../tests/instances/powersim.mtx.hgr", 0);
    //StaticHypergraph hg = io::readHypergraphFile("../tests/instances/karate_club.graph.hgr", 0);
    //StaticHypergraph hg = io::readHypergraphFile("../tests/instances/af_4_k101.mtx.hgr", 0);
    StaticHypergraph hg = io::readHypergraphFile("../tests/instances/nlpkkt120.mtx.hgr", 0);
    //StaticHypergraph hg = io::readHypergraphFile("../tests/instances/sat14_q_query_3_L80_coli.sat.cnf.dual.hgr", 0);

    utils::Timer::instance().start_timer("comm_init", "Initialize Community Hypergraph");
    CommunityHypergraph chg(hg);
    utils::Timer::instance().stop_timer("comm_init");
    parallel::scalable_vector<HypernodeID> communities = community_detection::hypergraph_louvain(chg, true);
    // Volume min = *std::min_element(hlmm.ratios.begin(), hlmm.ratios.end());
    // Volume max = *std::max_element(hlmm.ratios.begin(), hlmm.ratios.end());
    // Volume avg =  std::accumulate(hlmm.ratios.begin(), hlmm.ratios.end(), 0.0L) / hlmm.ratios.size();
    //chg = chg.contract(communities);
    LOG << "Number of nodes: " << chg.initialNumNodes();
    //LOG << "Prune Counter: " << hlmm.pruneCounter;
    // LOG << "Average: " << avg;
    // LOG << "Min: " << min;
    // LOG << "Max: " << max;
    LOG << utils::Timer::instance(true);
    ASSERT_TRUE(true);
}
}
}