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
    StaticHypergraph hg = io::readHypergraphFile("../tests/instances/af_4_k101.mtx.hgr", 0);
    //StaticHypergraph hg = io::readHypergraphFile("../tests/instances/nlpkkt120.mtx.hgr", 0);
    //StaticHypergraph hg = io::readHypergraphFile("../tests/instances/sat14_q_query_3_L80_coli.sat.cnf.dual.hgr", 0);
    Context context;
    context.preprocessing.community_detection.max_pass_iterations = 100;
    context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
    context.shared_memory.num_threads = 1;
    utils::Timer::instance().start_timer("comm_init", "Initialize Community Hypergraph");
    CommunityHypergraph chg(hg);
    utils::Timer::instance().stop_timer("comm_init");
    parallel::scalable_vector<HypernodeID> communities = community_detection::hypergraph_louvain(chg, context);
    LOG << "community count " << *std::max_element(communities.begin(), communities.end());
    //chg = chg.contract(communities);
    LOG << "Number of nodes: " << chg.initialNumNodes();
    LOG << utils::Timer::instance(true);
    ASSERT_TRUE(true);
}
}
}