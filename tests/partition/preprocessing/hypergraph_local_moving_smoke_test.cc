#include "gmock/gmock.h"

#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_modularity.h"
#include <mt-kahypar/io/hypergraph_io.h>
#include "mt-kahypar/partition/preprocessing/community_detection/local_moving_modularity.h"
#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_louvain.h"
#include "mt-kahypar/datastructures/graph.h"


using ::testing::Test;

namespace mt_kahypar {
namespace ds {

// TEST(AHypergraphLocalMovingSpeed, TestsTheSpeedOfTheDeltaCalculation) {
//     //StaticHypergraph hg = io::readHypergraphFile("../tests/instances/sat14_atco_enc1_opt2_10_16.cnf.primal.hgr", 0);
//     //StaticHypergraph hg = io::readHypergraphFile("../tests/instances/powersim.mtx.hgr", 0);
//     //StaticHypergraph hg = io::readHypergraphFile("../tests/instances/karate_club.graph.hgr", 0);
//     StaticHypergraph hg = io::readHypergraphFile("../tests/instances/af_4_k101.mtx.hgr", 0);
//     //StaticHypergraph hg = io::readHypergraphFile("../tests/instances/nlpkkt120.mtx.hgr", 0);
//     //StaticHypergraph hg = io::readHypergraphFile("../tests/instances/sat14_q_query_3_L80_coli.sat.cnf.dual.hgr", 0);
//     Context context;
//     context.preprocessing.community_detection.max_pass_iterations = 100;
//     context.preprocessing.community_detection.min_vertex_move_fraction = 0.0001;
//     context.shared_memory.num_threads = 1;
//     context.preprocessing.community_detection.hyperedge_size_caching_threshold = 0;
//     utils::Timer::instance().start_timer("comm_init", "Initialize Community Hypergraph");
//     CommunityHypergraph chg(hg, context);
//     utils::Timer::instance().stop_timer("comm_init");
//     parallel::scalable_vector<HypernodeID> communities = community_detection::hypergraph_louvain(chg, context);
//     LOG << "community count " << *std::max_element(communities.begin(), communities.end());
//     LOG << "Number of nodes: " << chg.initialNumNodes();
//     LOG << utils::Timer::instance(true);
//     ASSERT_TRUE(true);
// }

TEST(AHypergraphLocalMovingCompare, ComparesPartitions) {
    std::string name = "opt1.mtx.hgr";
    Context context;
    std::vector<int> communities;
    io::readPartitionFile("../../partitions/graph/" + name + ".part2.epsilon0.03.seed0.KaHyPar", communities);
    LOG << communities.size();
    parallel::scalable_vector<HypernodeID> clustering;
    for (const auto i : communities) {
        clustering.push_back(i);
    }

    StaticHypergraph hg = io::readHypergraphFile("/home/robert/benchmarks/set_a/" + name,0);

    // //Graph graph(hg, LouvainEdgeWeight::uniform);
    // LOG << graph.numNodes();
    // LOG << metrics::modularity(graph, clustering);
    CommunityHypergraph graph(hg, context);
    LOG << graph.initialNumNodes();
    LOG << "graph" << std::fixed << std::setprecision(15) << metrics::hyp_modularity(graph, clustering);

    std::vector<int> communities2;
    io::readPartitionFile("../../partitions/p16/" + name + ".part2.epsilon0.03.seed0.KaHyPar", communities2);
    LOG << V(communities2.size());
    parallel::scalable_vector<HypernodeID> clustering2;
    for (const auto i : communities2) {
        clustering2.push_back(i);
    }

    StaticHypergraph hg2 = io::readHypergraphFile("/home/robert/benchmarks/set_a/" + name,0);
    //hg2.setCommunityIDs(std::move(clustering2));
    //Graph graph2(hg2, LouvainEdgeWeight::uniform);
    CommunityHypergraph graph2(hg2, context);
    LOG << graph2.initialNumNodes();
    LOG << "hypergraph" << std::fixed << std::setprecision(15) << metrics::hyp_modularity(graph2, clustering2);
}
}
}