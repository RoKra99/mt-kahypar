#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_modularity.h";

namespace mt_kahypar::community_detection {

ds::Clustering hypergraph_local_moving_contract_recurse(ds::CommunityHypergraph& chg, HypergraphLocalMovingModularity& hlmm);

ds::Clustering hypergraph_louvain(ds::CommunityHypergraph& chg);
}