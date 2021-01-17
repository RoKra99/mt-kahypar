#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_modularity.h"

namespace mt_kahypar::community_detection {

parallel::scalable_vector<HypernodeID> hypergraph_local_moving_contract_recurse(ds::CommunityHypergraph& chg, HypergraphLocalMovingModularity& hlmm);

parallel::scalable_vector<HypernodeID> hypergraph_louvain(ds::CommunityHypergraph& chg, const bool deactivate_random = false);
}