#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_modularity.h"
#include "mt-kahypar/partition/preprocessing/community_detection/hypergraph_local_moving_map_equation.h"

namespace mt_kahypar::community_detection {

template<typename HypergraphLocalMoving>
parallel::scalable_vector<HypernodeID> hypergraph_local_moving_contract_recurse(ds::CommunityHypergraph& chg, HypergraphLocalMoving& hlmm);

parallel::scalable_vector<HypernodeID> hypergraph_louvain(ds::CommunityHypergraph& chg, const Context& context, const bool deactivate_random = false);
}