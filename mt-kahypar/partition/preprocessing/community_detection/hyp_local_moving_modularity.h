#pragma once

#include "mt-kahypar/datastructures/community_hypergraph.h"
#include "mt-kahypar/definitions.h"
#include <algorithm>

namespace mt_kahypar::metrics {
double hyp_modularity(const ds::CommunityHypergraph& hypergraph, const ds::Clustering& communities);
}

namespace mt_kahypar::community_detection {

class HypLocalMovingModularity {

    

    public:
    private:
};
}