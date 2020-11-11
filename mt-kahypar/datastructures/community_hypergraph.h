#pragma once

#include <atomic>
#include <type_traits>
#include <mutex>

#include "tbb/parallel_invoke.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/stl/thread_locals.h"

namespace mt_kahypar {
namespace ds {

class CommunityHypergraph {
public:
    CommunityHypergraph() = default;

    explicit CommunityHypergraph(Hypergraph& hypergraph) : _hg(&hypergraph), _vol_v(0) {
        tbb::parallel_invoke([&] {
            _node_volumes.resize("Preprocessing", "node_volumes", hypergraph.initialNumNodes(), true, true);
                             }, [&] {
                                 _community_volumes.resize("Preprocessing", "community_volumes", hypergraph.initialNumNodes(), true, true);
                             });
        computeAndSetInitialVolumes();
    }

    HyperedgeWeight nodeVolume(const NodeID n) {
        return _node_volumes[n];
    }

    HyperedgeWeight communityVolume(const PartitionID p) {
        return _community_volumes[p];
    }

    HyperedgeWeight totalVolume() {
        return _vol_v;
    }

private:

    // ! computes the volumes (weighted degrees) of each node.
    // ! The volumes of the communities are the same since it starts with singleton communities
    // ! The total volume is the sum of all node volumes
    void computeAndSetInitialVolumes() {
        for (const HypernodeID n : _hg->nodes()) {
            for (const HyperedgeID he : _hg->incidentEdges(n)) {
                const HyperedgeWeight weight = _hg->edgeWeight(he);
                _node_volumes[n] += weight;
                _community_volumes[n] += weight;
                _vol_v += weight;
            }
        }
    }

    // ! Hypergraph this datastructure is wrapped around
    Hypergraph* _hg;

    // ! volume for each node
    Array<HyperedgeWeight> _node_volumes;

    // ! volume for each community
    Array <HyperedgeWeight> _community_volumes;

    // ! total volume
    HyperedgeWeight _vol_v;


};
}
}