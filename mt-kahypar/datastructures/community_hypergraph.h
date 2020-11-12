#pragma once

#include <atomic>
#include <type_traits>
#include <mutex>

#include "tbb/parallel_invoke.h"
#include "tbb/parallel_for.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/stl/thread_locals.h"

namespace mt_kahypar {
namespace ds {

class CommunityHypergraph {
public:

    using HyperedgeIterator = typename Hypergraph::HyperedgeIterator;
    using IncidenceIterator = typename Hypergraph::IncidenceIterator;

    CommunityHypergraph() = default;

    explicit CommunityHypergraph(Hypergraph& hypergraph) : _hg(&hypergraph), _vol_v(0), _total_edge_weight(0) {
        tbb::parallel_invoke([&] {
            _node_volumes.resize("Preprocessing", "node_volumes", hypergraph.initialNumNodes(), true, true);
                             }, [&] {
                                 _community_volumes.resize("Preprocessing", "community_volumes", hypergraph.initialNumNodes(), true, true);
                             });
        // Initialize the communities as Singletons
        for (const HypernodeID& hn : _hg->nodes()) {
            _hg->setCommunityID(hn, hn);
        };
        computeAndSetTotalEdgeWeight();
        computeAndSetInitialVolumes();
    }

    ~CommunityHypergraph() {
        freeInternalData();
    }
    // ######################## Volumes ########################

    // ! returns the volume (weighted degreee) of the given node
    HyperedgeWeight nodeVolume(const NodeID hn) const {
        return _node_volumes[hn];
    }

    // ! returns the volume (weighted degreee) of the given community
    HyperedgeWeight communityVolume(const PartitionID p) const {
        return _community_volumes[p];
    }

    // ! returns the volume (weighted degreee) of the whole hypergraph
    HyperedgeWeight totalVolume() const {
        return _vol_v;
    }


    // ######################## Iterators ########################

    // ! Returns an iterator over the set of active edges of the hypergraph
    IteratorRange<HyperedgeIterator> edges() const {
        return _hg->edges();
    }

    // ! Returns a range to loop over the pins of hyperedge e.
    IteratorRange<IncidenceIterator> pins(const HyperedgeID he) const {
        return _hg->pins(he);
    }


    // ######################## Community ########################

    // ! Community id which hypernode u is assigned to
    PartitionID communityID(const HypernodeID u) const {
        return _hg->communityID(u);
    }


    // ######################## Hyperedges ########################

    // ! Weight of a hyperedge
    HypernodeWeight edgeWeight(const HyperedgeID e) const {
        return _hg->edgeWeight(e);
    }

    // ! Sum of all edgeweights
    HyperedgeWeight totalEdgeWeight() const {
        return _total_edge_weight;
    }

    // ! Number of pins of a hyperedge
    HypernodeID edgeSize(const HyperedgeID e) const {
        return _hg->edgeSize(e);
    }

private:

    void freeInternalData() {
        parallel::parallel_free(_node_volumes, _community_volumes);
        _vol_v = 0;
        _total_edge_weight = 0;
    }

    // ! computes the volumes (weighted degrees) of each node.
    // ! The volumes of the communities are the same since it starts with singleton communities
    // ! The total volume is the sum of all node volumes
    void computeAndSetInitialVolumes() {
        for (const HypernodeID& hn : _hg->nodes()) {
            for (const HyperedgeID& he : _hg->incidentEdges(hn)) {
                const HyperedgeWeight weight = _hg->edgeWeight(he);
                _node_volumes[hn] += weight;
                _community_volumes[hn] += weight;
                _vol_v += weight;
            }
        }
    }

    // ! computes and sets the total edgeweight
    void computeAndSetTotalEdgeWeight() {
        for (const HyperedgeID& hn : _hg->edges()) {
            _total_edge_weight += edgeWeight(hn);
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

    // ! sum of all edgeweights
    HyperedgeWeight _total_edge_weight;


};
}
}