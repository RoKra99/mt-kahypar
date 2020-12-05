#pragma once

#include "tbb/parallel_invoke.h"
#include "tbb/parallel_for.h"

#include "mt-kahypar/definitions.h"

namespace mt_kahypar {
namespace ds {

class CommunityHypergraph {
public:

    using CommunityVolumes = Array<HyperedgeWeight>;
    using EdgeSizes = std::vector<PartitionID>;
    using HyperedgeIterator = typename Hypergraph::HyperedgeIterator;
    using IncidenceIterator = typename Hypergraph::IncidenceIterator;
    using IncidentNetsIterator = typename Hypergraph::IncidentNetsIterator;
    using CommunityVolumeIterator = typename CommunityVolumes::const_iterator;
    using EdgeSizeIterator = typename EdgeSizes::const_iterator;

    CommunityHypergraph() = default;

    explicit CommunityHypergraph(Hypergraph& hypergraph) : _hg(&hypergraph), _vol_v(0), _total_edge_weight(0), _max_d_edge_weight_(0) {
        tbb::parallel_invoke([&] {
            _node_volumes.resize("Preprocessing", "node_volumes", hypergraph.initialNumNodes(), true, true);
                             }, [&] {
                                 _community_volumes.resize("Preprocessing", "community_volumes", hypergraph.initialNumNodes(), true, true);
                             }, [&] {
                                 _d_edge_weights.resize("Preprocessing", "d_edge_weights", hypergraph.maxEdgeSize() + 1, true, true);
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
    HyperedgeWeight nodeVolume(const HypernodeID hn) const {
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

    // ! sets the volume of a community to the spectified value
    void setCommunityVolume(const PartitionID community_id, const HyperedgeWeight volume) {
        _community_volumes[community_id] = volume;
    }

    // ! adds the volume of the node to the volume of the community
    void addCommunityVolume(const HypernodeID v, const PartitionID community_id) {
        _community_volumes[community_id] += _node_volumes[v];
    }

    // ! subtracts the volume of the node from the volume of the community
    void subtractCommunityVolume(const HypernodeID v, const PartitionID community_id) {
        _community_volumes[community_id] -= _node_volumes[v];
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

    // ! Returns a range to loop over the incident nets of hypernode u.
    IteratorRange<IncidentNetsIterator> incidentEdges(const HypernodeID u) const {
        return _hg->incidentEdges(u);
    }

    // ! Returns a range of the community volumes
    IteratorRange<CommunityVolumeIterator> communityVolumes() const {
        return IteratorRange<CommunityVolumeIterator>(_community_volumes.cbegin(), _community_volumes.cend());
    }

    // ! Returns a range of the occuring edge sizes
    IteratorRange<EdgeSizeIterator> edgeSizes() const {
        return IteratorRange<EdgeSizeIterator>(_valid_edge_sizes.cbegin(), _valid_edge_sizes.cend());
    }


    // ######################## Community ########################

    // ! Community id which hypernode u is assigned to
    PartitionID communityID(const HypernodeID u) const {
        return _hg->communityID(u);
    }

    // ! Assign a community to a hypernode
    void setCommunityID(const HypernodeID u, const PartitionID community_id) {
        _hg->setCommunityID(u, community_id);
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

    // ! Sum of all edweights of size d
    HyperedgeWeight edgeWeightBySize(const HypernodeID d) const {
        return _d_edge_weights[d];
    }

    // ! Number of pins of a hyperedge
    HypernodeID edgeSize(const HyperedgeID e) const {
        return _hg->edgeSize(e);
    }

    // ! Maximum size of a hyperedge
    HypernodeID maxEdgeSize() const {
        return _hg->maxEdgeSize();
    }

    // ! Maximum edgeweight accumulated by edgesize
    HypernodeID maxAccumulatedEdgeWeight() {
        return _max_d_edge_weight_;
    }


    // ######################## Nodes ########################

    // ! Initial number of hypernodes
    HypernodeID initialNumNodes() const {
        return _hg->initialNumNodes();
    }

private:

    void freeInternalData() {
        parallel::parallel_free(_node_volumes, _community_volumes, _d_edge_weights);
        _vol_v = 0;
        _total_edge_weight = 0;
        _max_d_edge_weight_ = 0;
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
            _d_edge_weights[edgeSize(hn)] += edgeWeight(hn);
        }
        for (size_t i = 0; i < _d_edge_weights.size(); ++i) {
            if (_d_edge_weights[i] > 0) {
                _valid_edge_sizes.emplace_back(i);
                _max_d_edge_weight_ = _d_edge_weights[i] > _max_d_edge_weight_ ? _d_edge_weights[i] : _max_d_edge_weight_;
            }
        }
    }

    // ! Hypergraph this datastructure is wrapped around
    Hypergraph* _hg;

    // ! volume for each node
    Array<HyperedgeWeight> _node_volumes;

    // ! volume for each community
    Array<HyperedgeWeight> _community_volumes;

    // ! total volume
    HyperedgeWeight _vol_v;

    // ! summed edgeweight for each edgesize
    Array<HyperedgeWeight> _d_edge_weights;

    // ! contains the indexes to all edgeSizes which occur in the graph
    std::vector<PartitionID> _valid_edge_sizes;

    // ! sum of all edgeweights
    HyperedgeWeight _total_edge_weight;

    // ! maximum value of _d_edge_weight 
    HyperedgeWeight _max_d_edge_weight_;


};
}
}