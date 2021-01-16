#pragma once

#include "tbb/parallel_invoke.h"
#include "tbb/parallel_for.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/community_count.h"
#include "mt-kahypar/datastructures/hash_maps.h"
#include "mt-kahypar/datastructures/hash_functions.h"

#include "static_hypergraph.h"

namespace mt_kahypar {
namespace ds {

class CommunityHypergraph {

    //using Map = RobinHoodMap<PartitionID, size_t>;
    using Map = HashMap<PartitionID, size_t, xxhash<uint32_t>>;

    //!Contains buffers that are needed during multilevel contractions.
    //!Struct is allocated on top level Communityhypergraph and passed to each contracted
    //!hypergraph such that memory can be reused in consecutive contractions.
    struct TmpCommunityHypergraphBuffer {
        explicit TmpCommunityHypergraphBuffer(const HypernodeID num_hypernodes, const HyperedgeID num_hyperedges): tmp_community_counts(num_hyperedges) {
            tmp_node_volumes.resize("Preprocessing", "tmp_community_volumes", num_hypernodes);
        }

        Array<parallel::AtomicWrapper<HyperedgeWeight>> tmp_node_volumes;
        parallel::scalable_vector<std::unique_ptr<CommunityCount<Map>>> tmp_community_counts;
    };

public:

    using EdgeSizes = parallel::scalable_vector<PartitionID>;
    using HyperedgeIterator = typename Hypergraph::HyperedgeIterator;
    using IncidenceIterator = typename Hypergraph::IncidenceIterator;
    using IncidentNetsIterator = typename Hypergraph::IncidentNetsIterator;
    using EdgeSizeIterator = typename EdgeSizes::const_iterator;
    static constexpr size_t EDGESIZE_THRESHHOLD = 0;

    CommunityHypergraph() :
        _community_counts(),
        _hg(nullptr),
        _node_volumes(),
        _vol_v(0),
        _d_edge_weights(),
        _valid_edge_sizes(),
        _total_edge_weight(0),
        _tmp_community_hypergraph_buffer(nullptr) {}

    explicit CommunityHypergraph(Hypergraph& hypergraph) : _community_counts(hypergraph.initialNumEdges()), _hg(&hypergraph), _vol_v(0), _total_edge_weight(0), _tmp_community_hypergraph_buffer(nullptr) {
        tbb::parallel_invoke([&] {
            _node_volumes.resize("Preprocessing", "node_volumes", hypergraph.initialNumNodes(), true, true);
            }, [&] {
                _d_edge_weights.resize("Preprocessing", "d_edge_weights", hypergraph.maxEdgeSize() + 1, true, true);
            });

        computeAndSetTotalEdgeWeight();
        computeAndSetInitialVolumes();
        // TODO: Add to register_memory_pool
        computeAndSetCommunityCounts();
        allocateTmpCommunityHypergraphBuffer();
    }


    CommunityHypergraph(CommunityHypergraph&& other) = default;
    CommunityHypergraph& operator=(CommunityHypergraph&& other) = default;

    ~CommunityHypergraph() {
        freeInternalData();
    }
    // ######################## Volumes ########################

    // ! returns the volume (weighted degreee) of the given node
    HyperedgeWeight nodeVolume(const HypernodeID hn) const {
        return _node_volumes[hn];
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

    // ! Returns a range to loop over the incident nets of hypernode u.
    IteratorRange<IncidentNetsIterator> incidentEdges(const HypernodeID u) const {
        return _hg->incidentEdges(u);
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

    // ! Minimum size of a Hyperedge
    HypernodeID minEdgeSize() const {
        ASSERT(!_valid_edge_sizes.empty());
        return _valid_edge_sizes[0];
    }

    // ! Number of different edge sizes that occur in the Hypergraph
    size_t numberOfUniqueEdgeSizes() const {
        return _valid_edge_sizes.size();
    }

    // ! Initial number of hyperedges
    HyperedgeID initialNumEdges() const {
        return _hg->initialNumEdges();
    }


    // ######################## Nodes ########################

    // ! Initial number of hypernodes
    HypernodeID initialNumNodes() const {
        return _hg->initialNumNodes();
    }

    CommunityHypergraph contract(StaticHypergraph& hypergraph, parallel::scalable_vector<HypernodeID>& communities);

    CommunityHypergraph mapContractedVolumes(StaticHypergraph& hypergraph);

    // ! contains the cut communities for each hyperedge
    // ! TODO: Get this to work with Array
    parallel::scalable_vector<std::unique_ptr<CommunityCount<Map>>> _community_counts;

private:

    void freeInternalData() {
        parallel::parallel_free(_node_volumes, _d_edge_weights);
        _vol_v = 0;
        _total_edge_weight = 0;
        if (_tmp_community_hypergraph_buffer) {
            delete(_tmp_community_hypergraph_buffer);
        }
    }

    // ! computes the volumes (weighted degrees) of each node.
    // ! The volumes of the communities are the same since it starts with singleton communities
    // ! The total volume is the sum of all node volumes
    void computeAndSetInitialVolumes() {
        for (const HypernodeID& hn : _hg->nodes()) {
            for (const HyperedgeID& he : _hg->incidentEdges(hn)) {
                const HyperedgeWeight weight = _hg->edgeWeight(he);
                _node_volumes[hn] += weight;
                _vol_v += weight;
            }
        }
        ASSERT(_vol_v > 0);
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
            }
        }
    }

    void computeAndSetCommunityCounts() {
        for (const HyperedgeID& he : _hg->edges()) {
            _community_counts[he] = edgeSize(he) > EDGESIZE_THRESHHOLD ? std::make_unique<CommunityCount<Map>>(edgeSize(he), pins(he)) : std::unique_ptr<CommunityCount<Map>>(nullptr);
        }
    }

    //! Allocate the temporary contraction buffer
    void allocateTmpCommunityHypergraphBuffer() {
        if (!_tmp_community_hypergraph_buffer) {
            _tmp_community_hypergraph_buffer = new TmpCommunityHypergraphBuffer(
                _hg->_num_hypernodes, _hg->_num_hyperedges);
        }
    }

    // ! Hypergraph this datastructure is wrapped around
    Hypergraph* _hg;

    // ! volume for each node
    Array<HyperedgeWeight> _node_volumes;

    // ! total volume
    HyperedgeWeight _vol_v;

    // ! summed edgeweight for each edgesize
    Array<HyperedgeWeight> _d_edge_weights;

    // ! contains the indexes to all edgeSizes which occur in the graph
    parallel::scalable_vector<PartitionID> _valid_edge_sizes;

    // ! sum of all edgeweights
    HyperedgeWeight _total_edge_weight;

    TmpCommunityHypergraphBuffer* _tmp_community_hypergraph_buffer;

};
}
}