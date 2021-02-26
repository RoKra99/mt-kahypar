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

    // ! Contains buffers that are needed during multilevel contractions.
    // ! Struct is allocated on top level Communityhypergraph and passed to each contracted
    // ! hypergraph such that memory can be reused in consecutive contractions.
    struct TmpCommunityHypergraphBuffer {
        explicit TmpCommunityHypergraphBuffer(const HypernodeID num_hypernodes) {
            tmp_node_volumes.resize("Preprocessing", "tmp_community_volumes", num_hypernodes);
        }

        Array<parallel::AtomicWrapper<HyperedgeWeight>> tmp_node_volumes;
    };

public:

    using EdgeSizes = parallel::scalable_vector<size_t>;
    using HyperedgeIterator = typename Hypergraph::HyperedgeIterator;
    using IncidenceIterator = typename Hypergraph::IncidenceIterator;
    using IncidentNetsIterator = typename Hypergraph::IncidentNetsIterator;
    using EdgeSizeIterator = typename EdgeSizes::const_iterator;
    //using Map = RobinHoodMap<PartitionID, size_t>;
    using Map = HashMap<PartitionID, size_t, xxhash<uint32_t>>;
    using VectorIterator = typename std::vector<PartitionID>::const_iterator;
    using MapIterator = typename Map::Iterator;

    CommunityHypergraph(const Context& context) :
        _community_counts(),
        _hg(nullptr),
        _node_volumes(),
        _vol_v(0),
        _d_edge_weights(),
        _valid_edge_sizes(),
        _total_edge_weight(0),
        _context(context),
        _tmp_community_hypergraph_buffer(nullptr),
        _community_count_locks() {}

    explicit CommunityHypergraph(Hypergraph& hypergraph, const Context& context) :
        _community_counts(hypergraph.initialNumEdges()),
        _hg(&hypergraph),
        _vol_v(0),
        _total_edge_weight(0),
        _context(context),
        _tmp_community_hypergraph_buffer(nullptr) {
        tbb::parallel_invoke([&] {
            _node_volumes.resize("Preprocessing", "node_volumes", hypergraph.initialNumNodes(), true, true);
            }, [&] {
                _d_edge_weights.resize("Preprocessing", "d_edge_weights", hypergraph.maxEdgeSize() + 1, true, true);
            }, [&] {
                _community_count_locks.resize("Preprocessing", "community_count_locks", hypergraph.initialNumEdges());
            });

        computeAndSetTotalEdgeWeight();
        computeAndSetInitialVolumes();
        computeAndSetCommunityCounts();
        allocateTmpCommunityHypergraphBuffer();
        // size_t sum_of_d = 0;
        // for (const auto& d : _valid_edge_sizes) {
        //     sum_of_d += d;
        // }
        // std::cout << _valid_edge_sizes.size() << ',' << _valid_edge_sizes[_valid_edge_sizes.size() - 1] << ',' << sum_of_d << ',';
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

    // ! IteratorRange over all communities with single cuts of the hyperedge
    IteratorRange<VectorIterator> singleCuts(const HyperedgeID he) {
        _community_count_locks[he].lock();
        IteratorRange<VectorIterator> it = _community_counts[he]->singleCuts();
        _community_count_locks[he].unlock();
        return it;
    }

    // ! IteratorRange over all communities with multiple cuts
    // TODO: cbegin() and cend()
    IteratorRange<MapIterator> multiCuts(const HyperedgeID he) {
        _community_count_locks[he].lock();
        IteratorRange<MapIterator> it = _community_counts[he]->multiCuts();
        _community_count_locks[he].unlock();
        return it;
    }

    // ######################## Community ########################

    // ! increments the unique node count for the given community.
    // ! It will be added if it is not in the datastructure yet.
    void addCommunityToHyperedge(const HyperedgeID he, const PartitionID community) {
        _community_count_locks[he].lock();
        if (_community_counts[he]) {
            ASSERT(edgeSize(he) > _context.preprocessing.community_detection.hyperedge_size_caching_threshold);
            _community_counts[he]->addToCommunity(community);
        }
        _community_count_locks[he].unlock();
    }

    // ! expects the community to be in the datastructure already, else undefined
    // ! decrements the unique node count for the given community. 
    // ! Community will be removed if the count is zero afterwards
    void removeCommunityFromHyperedge(const HyperedgeID he, const PartitionID community) {
        _community_count_locks[he].lock();
        if (_community_counts[he]) {
            ASSERT(edgeSize(he) > _context.preprocessing.community_detection.hyperedge_size_caching_threshold);
            _community_counts[he]->removeFromCommunity(community);
        }
        _community_count_locks[he].unlock();
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

    // ! contracts the Hypergraph according to the given communities
    CommunityHypergraph contract(StaticHypergraph& hypergraph, parallel::scalable_vector<HypernodeID>& communities);

    // ! Weight of a vertex
    HypernodeWeight nodeWeight(const HypernodeID u) const {
        return _hg->nodeWeight(u);
    }

    // ! Total weight of hypergraph
    HypernodeWeight totalWeight() const {
        return _hg->totalWeight();
    }

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

    // ! initializes a community count datastructures for all edges above a certain edgesize
    void computeAndSetCommunityCounts() {
        for (const HyperedgeID& he : _hg->edges()) {
            _community_counts[he] = edgeSize(he) > _context.preprocessing.community_detection.hyperedge_size_caching_threshold
                ? std::make_unique<CommunityCount<Map>>(edgeSize(he), pins(he)) : std::unique_ptr<CommunityCount<Map>>(nullptr);
        }
    }

    // ! Allocate the temporary contraction buffer
    void allocateTmpCommunityHypergraphBuffer() {
        if (!_tmp_community_hypergraph_buffer) {
            _tmp_community_hypergraph_buffer = new TmpCommunityHypergraphBuffer(
                _hg->_num_hypernodes);
        }
    }

    // ! contains the cut communities for each hyperedge
    parallel::scalable_vector<std::unique_ptr<CommunityCount<Map>>> _community_counts;

    // ! Hypergraph this datastructure is wrapped around
    Hypergraph* _hg;

    // ! volume for each node
    Array<HyperedgeWeight> _node_volumes;

    // ! total volume
    HyperedgeWeight _vol_v;

    // ! summed edgeweight for each edgesize
    Array<HyperedgeWeight> _d_edge_weights;

    // ! contains the indexes to all edgeSizes which occur in the graph
    parallel::scalable_vector<size_t> _valid_edge_sizes;

    // ! sum of all edgeweights
    HyperedgeWeight _total_edge_weight;

    // ! contains Hyperparameters for the algorithms
    const Context& _context;

    // ! contains structures used in contaction to avoid allocations
    TmpCommunityHypergraphBuffer* _tmp_community_hypergraph_buffer;

    Array<SpinLock> _community_count_locks;

};
}
}