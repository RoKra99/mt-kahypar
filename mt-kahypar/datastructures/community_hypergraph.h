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
        explicit TmpCommunityHypergraphBuffer(const HypernodeID num_hypernodes, const HypernodeID num_pins) {
            tbb::parallel_invoke([&] {
                tmp_node_volumes.resize("Preprocessing", "tmp_community_volumes", num_hypernodes);
            }, [&] {
                multipin_mapping.resize("Preprocessing", "multipin_mapping", num_pins);
            });
        }

        Array<parallel::AtomicWrapper<HyperedgeWeight>> tmp_node_volumes;
        Array<HypernodeID> multipin_mapping;
    };

    struct EdgeSizeObject {
        const size_t index;
        const size_t remaining_d;
        const Volume weight;
    };

public:

    using EdgeSizes = parallel::scalable_vector<EdgeSizeObject>;
    using HyperedgeIterator = typename Hypergraph::HyperedgeIterator;
    using IncidenceIterator = typename Hypergraph::IncidenceIterator;
    using IncidentNetsIterator = typename Hypergraph::IncidentNetsIterator;
    using MultipinIterator = typename Array<Multipin>::const_iterator;
    using EdgeSizeIterator = typename EdgeSizes::const_iterator;
    using Map = HashMap<PartitionID, size_t, xxhash<uint32_t>>;
    //using Map = FixedSizeSparseMap<PartitionID, size_t>;
    using VectorIterator = typename std::vector<PartitionID>::const_iterator;
    using MapIterator = typename Map::Iterator;

    CommunityHypergraph(const Context& context, const bool use_multipins = false) :
        _use_multipins(use_multipins),
        _multipin_incidence_array(),
        _multipin_indexes(),
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

    explicit CommunityHypergraph(Hypergraph& hypergraph, const Context& context, const bool use_multipins = false) :
        _use_multipins(use_multipins),
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
        if (use_multipins) {
            initializeMultipinDataStructures();
        }
        computeAndSetTotalEdgeWeight();
        computeAndSetInitialVolumes();
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

    // ! Returns a range to loop over the pins and multiplicities of hyperedge e.
    IteratorRange<MultipinIterator> multipins(const HyperedgeID e) const {
        ASSERT(!_hg->hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
        ASSERT(e < _multipin_indexes.size());
        ASSERT(e + 1 < _multipin_indexes.size());
        const size_t start = _multipin_indexes[e];
        const size_t end = _multipin_indexes[e + 1];
        return IteratorRange<MultipinIterator>(
            _multipin_incidence_array.cbegin() + start,
            _multipin_incidence_array.cbegin() + end);
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
        ASSERT(_community_counts[he]);
        _community_count_locks[he].lock();
        IteratorRange<VectorIterator> it = _community_counts[he]->singleCuts();
        _community_count_locks[he].unlock();
        return it;
    }

    // ! IteratorRange over all communities with multiple cuts
    IteratorRange<MapIterator> multiCuts(const HyperedgeID he) {
        ASSERT(_community_counts[he]);
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

private:

    void freeInternalData() {
        //parallel::parallel_free(_node_volumes, _d_edge_weights, _multipin_incidence_array, _multipin_indexes, _community_count_locks);
        if (_tmp_community_hypergraph_buffer) {
            delete(_tmp_community_hypergraph_buffer);
        }
    }

    // ! computes the volumes (weighted degrees) of each node.
    // ! The volumes of the communities are the same since it starts with singleton communities
    // ! The total volume is the sum of all node volumes
    void computeAndSetInitialVolumes() {
        tbb::enumerable_thread_specific<HyperedgeWeight> local_total_volume(0);
        _hg->doParallelForAllNodes([&](const HypernodeID hn) {
            for (const HyperedgeID& he : _hg->incidentEdges(hn)) {
                const HyperedgeWeight weight = _hg->edgeWeight(he);
                _node_volumes[hn] += weight;
                local_total_volume.local() += weight;
            }
        });
        _vol_v = local_total_volume.combine(std::plus<HyperedgeWeight>());
        ASSERT(_vol_v > 0);
    }

    // ! computes and sets the total edgeweight
    void computeAndSetTotalEdgeWeight() {
        tbb::enumerable_thread_specific<HyperedgeWeight> local_total_edge_weight(0);
        _hg->doParallelForAllEdges([&](const HyperedgeID he) {
            const HyperedgeWeight weight = _hg->edgeWeight(he);
            _d_edge_weights[edgeSize(he)] += weight;
            local_total_edge_weight.local() += weight;
        });
        _total_edge_weight = local_total_edge_weight.combine(std::plus<HyperedgeWeight>());
        size_t index = 0;
        size_t biggest_d_yet = 1;
        for (size_t i = 0; i < _d_edge_weights.size(); ++i) {
            if (_d_edge_weights[i] > 0) {
                const size_t remaining_d = i - biggest_d_yet;
                biggest_d_yet = i;
                _valid_edge_sizes.push_back({ index, remaining_d, static_cast<Volume>(_d_edge_weights[i]) });
                ++index;
            }
        }
    }

    // ! initializes a community count datastructures for all edges above a certain edgesize
    void computeAndSetCommunityCounts() {
        if (_use_multipins) {
            _hg->doParallelForAllEdges([&](const HyperedgeID he) {
                _community_counts[he] = edgeSize(he) > _context.preprocessing.community_detection.hyperedge_size_caching_threshold
                    ? std::make_unique<CommunityCount<Map>>(edgeSize(he), multipins(he)) : std::unique_ptr<CommunityCount<Map>>(nullptr);
            });
        } else {
            _hg->doParallelForAllEdges([&](const HyperedgeID he) {
                _community_counts[he] = edgeSize(he) > _context.preprocessing.community_detection.hyperedge_size_caching_threshold
                    ? std::make_unique<CommunityCount<Map>>(edgeSize(he), pins(he)) : std::unique_ptr<CommunityCount<Map>>(nullptr);
            });
        }
    }

    // ! Allocate the temporary contraction buffer
    void allocateTmpCommunityHypergraphBuffer() {
        if (!_tmp_community_hypergraph_buffer) {
            _tmp_community_hypergraph_buffer = new TmpCommunityHypergraphBuffer(
                _hg->_num_hypernodes, _hg->_num_pins);
        }
    }

    // ! Initializes the incidence and index array for multipins
    void initializeMultipinDataStructures() {
        tbb::parallel_invoke([&] {
            _multipin_incidence_array.resize("Preprocessing", "multipin_incidence_array", _hg->initialNumPins());
        }, [&] {
            _multipin_indexes.resize("Preprocessing", "multipin_indexes", _hg->initialNumEdges() + 1);
        });

        tbb::parallel_for(0U, _hg->initialNumPins(), [&](const HypernodeID hn) {
            _multipin_incidence_array[hn].id = _hg->_incidence_array[hn];
            _multipin_incidence_array[hn].multiplicity = 1U;
        });

        tbb::parallel_for(0U, _hg->initialNumEdges(), [&](const HyperedgeID he) {
            _multipin_indexes[he] = _hg->hyperedge(he).firstEntry();
        });
        _multipin_indexes[_hg->initialNumEdges()] = _hg->initialNumPins();
    }

    const bool _use_multipins;

    Array<Multipin> _multipin_incidence_array;

    Array<size_t> _multipin_indexes;

    // ! contains the cut communities for each hyperedge
    parallel::scalable_vector<std::unique_ptr<CommunityCount<Map>>> _community_counts;

    // ! Hypergraph this datastructure is wrapped around
    Hypergraph* _hg;

    // ! volume for each node
    Array<HyperedgeWeight> _node_volumes;

    // ! total volume
    HyperedgeWeight _vol_v;

    // ! summed edgeweight for each edgesize
    Array<parallel::AtomicWrapper<HyperedgeWeight>> _d_edge_weights;

    // ! contains the indexes to all edgeSizes which occur in the graph
    EdgeSizes _valid_edge_sizes;

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