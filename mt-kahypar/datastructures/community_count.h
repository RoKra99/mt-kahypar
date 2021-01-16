#pragma once

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/datastructures/hash_maps.h"

namespace mt_kahypar {
namespace ds {

template <typename HashMap>
class CommunityCount {

public:

    using IncidenceIterator = typename Hypergraph::IncidenceIterator;
    using VectorIterator = typename std::vector<PartitionID>::const_iterator;
    using MapIterator = typename HashMap::Iterator;

    CommunityCount() = default;

    // TODO: Initialization size of Hashmap
    explicit CommunityCount(HypernodeID pincount, IteratorRange<IncidenceIterator> pins) : _communities(pincount / 2) {
        _single_cut_communities.reserve(pincount);
        // since there are no duplicate pins in the initial Hypergraph
        // and every pin is in its own community
        for (const HypernodeID hn : pins) {
            _single_cut_communities.push_back(hn);
        }
    }

    ~CommunityCount() = default;

    // ! expects the community to be in the datastructure already, else undefined
    // ! decrements the unique node count for the given community. 
    // ! Community will be removed if the count is zero afterwards
    void removeFromCommunity(PartitionID id) {
        if (!_communities.decrement(id)) {
            removeFromSingleCut(id);
        }
    }

    // ! increments the unique node count for the given community.
    // ! It will be added if it is not in the datastructure yet.
    void addToCommunity(PartitionID id) {
        if (!_communities.increment(id)) {
            if (removeFromSingleCut(id)) {
                _communities.insert(id, 2);
            } else {
                _single_cut_communities.push_back(id);
            }
        }
    }

    void contract(const parallel::scalable_vector<HypernodeID>& communities) {
        // mapping the single cut communities to the new community ids
        for (size_t i = 0; i < _single_cut_communities.size(); ++i) {
            _single_cut_communities[i] = communities[_single_cut_communities[i]];
        }

        // mapping the mulitcuts to the new community ids
        // every former muli cut is now jsut a single cut due to contraction
        for (const auto& e : _communities) {
            _single_cut_communities.push_back(communities[e.first]);
        }

        _communities.clear();
    }

    // ! IteratorRange over all communities with single cuts
    IteratorRange<VectorIterator> singleCuts() {
        return IteratorRange<VectorIterator>(_single_cut_communities.cbegin(), _single_cut_communities.cend());
    }

    // ! IteratorRange over all communities with multiple cuts
    // TODO: cbegin() and cend()
    IteratorRange<MapIterator> multiCut() {
        return IteratorRange<MapIterator>(_communities.begin(), _communities.end());
    }

private:
    // TODO: This is O(n), improve?
    // ! Removes the community from the single cut datastructure
    bool removeFromSingleCut(PartitionID id) {
        for (size_t i = 0; i < _single_cut_communities.size(); ++i) {
            if (_single_cut_communities[i] == id) {
                _single_cut_communities[i] = _single_cut_communities[_single_cut_communities.size() - 1];
                _single_cut_communities.pop_back();
                return true;
            }
        }
        return false;
    }

    std::vector<PartitionID> _single_cut_communities;
    HashMap _communities;

};

}
}