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

    explicit CommunityCount(const HypernodeID pincount, IteratorRange<IncidenceIterator> pins, const bool multipins = false) : _communities(pincount / 2) {
        _single_cut_communities.reserve(pincount);
        if (!multipins) {
            // since there are no duplicate pins in the initial Hypergraph
            // and every pin is in its own community
            for (const HypernodeID hn : pins) {
                _single_cut_communities.push_back(hn);
            }
        } else {
            HypernodeID prev = std::numeric_limits<HypernodeID>::max();
            for (const HypernodeID hn : pins) {
                if (prev != hn) {
                    _single_cut_communities.push_back(hn);
                    prev = hn;
                }
            }
        }
        end_of_single_cuts = _single_cut_communities.size();
    }

    ~CommunityCount() = default;

    // ! expects the community to be in the datastructure already, else undefined
    // ! decrements the unique node count for the given community. 
    // ! Community will be removed if the count is zero afterwards
    void removeFromCommunity(const PartitionID id) {
        if (!_communities.decrement(id)) {
            removeFromSingleCut(id);
        }
    }

    // ! increments the unique node count for the given community.
    // ! It will be added if it is not in the datastructure yet.
    void addToCommunity(const PartitionID id) {
        if (!_communities.increment(id)) {
            if (removeFromSingleCut(id)) {
                _communities.insert(id, 2);
            } else {
                ASSERT(end_of_single_cuts < _single_cut_communities.size(), V(end_of_single_cuts) << V(_single_cut_communities.size()));
                _single_cut_communities[end_of_single_cuts] = id;
                ++end_of_single_cuts;
            }
        }
    }

    // ! IteratorRange over all communities with single cuts
    IteratorRange<VectorIterator> singleCuts() const {
        ASSERT(_single_cut_communities.cbegin() + end_of_single_cuts <= _single_cut_communities.cend(), V(end_of_single_cuts) << V(_single_cut_communities.size()));
        return IteratorRange<VectorIterator>(_single_cut_communities.cbegin(), _single_cut_communities.cbegin() + end_of_single_cuts);
    }

    // ! IteratorRange over all communities with multiple cuts
    // TODO: cbegin() and cend()
    IteratorRange<MapIterator> multiCuts() const {
        return IteratorRange<MapIterator>(_communities.begin(), _communities.end());
    }

private:
    // TODO: This is O(n), improve?
    // ! Removes the community from the single cut datastructure
    bool removeFromSingleCut(const PartitionID id) {
        for (size_t i = 0; i < end_of_single_cuts; ++i) {
            if (_single_cut_communities[i] == id) {
                ASSERT(end_of_single_cuts > 0 && end_of_single_cuts <= _single_cut_communities.size(), V(end_of_single_cuts) << V(_single_cut_communities.size()));
                _single_cut_communities[i] = _single_cut_communities[--end_of_single_cuts];
                return true;
            }
        }
        return false;
    }

    //! index of the first invalid entry
    size_t end_of_single_cuts;

    //! communities which only contain a single node of this hyperedge
    std::vector<PartitionID> _single_cut_communities;

    //! communities which contain multiple nodes of this hyperedge
    HashMap _communities;

};

}
}