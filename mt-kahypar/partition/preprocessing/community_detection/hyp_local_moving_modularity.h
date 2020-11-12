#pragma once

#include "mt-kahypar/datastructures/community_hypergraph.h"
#include "mt-kahypar/definitions.h"
#include <algorithm>

namespace mt_kahypar::metrics {
double hyp_modularity(const ds::CommunityHypergraph& hypergraph);
}

namespace mt_kahypar::community_detection {

class HypLocalMovingModularity {



public:
    HypLocalMovingModularity(ds::CommunityHypergraph& hypergraph) : _hg(&hypergraph), weight_to_community(hypergraph.initialNumNodes(), std::numeric_limits<double>::max()) {}

    // ! try to move the given node into a neighbouring community
    // ! returns the change in modularity (0.0 if the node was not moved)
    double tryMove(const HypernodeID v) {
        const PartitionID comm_v = _hg->communityID(v);
        // calculate delta for all neighbouring communities
        for (const HyperedgeID& he : _hg->incidentEdges(v)) {
            for (const HypernodeID& hn : _hg->pins(he)) {
                const PartitionID comm_hn = _hg->communityID(hn);
                // if the node hn belonges to a different community which we see for the first time
                if (hn != v && comm_hn == comm_v && weight_to_community[comm_hn] < std::numeric_limits<double>::max()) {
                    weight_to_community[comm_hn] = calculateDelta(v, comm_hn);
                    neigh_communities.emplace_back(comm_hn);
                }
            }
        }
        // find community with optimal(minimal) delta and reset the clearlist
        double best_delta = 0;
        PartitionID best_comm = comm_v;
        for (const PartitionID& p : neigh_communities) {
            if (weight_to_community[p] < best_delta) {
                best_delta = weight_to_community[p];
                best_comm = p;
            }
            weight_to_community[p] = std::numeric_limits<double>::max();
        }
        neigh_communities.clear();

        // make the actual move only when communities are different
        if (best_comm != comm_v) {
            _hg->setCommunityID(v, best_comm);
            _hg->addCommunityVolume(v, best_comm);
            _hg->subtractCommunityVolume(v, comm_v);
        }
        return best_delta;
    }



private:

    // ! calculates the change in modularity if the given node would be moved to the given community
    double calculateDelta(const HypernodeID v, const PartitionID destination) {
        const PartitionID comm_v = _hg->communityID(v);
        if (comm_v == destination) {
            return 0.0;
        }
        double edge_contribution = 0.0;
        for (const HyperedgeID& he : _hg->incidentEdges(v)) {
            bool pins_in_comm_v = false;
            bool pins_in_destination = false;
            for (const HypernodeID& hn : _hg->pins(he)) {
                const PartitionID comm_hn = _hg->communityID(hn);
                if (hn != v) {
                    pins_in_comm_v |= comm_hn == comm_v;
                    pins_in_destination |= comm_hn == destination;
                }
            }
            edge_contribution -= pins_in_comm_v ? _hg->edgeWeight(he) : 0;
            edge_contribution += pins_in_destination ? _hg->edgeWeight(he) : 0;
        }

        const HyperedgeWeight vol_total = _hg->totalVolume();
        const HyperedgeWeight vol_v = _hg->nodeVolume(v);
        const HyperedgeWeight vol_c = _hg->communityVolume(comm_v);
        const HyperedgeWeight vol_d = _hg->communityVolume(destination);
        double exp_edge_contribution = 0.0;
        for (HyperedgeID d = 0; d < _hg->maxEdgeSize(); ++d) {
            exp_edge_contribution += (static_cast<double>(_hg->dEdgeWeight(d)) / pow(vol_total, d)) * (pow(vol_total - vol_c + vol_v, d) - pow(vol_total - vol_c, d)
                                                                                                       + pow(vol_total - vol_d - vol_v, d) - pow(vol_total - vol_d, d));
        }
        return (edge_contribution + exp_edge_contribution) / _hg->totalEdgeWeight();
    }

    ds::CommunityHypergraph* _hg;

    // clearlist to find neighbouring community with best modularity value
    std::vector<double> weight_to_community;
    std::vector<PartitionID> neigh_communities;
};
}