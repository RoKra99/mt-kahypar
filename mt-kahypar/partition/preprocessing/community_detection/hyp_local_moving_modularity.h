#pragma once

#include "mt-kahypar/datastructures/community_hypergraph.h"
#include "mt-kahypar/definitions.h"

namespace mt_kahypar::metrics {
Volume hyp_modularity(const ds::CommunityHypergraph& hypergraph);
}

namespace mt_kahypar::community_detection {


struct CommunityMove {
    HypernodeID node_to_move;
    PartitionID destination_community;
    Volume delta;
};

class HypLocalMovingModularity {


public:
    HypLocalMovingModularity(ds::CommunityHypergraph& hypergraph) : _hg(&hypergraph) {
        tbb::parallel_invoke([&] {
            _weight_to_community.resize("Preprocessing", "clearlist_array", hypergraph.initialNumNodes());
        }, [&] {
            community_edge_contribution.resize("Preprocessing", "clearlist_edge_contribution", hypergraph.initialNumNodes(), true, true);
        }, [&] {
            pins_in_community.resize("Preprocessing", "clearlist_pins_in_community", hypergraph.initialNumNodes());
        });
        
        tbb::parallel_invoke([&] {
            _weight_to_community.assign(hypergraph.initialNumNodes(), std::numeric_limits<Volume>::max());
        }, [&] {
            pins_in_community.assign(hypergraph.initialNumNodes(), false);
        });
        
    }

    ~HypLocalMovingModularity() {
        parallel::free(_weight_to_community);
    }

    CommunityMove calculateBestMove(const HypernodeID v) {
        const PartitionID comm_v = _hg->communityID(v);
        // the sum of edgeweights which only have v in that community
        HyperedgeWeight edge_contribution_c = 0;
        // sum of all edgeweights incident to v
        HyperedgeWeight sum_of_edgeweights = 0;
        std::vector<PartitionID> community_neighbours_of_edge;
        std::vector<PartitionID> community_neighbours_of_node;

        for (const HyperedgeID& he : _hg->incidentEdges(v)) {
            std::cout << "Hyperedge nr. " << he << std::endl;
            const HyperedgeWeight edge_weight = _hg->edgeWeight(he);
            for (const HypernodeID& hn : _hg->pins(he)) {
                std::cout << hn << std::endl;
                const PartitionID comm_hn = _hg->communityID(hn);
                if (hn != v && !pins_in_community[comm_hn]) {
                    pins_in_community[comm_hn] = true;
                    community_neighbours_of_edge.emplace_back(comm_hn); //TODO: what happens with comm_v?
                }
            }

            if (!pins_in_community[comm_v]) {
                edge_contribution_c += edge_weight;
            }
            pins_in_community[comm_v] = false;

            for (const PartitionID& community : community_neighbours_of_edge) {
                if (!community_edge_contribution[community]) {
                    community_neighbours_of_node.emplace_back(community);
                }
                community_edge_contribution[community] -= edge_weight;
                pins_in_community[community] = false;
            }
            community_neighbours_of_edge.clear();
            sum_of_edgeweights += edge_weight;
        }

        PartitionID best_community = comm_v;
        Volume best_delta = 0.0;
        const HyperedgeWeight vol_total = _hg->totalVolume();
        const HyperedgeWeight vol_v = _hg->nodeVolume(v);
        const HyperedgeWeight vol_c = _hg->communityVolume(comm_v);
        const HyperedgeWeight total_edge_weight = _hg->totalEdgeWeight();
        for (const PartitionID community : community_neighbours_of_node) {
            community_edge_contribution[community] += sum_of_edgeweights - edge_contribution_c;
            const HyperedgeWeight vol_d = _hg->communityVolume(community);
            Volume exp_edge_contribution = 0.0;
            for (const size_t d : _hg->edgeSizes()) {
                std::cout << "EdgeSize: " << d << std::endl;
                exp_edge_contribution += (static_cast<Volume>(_hg->dEdgeWeight(d)) / powl(vol_total, d)) * (powl(vol_total - vol_c + vol_v, d) - powl(vol_total - vol_c, d)
                                                                                                           + powl(vol_total - vol_d - vol_v, d) - powl(vol_total - vol_d, d));
            }
            Volume delta = (static_cast<Volume>(community_edge_contribution[community]) + exp_edge_contribution) / total_edge_weight;
            if (delta < best_delta) {
                best_delta = delta;
                best_community = community;
            }
            community_edge_contribution[community] = 0.0;
        }
        
        community_neighbours_of_node.clear();
        CommunityMove cm;
        cm.destination_community = best_community;
        cm.delta = best_community == comm_v ? 0.0 : best_delta;
        cm.node_to_move = v;
        return cm;
    }

    // ! executes the given move
    void makeMove(const CommunityMove& move) {
        if (move.delta < 0.0) {
            _hg->addCommunityVolume(move.node_to_move, move.destination_community);
            _hg->subtractCommunityVolume(move.node_to_move, _hg->communityID(move.node_to_move));
            _hg->setCommunityID(move.node_to_move, move.destination_community);
        }
    }

    // ! try to move the given node into a neighbouring community
    // ! returns the change in modularity (0.0 if the node was not moved)
    Volume tryMove(const HypernodeID v) {
        const PartitionID comm_v = _hg->communityID(v);
        // calculate delta for all neighbouring communities
        for (const HyperedgeID& he : _hg->incidentEdges(v)) {
            for (const HypernodeID& hn : _hg->pins(he)) {
                const PartitionID comm_hn = _hg->communityID(hn);
                // if the node hn belonges to a different community which we see for the first time
                if (hn != v && comm_hn != comm_v && !(_weight_to_community[comm_hn] < std::numeric_limits<Volume>::max())) {
                    _weight_to_community[comm_hn] = calculateDelta(v, comm_hn);
                    _neigh_communities.emplace_back(comm_hn);
                }
            }
        }
        // find community with optimal(minimal) delta and reset the clearlist
        Volume best_delta = 0;
        PartitionID best_comm = comm_v;
        for (const PartitionID& p : _neigh_communities) {
            if (_weight_to_community[p] < best_delta) {
                best_delta = _weight_to_community[p];
                best_comm = p;
            }
            _weight_to_community[p] = std::numeric_limits<Volume>::max();
        }
        _neigh_communities.clear();

        // make the actual move only when communities are different
        if (best_comm != comm_v) {
            _hg->setCommunityID(v, best_comm);
            _hg->addCommunityVolume(v, best_comm);
            _hg->subtractCommunityVolume(v, comm_v);
        }
        return best_delta;
    }

    // ! just for testing
    Volume weightTOCommunity(const PartitionID community_id) {
        return _weight_to_community[community_id];
    }



private:

    // ! calculates the change in modularity if the given node would be moved to the given community
    Volume calculateDelta(const HypernodeID v, const PartitionID destination) {
        const PartitionID comm_v = _hg->communityID(v);
        // move to the same community
        if (comm_v == destination) {
            return 0.0;
        }
        // edge contibution
        Volume edge_contribution = 0.0;
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
            // only contribute if the edge has no other pins in the community
            edge_contribution -= pins_in_comm_v ? 0 : _hg->edgeWeight(he);
            edge_contribution += pins_in_destination ? 0 : _hg->edgeWeight(he);
        }

        // expected edge contribution
        const HyperedgeWeight vol_total = _hg->totalVolume();
        const HyperedgeWeight vol_v = _hg->nodeVolume(v);
        const HyperedgeWeight vol_c = _hg->communityVolume(comm_v);
        const HyperedgeWeight vol_d = _hg->communityVolume(destination);
        Volume exp_edge_contribution = 0.0;
        for (const size_t d : _hg->edgeSizes()) {
            exp_edge_contribution += (static_cast<Volume>(_hg->dEdgeWeight(d)) / pow(vol_total, d)) * (powl(vol_total - vol_c + vol_v, d) - powl(vol_total - vol_c, d)
                                                                                                       + powl(vol_total - vol_d - vol_v, d) - powl(vol_total - vol_d, d));
        }
        return (edge_contribution + exp_edge_contribution) / _hg->totalEdgeWeight();
    }

    ds::CommunityHypergraph* _hg;

    // clearlist to find neighbouring community with best modularity value
    ds::Array<Volume> _weight_to_community;
    std::vector<PartitionID> _neigh_communities;

    ds::Array<HyperedgeWeight> community_edge_contribution;
    ds::Array<bool> pins_in_community;
};
}