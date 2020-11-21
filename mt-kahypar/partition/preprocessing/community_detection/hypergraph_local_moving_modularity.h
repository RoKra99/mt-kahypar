#pragma once

#include "mt-kahypar/datastructures/community_hypergraph.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/exponentiations.h"
#include "mt-kahypar/utils/timer.h"

namespace mt_kahypar::metrics {
Volume hyp_modularity(const ds::CommunityHypergraph& hypergraph);
}

namespace mt_kahypar::community_detection {


struct CommunityMove {
    HypernodeID node_to_move;
    PartitionID destination_community;
    Volume delta;
};

class HypergraphLocalMovingModularity {


public:
    HypergraphLocalMovingModularity(ds::CommunityHypergraph& hypergraph) : _hg(&hypergraph) {
        tbb::parallel_invoke([&] {
            _community_edge_contribution.resize("Preprocessing", "clearlist_edge_contribution", hypergraph.initialNumNodes(), true, true);
                             }, [&] {
                                 _pins_in_community.resize("Preprocessing", "clearlist_pins_in_community", hypergraph.initialNumNodes());
                             }, [&] {
                                 _powers_of_total_volume.resize("Preprocessing", "powers_of_total_volume", hypergraph.maxEdgeSize() + 1);
                             }, [&] {
                                 _powers_of_source_community.resize("Preprocessing", "powers_of_source_community", hypergraph.maxEdgeSize() + 1);
                             });

        _pins_in_community.assign(hypergraph.initialNumNodes(), false);
        computeAndSetPowers();
    }

    ~HypergraphLocalMovingModularity() {
        parallel::parallel_free(_community_edge_contribution, _pins_in_community, _powers_of_total_volume, _powers_of_source_community);
    }

    // ! calculates the best modularity move for the given node
    CommunityMove calculateBestMove(const HypernodeID v) {
        static constexpr bool debug = false;
        const PartitionID comm_v = _hg->communityID(v);
        // the sum of edgeweights which only have v in that community
        HyperedgeWeight edge_contribution_c = 0;
        // sum of all edgeweights incident to v
        HyperedgeWeight sum_of_edgeweights = 0;
        std::vector<PartitionID> community_neighbours_of_edge;
        std::vector<PartitionID> community_neighbours_of_node;
        utils::Timer::instance().start_timer("edge_contribution", "EdgeContribution");
        for (const HyperedgeID& he : _hg->incidentEdges(v)) {
            DBG << "Hyperedge nr. " << he;
            const HyperedgeWeight edge_weight = _hg->edgeWeight(he);
            for (const HypernodeID& hn : _hg->pins(he)) {
                DBG << "Pin in edge: " << hn;
                const PartitionID comm_hn = _hg->communityID(hn);
                if (hn != v && !_pins_in_community[comm_hn]) {
                    _pins_in_community[comm_hn] = true;
                    community_neighbours_of_edge.emplace_back(comm_hn); //TODO: what happens with comm_v?
                }
            }

            if (!_pins_in_community[comm_v]) {
                edge_contribution_c += edge_weight;
            }
            _pins_in_community[comm_v] = false;

            for (const PartitionID& community : community_neighbours_of_edge) {
                if (!_community_edge_contribution[community]) {
                    community_neighbours_of_node.emplace_back(community);
                }
                _community_edge_contribution[community] -= edge_weight;
                _pins_in_community[community] = false;
            }
            community_neighbours_of_edge.clear();
            sum_of_edgeweights += edge_weight;
        }
        utils::Timer::instance().stop_timer("edge_contribution");

        PartitionID best_community = comm_v;
        Volume best_delta = 0.0;
        const HyperedgeWeight vol_total = _hg->totalVolume();
        const HyperedgeWeight vol_v = _hg->nodeVolume(v);
        const HyperedgeWeight vol_c = _hg->communityVolume(comm_v);

        //precalculate the powers for the source community
        for (const size_t d : _hg->edgeSizes()) { 
            _powers_of_source_community[d] = math::fast_power(vol_total - vol_c + vol_v, d) - math::fast_power(vol_total - vol_c, d);
        }

        utils::Timer::instance().start_timer("exp_edge_contribution", "ExpectedEdgeContribution");
        for (const PartitionID community : community_neighbours_of_node) {
            _community_edge_contribution[community] += sum_of_edgeweights - edge_contribution_c;
            const HyperedgeWeight vol_d = _hg->communityVolume(community);
            Volume exp_edge_contribution = 0.0L;
            // only the calculations for vol_d change here
            for (const size_t d : _hg->edgeSizes()) {
                exp_edge_contribution += _powers_of_total_volume[d] * (_powers_of_source_community[d] + math::fast_power(vol_total - vol_d - vol_v, d) - math::fast_power(vol_total - vol_d, d));
            }
            Volume delta = (static_cast<Volume>(_community_edge_contribution[community]) + exp_edge_contribution);
            if (delta < best_delta) {
                best_delta = delta;
                best_community = community;
            }
            _community_edge_contribution[community] = 0.0L;
        }
        utils::Timer::instance().stop_timer("exp_edge_contribution");

        community_neighbours_of_node.clear();
        CommunityMove cm;
        cm.destination_community = best_community;
        cm.delta = best_community == comm_v ? 0.0L : best_delta;
        cm.node_to_move = v;
        return cm;
    }

    // ! executes the given move
    void makeMove(const CommunityMove& move) {
        static constexpr bool debug = false;
        DBG << "Move from " << _hg->communityID(move.node_to_move) << " to " << move.destination_community;
        if (move.delta < 0.0L) {
            _hg->addCommunityVolume(move.node_to_move, move.destination_community);
            _hg->subtractCommunityVolume(move.node_to_move, _hg->communityID(move.node_to_move));
            _hg->setCommunityID(move.node_to_move, move.destination_community);
        }
    }

private:

    void computeAndSetPowers() {
        const HyperedgeWeight vol_total = _hg->totalVolume();
        for (const size_t d : _hg->edgeSizes()) {
            _powers_of_total_volume[d] = static_cast<Volume>(_hg->edgeWeightBySize(d)) / math::fast_power(vol_total, d);
        }
    }

    // ! Hypergraph on which the local moving is performed
    ds::CommunityHypergraph* _hg;

    // ! for clearlists
    ds::Array<HyperedgeWeight> _community_edge_contribution;
    ds::Array<bool> _pins_in_community;

    // ! contains the factor |E_d| / vol(V)^d for all valid edgesizes d
    ds::Array<Volume> _powers_of_total_volume;

    // ! contains (vol_V - vol(C)+vol(v))^d - (vol(V)-vol(C))^d for all valid edgesizes d
    // TODO: Note the values in here are not cleared after each call to tryMove
    ds::Array<HyperedgeWeight> _powers_of_source_community;
};
}