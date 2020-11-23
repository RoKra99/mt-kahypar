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
        utils::Timer::instance().start_timer("calculate_best_move", "Calculate best move");
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
        utils::Timer::instance().start_timer("exp_edge_contribution", "ExpectedEdgeContribution");
        PartitionID best_community = comm_v;
        Volume best_delta = 0.0;
        const HyperedgeWeight vol_total = _hg->totalVolume();
        const HyperedgeWeight vol_v = _hg->nodeVolume(v);
        const HyperedgeWeight vol_c = _hg->communityVolume(comm_v);

        //precalculate the powers for the source community
        const Volume source_fraction_minus = 1.0L - static_cast<Volume>(vol_c - vol_v)/vol_total;
        const Volume source_fraction = 1.0L - static_cast<Volume>(vol_c) / vol_total;
        for (const size_t d : _hg->edgeSizes()) {
            _powers_of_source_community[d] = math::fast_power(source_fraction_minus, d) - math::fast_power(source_fraction, d);
        }
        // std::vector<Volume> destination_powers(_hg->maxEdgeSize() + 1, 0.0L);
        // std::vector<Volume> destination_powers_minus(_hg->maxEdgeSize() + 1, 0.0L);
        // destination_powers[0] = 1.0L;
        // destination_powers_minus[0] = 1.0L;
        // actual calculation for the expected edge contribution of each neighbor community
        for (const PartitionID community : community_neighbours_of_node) {
            _community_edge_contribution[community] += sum_of_edgeweights - edge_contribution_c;
            const HyperedgeWeight vol_destination = _hg->communityVolume(community);
            const Volume destination_fraction = 1.0L - static_cast<Volume>(vol_destination + vol_v) / vol_total;
            const Volume destination_fraction_minus = 1.0L - static_cast<Volume>(vol_destination) / vol_total;
            Volume exp_edge_contribution = 0.0L;
            
            size_t biggest_d_yet = 1;
            Volume d_minus_prev = destination_fraction_minus;
            Volume d_prev = destination_fraction;
            for (const size_t d : _hg->edgeSizes()) {
                const size_t remaining_d = d - biggest_d_yet;
                d_minus_prev *= math::fast_power(destination_fraction_minus, remaining_d);
                d_prev *= math::fast_power(destination_fraction, remaining_d);
                exp_edge_contribution += static_cast<Volume>(_hg->edgeWeightBySize(d)) * (_powers_of_source_community[d] + d_prev - d_minus_prev);
                biggest_d_yet = d;
            }
            // for (size_t i = 1; i < _hg->maxEdgeSize() + 1; ++i) {
            //     destination_powers[i] = 0.0L;
            //     destination_powers_minus[i] = 0.0L;
            // }
            
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
        utils::Timer::instance().stop_timer("calculate_best_move");
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
    ds::Array<Volume> _powers_of_source_community;
};
}