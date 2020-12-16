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
                _powers_of_source_community.resize("Preprocessing", "powers_of_source_community", hypergraph.maxEdgeSize() + 1);
            });

        _pins_in_community.assign(hypergraph.initialNumNodes(), false);

        _reciprocal_vol_total = 1.0L / _hg->totalVolume();
    }

    ~HypergraphLocalMovingModularity() {
        parallel::parallel_free(_community_edge_contribution, _pins_in_community, _powers_of_source_community);
    }

    // ! calculates the best modularity move for the given node
    CommunityMove calculateBestMove(const HypernodeID v) {
        ASSERT(_community_neighbours_of_edge.empty());
        ASSERT(_community_neighbours_of_node.empty());
        ASSERT(communityEdgeContributionisEmpty());
        utils::Timer::instance().start_timer("calculate_best_move", "Calculate best move");
        const PartitionID comm_v = _hg->communityID(v);
        // the sum of edgeweights which only have v in that community
        HyperedgeWeight edge_contribution_c = 0;
        // sum of all edgeweights incident to v
        HyperedgeWeight sum_of_edgeweights = 0;
        utils::Timer::instance().start_timer("edge_contribution", "EdgeContribution");
        for (const HyperedgeID& he : _hg->incidentEdges(v)) {
            const HyperedgeWeight edge_weight = _hg->edgeWeight(he);
            ASSERT(_community_neighbours_of_edge.empty());
            for (const HypernodeID& hn : _hg->pins(he)) {
                const PartitionID comm_hn = _hg->communityID(hn);
                if (hn != v && !_pins_in_community[comm_hn]) {
                    _pins_in_community[comm_hn] = true;
                    _community_neighbours_of_edge.emplace_back(comm_hn);
                }
            }

            if (!_pins_in_community[comm_v]) {
                edge_contribution_c += edge_weight;
            }
            _pins_in_community[comm_v] = false;

            for (const PartitionID& community : _community_neighbours_of_edge) {
                if (!_community_edge_contribution[community]) {
                    _community_neighbours_of_node.emplace_back(community);
                }
                _community_edge_contribution[community] -= edge_weight;
                _pins_in_community[community] = false;
            }
            _community_neighbours_of_edge.clear();
            sum_of_edgeweights += edge_weight;
        }
        utils::Timer::instance().stop_timer("edge_contribution");



        utils::Timer::instance().start_timer("exp_edge_contribution", "ExpectedEdgeContribution");

        const HyperedgeWeight vol_v = _hg->nodeVolume(v);
        const HyperedgeWeight vol_c = _hg->communityVolume(comm_v);
        const HyperedgeWeight vol_c_minus_vol_v = vol_c - vol_v;

        const Volume source_fraction_minus = 1.0L - static_cast<Volume>(vol_c_minus_vol_v) * _reciprocal_vol_total;
        const Volume source_fraction = 1.0L - static_cast<Volume>(vol_c) * _reciprocal_vol_total;
        size_t biggest_d_yet = 1;
        Volume power_d_fraction_minus = source_fraction_minus;
        Volume power_d_fraction = source_fraction;
        bool calculated_c = false;

        PartitionID best_community = comm_v;
        Volume best_delta = 0.0L;
        const HyperedgeWeight sum_of_edgeweights_minus_edgecontribution_c = sum_of_edgeweights - edge_contribution_c;

        // expected edgecontribution starts here
        for (const PartitionID community : _community_neighbours_of_node) {
            //++overall_checks;
            _community_edge_contribution[community] += sum_of_edgeweights_minus_edgecontribution_c;
            const HyperedgeWeight vol_destination_minus = _hg->communityVolume(community);
            const HyperedgeWeight vol_destination = vol_destination_minus + vol_v;
            const HyperedgeWeight destination_edge_contribution = _community_edge_contribution[community];

            // delta will not be < 0
            if ((destination_edge_contribution >= 0 || best_delta < destination_edge_contribution)
                && vol_c_minus_vol_v <= vol_destination_minus) {
                //++pruned_by_old;
                _community_edge_contribution[community] = 0;
                continue;
            }

            // // pruning via the geometric series
            const Volume destination_fraction = 1.0L - static_cast<Volume>(vol_destination) * _reciprocal_vol_total;
            const Volume destination_fraction_minus = 1.0L - static_cast<Volume>(vol_destination_minus) * _reciprocal_vol_total;

            // precalculate the powers for the source community only once
            // and only if not every possible move is pruned beforehand
            if (!calculated_c) {
                for (const size_t d : _hg->edgeSizes()) {
                    const size_t remaining_d = d - biggest_d_yet;
                    power_d_fraction_minus *= math::fast_power(source_fraction_minus, remaining_d);
                    power_d_fraction *= math::fast_power(source_fraction, remaining_d);
                    _powers_of_source_community[d] = power_d_fraction_minus - power_d_fraction;
                    biggest_d_yet = d;
                }
                calculated_c = true;
            }

            Volume exp_edge_contribution = 0.0L;
            // if this is equal the expected_edge_contribution will be 0
            if (vol_c_minus_vol_v != vol_destination_minus) {
                biggest_d_yet = 1;
                power_d_fraction_minus = destination_fraction_minus;
                power_d_fraction = destination_fraction;
                //actual calculation of the expected edge contribution for the given community
                for (const size_t d : _hg->edgeSizes()) {
                    const size_t remaining_d = d - biggest_d_yet;
                    power_d_fraction_minus *= math::fast_power(destination_fraction_minus, remaining_d);
                    power_d_fraction *= math::fast_power(destination_fraction, remaining_d);
                    exp_edge_contribution += static_cast<Volume>(_hg->edgeWeightBySize(d)) * (_powers_of_source_community[d] + power_d_fraction - power_d_fraction_minus);
                    biggest_d_yet = d;
                }
                ASSERT((vol_c_minus_vol_v > vol_destination_minus && exp_edge_contribution < 0.0L)
                    || (vol_c_minus_vol_v < vol_destination_minus&& exp_edge_contribution > 0.0L)
                    || (vol_c_minus_vol_v == vol_destination_minus));
            }

            Volume delta = static_cast<Volume>(destination_edge_contribution) + exp_edge_contribution;
            if (delta < best_delta) {
                best_delta = delta;
                best_community = community;
            }
            _community_edge_contribution[community] = 0;
        }
        utils::Timer::instance().stop_timer("exp_edge_contribution");

        _community_neighbours_of_node.clear();
        CommunityMove cm;
        cm.destination_community = best_community;
        cm.delta = best_community == comm_v ? 0.0L : best_delta;
        cm.node_to_move = v;
        utils::Timer::instance().stop_timer("calculate_best_move");
        return cm;
    }

    // ! executes the given move
    bool makeMove(const CommunityMove& move) {
        if (move.delta < 0.0L) {
            ASSERT(move.destination_community != _hg->communityID(move.node_to_move));
            _hg->addCommunityVolume(move.node_to_move, move.destination_community);
            _hg->subtractCommunityVolume(move.node_to_move, _hg->communityID(move.node_to_move));
            _hg->setCommunityID(move.node_to_move, move.destination_community);
            return true;
        }
        return false;
    }

    //size_t overall_checks = 0;
    //size_t pruned_by_old = 0;

private:

    // ! only for testing
    bool communityEdgeContributionisEmpty() {
        bool result = true;
        for (const HyperedgeWeight hw : _community_edge_contribution) {
            result &= hw == 0;
        }
        return result;
    }

    // ! reciprocal total volume
    Volume _reciprocal_vol_total = 0.0L;

    // ! Hypergraph on which the local moving is performed
    ds::CommunityHypergraph* _hg;

    // ! for clearlists
    ds::Array<HyperedgeWeight> _community_edge_contribution;
    ds::Array<bool> _pins_in_community;

    // ! contains (vol_V - vol(C)+vol(v))^d - (vol(V)-vol(C))^d for all valid edgesizes d
    // ! Note the values in here are not cleared after each call to calculateBestMove
    ds::Array<Volume> _powers_of_source_community;

    // used in clearlist for calculating the edge contribution
    std::vector<PartitionID> _community_neighbours_of_edge;
    std::vector<PartitionID> _community_neighbours_of_node;

};
}