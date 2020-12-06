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
    }

    ~HypergraphLocalMovingModularity() {
        parallel::parallel_free(_community_edge_contribution, _pins_in_community, _powers_of_source_community);
    }

    // ! calculates the best modularity move for the given node
    CommunityMove calculateBestMove(const HypernodeID v) {
        ASSERT(communityEdgeContributionisEmpty());
        utils::Timer::instance().start_timer("calculate_best_move", "Calculate best move");
        const PartitionID comm_v = _hg->communityID(v);
        // the sum of edgeweights which only have v in that community
        HyperedgeWeight edge_contribution_c = 0;
        // sum of all edgeweights incident to v
        HyperedgeWeight sum_of_edgeweights = 0;
        std::vector<PartitionID> community_neighbours_of_edge;
        // -1 if invalid
        std::vector<PartitionID> community_neighbours_of_node;
        utils::Timer::instance().start_timer("edge_contribution", "EdgeContribution");
        for (const HyperedgeID& he : _hg->incidentEdges(v)) {
            const HyperedgeWeight edge_weight = _hg->edgeWeight(he);
            ASSERT(community_neighbours_of_edge.empty());
            for (const HypernodeID& hn : _hg->pins(he)) {
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
        const HyperedgeWeight vol_total = _hg->totalVolume();
        const HyperedgeWeight vol_v = _hg->nodeVolume(v);
        const HyperedgeWeight vol_c = _hg->communityVolume(comm_v);
        const Volume reciprocal_vol_total = 1.0L / vol_total;

        //precalculate the powers for the source community
        const Volume source_fraction_minus = 1.0L - static_cast<Volume>(vol_c - vol_v) * reciprocal_vol_total;
        const Volume source_fraction = 1.0L - static_cast<Volume>(vol_c) * reciprocal_vol_total;
        size_t biggest_d_yet = 1;
        Volume d_minus_prev = source_fraction_minus;
        Volume d_prev = source_fraction;

        for (const size_t d : _hg->edgeSizes()) {
            const size_t remaining_d = d - biggest_d_yet;
            d_minus_prev *= math::fast_power(source_fraction_minus, remaining_d);
            d_prev *= math::fast_power(source_fraction, remaining_d);
            _powers_of_source_community[d] = d_minus_prev - d_prev;
            biggest_d_yet = d;
        }

        const HyperedgeWeight sum_of_edgeweights_minus_edgecontribution_c = sum_of_edgeweights - edge_contribution_c;
        const HyperedgeWeight vol_c_minus_vol_v = vol_c - vol_v;
        HyperedgeWeight max_edgeWeight = _hg->maxAccumulatedEdgeWeight();
        HyperedgeWeight nr_unique_edgesizes = _hg->numberOfUniqueEdgeSizes();
        Volume geometric_approximation_c_minus = 0.0L;
        Volume geometric_approximation_c = 0.0L;
        // for pruning with the geometric series
        if (vol_c_minus_vol_v) {
            geometric_approximation_c_minus = geometric_approximation(vol_c_minus_vol_v, source_fraction_minus, nr_unique_edgesizes);
            geometric_approximation_c = geometric_approximation(vol_c, source_fraction, nr_unique_edgesizes);
        }

        PartitionID best_community = comm_v;
        Volume best_delta = 0.0;

        // actual calculation for the expected edge contribution of each neighbour community
        for (const PartitionID community : community_neighbours_of_node) {
            _community_edge_contribution[community] += sum_of_edgeweights_minus_edgecontribution_c;
            const HyperedgeWeight vol_destination_minus = _hg->communityVolume(community);
            const HyperedgeWeight vol_destination = vol_destination_minus + vol_v;

            // delta will not be < 0
            if ((_community_edge_contribution[community] >= 0 || best_delta < _community_edge_contribution[community])
                && vol_c_minus_vol_v <= vol_destination) {
                _community_edge_contribution[community] = 0;
                continue;
            }

            const Volume destination_fraction = 1.0L - static_cast<Volume>(vol_destination) * reciprocal_vol_total;
            const Volume destination_fraction_minus = 1.0L - static_cast<Volume>(vol_destination_minus) * reciprocal_vol_total;
            Volume geometric_approximation_exp = std::numeric_limits<Volume>::max();
            if (vol_c_minus_vol_v > vol_destination_minus) {
                Volume geometric_approximation_d = geometric_approximation(vol_destination, destination_fraction, nr_unique_edgesizes);
                Volume geometric_approximation_d_minus = geometric_approximation(vol_destination_minus, destination_fraction_minus, nr_unique_edgesizes);

                // trying to avoid overflow by changing the multiplication order
                geometric_approximation_exp = vol_total
                    * (geometric_approximation_c_minus - geometric_approximation_c + geometric_approximation_d - geometric_approximation_d_minus)
                    * max_edgeWeight;

                if (_community_edge_contribution[community] + geometric_approximation_exp > 0.0L
                    || _community_edge_contribution[community] + geometric_approximation_exp > best_delta) {
                    //TESTING
                    //++pruned;
                    _community_edge_contribution[community] = 0;
                    continue;
                }
            }

            Volume exp_edge_contribution = 0.0L;
            // if this is equal the expected_edge_contribution will be 0
            if (vol_c_minus_vol_v != vol_destination_minus) {
                biggest_d_yet = 1;
                d_minus_prev = destination_fraction_minus;
                d_prev = destination_fraction;
                for (const size_t d : _hg->edgeSizes()) {
                    const size_t remaining_d = d - biggest_d_yet;
                    d_minus_prev *= math::fast_power(destination_fraction_minus, remaining_d);
                    d_prev *= math::fast_power(destination_fraction, remaining_d);
                    exp_edge_contribution += static_cast<Volume>(_hg->edgeWeightBySize(d)) * (_powers_of_source_community[d] + d_prev - d_minus_prev);
                    biggest_d_yet = d;
                }
                ASSERT((vol_c_minus_vol_v > vol_destination_minus && exp_edge_contribution < 0.0L)
                    || (vol_c_minus_vol_v < vol_destination_minus&& exp_edge_contribution > 0.0L)
                    || (vol_c_minus_vol_v == vol_destination_minus));
            }

            // if (geometric_approximation_exp < std::numeric_limits<Volume>::max()) {
            //     Volume ratio = geometric_approximation_exp / exp_edge_contribution;
            //     ratios.emplace_back(ratio);
            // }

            Volume delta = static_cast<Volume>(_community_edge_contribution[community]) + exp_edge_contribution;
            if (delta < best_delta) {
                best_delta = delta;
                best_community = community;
            }
            _community_edge_contribution[community] = 0;
        }
        utils::Timer::instance().stop_timer("exp_edge_contribution");

        community_neighbours_of_node.clear();
        CommunityMove cm;
        cm.destination_community = best_community;
        cm.delta = best_community == comm_v ? 0.0L : best_delta;
        cm.node_to_move = v;
        utils::Timer::instance().stop_timer("calculate_best_move");
        ASSERT(communityEdgeContributionisEmpty());
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

    //TESTING
    // std::vector<Volume> ratios;
    // size_t pruned = 0;
private:

    // only for testing
    bool communityEdgeContributionisEmpty() {
        bool result = true;
        for (const HyperedgeWeight hw : _community_edge_contribution) {
            result &= hw == 0;
        }
        return result;
    }

    // TODO: currently using 2 as d_min (can we assume that 2 will always be the smallest edge size?)
    Volume geometric_approximation(HyperedgeWeight vol_community, Volume q, size_t d_max) const {
        ASSERT(vol_community > 0);
        return math::fast_power(q, 2) * (1.0L - math::fast_power(q, d_max)) / vol_community;
    }

    // ! Hypergraph on which the local moving is performed
    ds::CommunityHypergraph* _hg;

    // ! for clearlists
    ds::Array<HyperedgeWeight> _community_edge_contribution;
    ds::Array<bool> _pins_in_community;

    // ! contains (vol_V - vol(C)+vol(v))^d - (vol(V)-vol(C))^d for all valid edgesizes d
    // TODO: Note the values in here are not cleared after each call to calculateBestMove
    ds::Array<Volume> _powers_of_source_community;

};
}