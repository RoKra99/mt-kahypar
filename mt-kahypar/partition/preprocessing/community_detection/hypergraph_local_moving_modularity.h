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
        utils::Timer::instance().start_timer("calculate_best_move", "Calculate best move");
        static constexpr bool debug = false;
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
        //TODO: Note: The smaller the volume of D the "better"(smaller) the expected edge contribution 
        // utils::Timer::instance().start_timer("Pruning EdgeContribution", "pruning");
        // key = edge contribution, value = index in community_neighbours_of_node
        // std::unordered_map<HyperedgeWeight, PartitionID> best_volume_for_edge_contribution;
        // for (size_t i = 0; i < community_neighbours_of_node.size(); ++i) {
        //     const PartitionID communityID = community_neighbours_of_node[i];
        //     auto iter = best_volume_for_edge_contribution.emplace(_community_edge_contribution[communityID], i);

        //     // key was already in map
        //     if (!iter.second) {
        //         // old value is better
        //         if (_hg->communityVolume(community_neighbours_of_node[iter.first->second]) < _hg->communityVolume(communityID)) {
        //             community_neighbours_of_node[i] = -1;
        //         }
        //         else { // new value is better
        //             community_neighbours_of_node[iter.first->second] = -1;
        //             iter.first->second = i;
        //         }
        //     }
        // }


        // std::vector<size_t> sorted = sort_indexes(community_neighbours_of_node);
        // HyperedgeWeight current_edge_contribution = std::numeric_limits<HyperedgeWeight>::max();
        // PartitionID current_best_index = -1;
        // HyperedgeWeight current_best_volume = -1;
        // for (size_t i : sorted) {
        //     const PartitionID communityID = community_neighbours_of_node[i];
        //     const HyperedgeWeight edge_contribution = _community_edge_contribution[communityID];
        //     const HyperedgeWeight community_volume = _hg->communityVolume(communityID);

        //     if (current_edge_contribution == edge_contribution) {
        //         if (current_best_volume < community_volume) {
        //             community_neighbours_of_node[i] = -1;
        //         } else {
        //             community_neighbours_of_node[current_best_index] = -1;
        //             current_best_index = i;
        //             current_best_volume = community_volume;
        //         }
        //     } else {
        //         current_edge_contribution = edge_contribution;
        //         current_best_index = i;
        //         current_best_volume = community_volume;
        //     }
        // }



        //utils::Timer::instance().stop_timer("pruning");



        utils::Timer::instance().start_timer("exp_edge_contribution", "ExpectedEdgeContribution");
        PartitionID best_community = comm_v;
        Volume best_delta = 0.0;
        const HyperedgeWeight vol_total = _hg->totalVolume();
        const HyperedgeWeight vol_v = _hg->nodeVolume(v);
        const HyperedgeWeight vol_c = _hg->communityVolume(comm_v);

        //precalculate the powers for the source community
        const Volume source_fraction_minus = 1.0L - static_cast<Volume>(vol_c - vol_v) / vol_total;
        const Volume source_fraction = 1.0L - static_cast<Volume>(vol_c) / vol_total;
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

        // actual calculation for the expected edge contribution of each neighbour community
        for (const PartitionID community : community_neighbours_of_node) {
            //if (community == -1) continue;
            _community_edge_contribution[community] += sum_of_edgeweights_minus_edgecontribution_c;
            const HyperedgeWeight vol_destination = _hg->communityVolume(community);

            // delta will not be < 0
            if ((_community_edge_contribution[community] >= 0 || best_delta < _community_edge_contribution[community]) && vol_c_minus_vol_v <= vol_destination) continue;

            const Volume destination_fraction = 1.0L - static_cast<Volume>(vol_destination + vol_v) / vol_total;
            const Volume destination_fraction_minus = 1.0L - static_cast<Volume>(vol_destination) / vol_total;
            Volume exp_edge_contribution = 0.0L;

            // if this is equal the expected_edge_contribution will be 0
            if (vol_c_minus_vol_v != vol_destination) {
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
            }

            Volume delta = static_cast<Volume>(_community_edge_contribution[community]) + exp_edge_contribution;
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

    std::vector<size_t> sort_indexes(const std::vector<PartitionID>& v) {
        std::vector<size_t> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);

        std::stable_sort(idx.begin(), idx.end(), [&v, this](size_t i1, size_t i2) {
            return _community_edge_contribution[v[i1]] < _community_edge_contribution[v[i2]];
                         });
        return idx;
    }

    // ! Hypergraph on which the local moving is performed
    ds::CommunityHypergraph* _hg;

    // ! for clearlists
    ds::Array<HyperedgeWeight> _community_edge_contribution;
    ds::Array<bool> _pins_in_community;


    // ! contains (vol_V - vol(C)+vol(v))^d - (vol(V)-vol(C))^d for all valid edgesizes d
    // TODO: Note the values in here are not cleared after each call to tryMove
    ds::Array<Volume> _powers_of_source_community;
};
}