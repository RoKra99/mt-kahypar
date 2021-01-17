#pragma once

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/community_hypergraph.h"
#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/exponentiations.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/partition/context.h"


namespace mt_kahypar::community_detection {

struct CommunityMove {
    HypernodeID node_to_move;
    PartitionID destination_community;
    Volume delta;
};

class HypergraphLocalMovingModularity {
private:
    using AtomicHyperedgeWeight = parallel::AtomicWrapper<HyperedgeWeight>;
    using CommunityVolumes = parallel::scalable_vector<AtomicHyperedgeWeight>;
    using CommunityVolumeIterator = typename CommunityVolumes::const_iterator;

public:
    static constexpr bool debug = false;
    static constexpr bool enable_heavy_assert = false;

    HypergraphLocalMovingModularity(ds::CommunityHypergraph& hypergraph, const Context& context, const bool deactivate_random = false) :
        _context(context),
        _community_volumes(hypergraph.initialNumNodes()),
        _deactivate_random(deactivate_random) {
        tbb::parallel_invoke([&] {
            _community_edge_contribution.resize("Preprocessing", "clearlist_edge_contribution", hypergraph.initialNumNodes(), true, true);
            }, [&] {
                _powers_of_source_community.resize("Preprocessing", "powers_of_source_community", hypergraph.maxEdgeSize() + 1);
            });

        _reciprocal_vol_total = 1.0L / hypergraph.totalVolume();
    }

    ~HypergraphLocalMovingModularity() {
        parallel::parallel_free(_community_edge_contribution, _powers_of_source_community);
    }

    // ! calculates the best modularity move for the given node
    KAHYPAR_ATTRIBUTE_ALWAYS_INLINE CommunityMove calculateBestMove(const ds::CommunityHypergraph& chg, parallel::scalable_vector<HypernodeID>& communities, const HypernodeID v) {
        ASSERT(_community_neighbours_of_node.empty());
        HEAVY_PREPROCESSING_ASSERT(communityEdgeContributionisEmpty());
        utils::Timer::instance().start_timer("calculate_best_move", "Calculate best move");
        const PartitionID comm_v = communities[v];
        // the sum of edgeweights which only have v in that community
        HyperedgeWeight edge_contribution_c = 0;
        // sum of all edgeweights incident to v
        HyperedgeWeight sum_of_edgeweights = 0;
        utils::Timer::instance().start_timer("edge_contribution", "EdgeContribution");
        for (const HyperedgeID& he : chg.incidentEdges(v)) {
            const HyperedgeWeight edge_weight = chg.edgeWeight(he);
            sum_of_edgeweights += edge_weight;
            for (auto& community : chg._community_counts[he]->singleCuts()) {
                if (community != comm_v && !_community_edge_contribution[community]) {
                    _community_neighbours_of_node.emplace_back(community);
                }
                _community_edge_contribution[community] -= edge_weight;

            }

            for (auto& e : chg._community_counts[he]->multiCut()) {
                if (e.first != comm_v) {
                    if (!_community_edge_contribution[e.first]) {
                        _community_neighbours_of_node.emplace_back(e.first);
                    }
                    _community_edge_contribution[e.first] -= edge_weight;
                } else if (e.second == 1) {
                    _community_edge_contribution[comm_v] -= edge_weight;
                }
            }
        }
        edge_contribution_c = -_community_edge_contribution[comm_v];
        _community_edge_contribution[comm_v] = 0;
        utils::Timer::instance().stop_timer("edge_contribution");



        utils::Timer::instance().start_timer("exp_edge_contribution", "ExpectedEdgeContribution");

        const HyperedgeWeight vol_v = chg.nodeVolume(v);
        const HyperedgeWeight vol_c = _community_volumes[comm_v];
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
            ++overall_checks;
            _community_edge_contribution[community] += sum_of_edgeweights_minus_edgecontribution_c;
            const HyperedgeWeight vol_destination_minus = _community_volumes[community];
            const HyperedgeWeight vol_destination = vol_destination_minus + vol_v;
            const HyperedgeWeight destination_edge_contribution = _community_edge_contribution[community];

            // delta will not be < 0
            if ((destination_edge_contribution >= 0 || best_delta < destination_edge_contribution)
                && vol_c_minus_vol_v <= vol_destination_minus) {
                ++pruned_by_old;
                _community_edge_contribution[community] = 0;
                continue;
            }

            // // pruning via the geometric series
            const Volume destination_fraction = 1.0L - static_cast<Volume>(vol_destination) * _reciprocal_vol_total;
            const Volume destination_fraction_minus = 1.0L - static_cast<Volume>(vol_destination_minus) * _reciprocal_vol_total;

            // precalculate the powers for the source community only once
            // and only if not every possible move is pruned beforehand
            if (!calculated_c) {
                for (const size_t d : chg.edgeSizes()) {
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
                for (const size_t d : chg.edgeSizes()) {
                    const size_t remaining_d = d - biggest_d_yet;
                    power_d_fraction_minus *= math::fast_power(destination_fraction_minus, remaining_d);
                    power_d_fraction *= math::fast_power(destination_fraction, remaining_d);
                    exp_edge_contribution += static_cast<Volume>(chg.edgeWeightBySize(d)) * (_powers_of_source_community[d] + power_d_fraction - power_d_fraction_minus);
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
    KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool makeMove(ds::CommunityHypergraph& chg, parallel::scalable_vector<HypernodeID>& communities, const CommunityMove& move) {
        if (move.delta < 0.0L) {
            const PartitionID source_community = communities[move.node_to_move];
            ASSERT(move.destination_community != source_community);
            _community_volumes[move.destination_community] += chg.nodeVolume(move.node_to_move);
            _community_volumes[source_community] -= chg.nodeVolume(move.node_to_move);
            communities[move.node_to_move] = move.destination_community;
            for (const HyperedgeID& he : chg.incidentEdges(move.node_to_move)) {
                chg._community_counts[he]->addToCommunity(move.destination_community);
                chg._community_counts[he]->removeFromCommunity(source_community);
            }
            return true;
        }
        return false;
    }

    bool localMoving(ds::CommunityHypergraph& chg, parallel::scalable_vector<HypernodeID>& communities) {
        parallel::scalable_vector<HypernodeID> nodes(chg.initialNumNodes());

        //        for (HypernodeID i = 0; i < chg.initialNumNodes(); ++i) {
        tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID i) {
            communities[i] = i;
            nodes[i] = i;
            _community_volumes[i].store(chg.nodeVolume(i));
            });
        //}
        if (!_deactivate_random) {
            utils::Randomize::instance().parallelShuffleVector(nodes, 0UL, nodes.size());
        }
        bool changed_clustering = false;
        size_t nr_nodes_moved = chg.initialNumNodes();
        for (size_t round = 0;
            nr_nodes_moved >= _context.preprocessing.community_detection.min_vertex_move_fraction * chg.initialNumNodes()
            && round < _context.preprocessing.community_detection.max_pass_iterations; ++round) {
            nr_nodes_moved = 0;
            for (HypernodeID& hn : nodes) {
                if (makeMove(chg, communities, calculateBestMove(chg, communities, hn))) {
                    ++nr_nodes_moved;
                }
            }
            changed_clustering |= nr_nodes_moved > 0;
        }
        return changed_clustering;
    }

    void initializeCommunityVolumes(const ds::CommunityHypergraph& chg, const parallel::scalable_vector<HypernodeID>& communities) {
        tbb::parallel_for(0UL, _community_volumes.size(), [&](const size_t i) {
            _community_volumes[i].store(0);
            });
        tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID i) {
            _community_volumes[communities[i]] += chg.nodeVolume(i);
            });
    }

    //TODO: Just for testing
    IteratorRange<CommunityVolumeIterator> communityVolumes(const ds::CommunityHypergraph& chg) const {
        return IteratorRange<CommunityVolumeIterator>(_community_volumes.cbegin(), _community_volumes.cbegin() + chg.initialNumNodes());
    }

    size_t overall_checks = 0;
    size_t pruned_by_old = 0;

private:

    // ! only for testing
    bool communityEdgeContributionisEmpty() {
        bool result = true;
        for (const HyperedgeWeight hw : _community_edge_contribution) {
            result &= hw == 0;
        }
        return result;
    }

    const Context& _context;

    // ! reciprocal total volume
    Volume _reciprocal_vol_total = 0.0L;

    // ! for clearlists
    ds::Array<HyperedgeWeight> _community_edge_contribution;

    // ! contains (vol_V - vol(C)+vol(v))^d - (vol(V)-vol(C))^d for all valid edgesizes d
    // ! Note the values in here are not cleared after each call to calculateBestMove
    ds::Array<Volume> _powers_of_source_community;

    // used in clearlist for calculating the edge contribution
    parallel::scalable_vector<PartitionID> _community_neighbours_of_node;

    parallel::scalable_vector<AtomicHyperedgeWeight> _community_volumes;

    const bool _deactivate_random;

};
}

namespace mt_kahypar::metrics {
Volume hyp_modularity(const ds::CommunityHypergraph& hypergraph, const parallel::scalable_vector<HypernodeID>& communities, const community_detection::HypergraphLocalMovingModularity& hlmm);
}