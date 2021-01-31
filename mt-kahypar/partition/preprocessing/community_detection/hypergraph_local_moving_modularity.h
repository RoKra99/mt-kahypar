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
        _reciprocal_vol_total(1.0L / hypergraph.totalVolume()),
        _community_edge_contribution(hypergraph.initialNumNodes(), 0U),
        _pins_in_community(hypergraph.initialNumNodes(), false),
        _powers_of_source_community(hypergraph.maxEdgeSize() + 1, 0.L),
        _community_volumes(hypergraph.initialNumNodes()),
        _deactivate_random(deactivate_random) {}


    ~HypergraphLocalMovingModularity() = default;

    // ! calculates the best modularity move for the given node
    KAHYPAR_ATTRIBUTE_ALWAYS_INLINE CommunityMove calculateBestMove(ds::CommunityHypergraph& chg, parallel::scalable_vector<HypernodeID>& communities, const HypernodeID v) {
        ASSERT(_community_neighbours_of_node.local().empty());
        HEAVY_PREPROCESSING_ASSERT(communityEdgeContributionisEmpty());
        //utils::Timer::instance().start_timer("calculate_best_move", "Calculate best move");
        const PartitionID comm_v = communities[v];
        ASSERT(static_cast<size_t>(comm_v) < _community_edge_contribution.local().size());
        ASSERT(static_cast<size_t>(comm_v) < _community_edge_contribution.local().size());
        // sum of all edgeweights incident to v
        HyperedgeWeight sum_of_edgeweights = 0;
        //utils::Timer::instance().start_timer("edge_contribution", "EdgeContribution");
        for (const HyperedgeID& he : chg.incidentEdges(v)) {
            const HyperedgeWeight edge_weight = chg.edgeWeight(he);
            sum_of_edgeweights += edge_weight;
            // cuts for this hyperedge are cached
            if (chg.edgeSize(he) > _context.preprocessing.community_detection.hyperedge_size_caching_threshold) {
                for (auto& community : chg.singleCuts(he)) {
                    if (community != comm_v && !_community_edge_contribution.local()[community]) {
                        _community_neighbours_of_node.local().emplace_back(community);
                    }
                    _community_edge_contribution.local()[community] -= edge_weight;

                }

                for (const auto& e : chg.multiCuts(he)) {
                    //LOG << _community_edge_contribution.local().size() << e.first;
                    ASSERT(_community_edge_contribution.local().size() > static_cast<size_t>(e.first), "size: " << _community_edge_contribution.local().size() << ", e.first: " << e.first);
                    if (e.first != comm_v) {
                        if (!_community_edge_contribution.local()[e.first]) {
                            _community_neighbours_of_node.local().emplace_back(e.first);
                        }
                        _community_edge_contribution.local()[e.first] -= edge_weight;
                    } else if (e.second == 1) {
                        _community_edge_contribution.local()[comm_v] -= edge_weight;
                    }
                }
            } else { // hyperedge is not cached
                ASSERT(_community_neighbours_of_edge.local().empty());
                for (const HypernodeID& hn : chg.pins(he)) {
                    const PartitionID comm_hn = communities[hn];
                    if (hn != v && !_pins_in_community.local()[comm_hn]) {
                        _pins_in_community.local()[comm_hn] = true;
                        if (comm_hn != comm_v) {
                            _community_neighbours_of_edge.local().emplace_back(comm_hn);
                        }
                    }
                }

                if (!_pins_in_community.local()[comm_v]) {
                    _community_edge_contribution.local()[comm_v] -= edge_weight;
                }
                _pins_in_community.local()[comm_v] = false;

                for (const PartitionID& community : _community_neighbours_of_edge.local()) {
                    if (!_community_edge_contribution.local()[community]) {
                        _community_neighbours_of_node.local().emplace_back(community);
                    }
                    _community_edge_contribution.local()[community] -= edge_weight;
                    _pins_in_community.local()[community] = false;
                }
                _community_neighbours_of_edge.local().clear();
            }
        }
        const HyperedgeWeight edge_contribution_c = -_community_edge_contribution.local()[comm_v];
        _community_edge_contribution.local()[comm_v] = 0;
        //utils::Timer::instance().stop_timer("edge_contribution");



        //utils::Timer::instance().start_timer("exp_edge_contribution", "ExpectedEdgeContribution");

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
        for (const PartitionID community : _community_neighbours_of_node.local()) {
            ++overall_checks;
            _community_edge_contribution.local()[community] += sum_of_edgeweights_minus_edgecontribution_c;
            const HyperedgeWeight vol_destination_minus = _community_volumes[community];
            const HyperedgeWeight vol_destination = vol_destination_minus + vol_v;
            const HyperedgeWeight destination_edge_contribution = _community_edge_contribution.local()[community];

            // delta will not be < 0
            if ((destination_edge_contribution >= 0 || best_delta < destination_edge_contribution)
                && vol_c_minus_vol_v <= vol_destination_minus) {
                ++pruned_by_old;
                _community_edge_contribution.local()[community] = 0;
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
                    _powers_of_source_community.local()[d] = power_d_fraction_minus - power_d_fraction;
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
                    exp_edge_contribution += static_cast<Volume>(chg.edgeWeightBySize(d)) * (_powers_of_source_community.local()[d] + power_d_fraction - power_d_fraction_minus);
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
            _community_edge_contribution.local()[community] = 0;
        }
        //utils::Timer::instance().stop_timer("exp_edge_contribution");

        _community_neighbours_of_node.local().clear();
        CommunityMove cm;
        cm.destination_community = best_community;
        cm.delta = best_community == comm_v ? 0.0L : best_delta;
        cm.node_to_move = v;
        //utils::Timer::instance().stop_timer("calculate_best_move");
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
                chg.addCommunityToHyperedge(he, move.destination_community);
                chg.removeCommunityFromHyperedge(he, source_community);
            }
            return true;
        }
        return false;
    }

    bool localMoving(ds::CommunityHypergraph& chg, parallel::scalable_vector<HypernodeID>& communities) {
        parallel::scalable_vector<HypernodeID> nodes(chg.initialNumNodes());

        tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID i) {
            communities[i] = i;
            nodes[i] = i;
            _community_volumes[i].store(chg.nodeVolume(i));
            });
        bool changed_clustering = false;
        size_t nr_nodes_moved = chg.initialNumNodes();
        for (size_t round = 0;
            nr_nodes_moved >= _context.preprocessing.community_detection.min_vertex_move_fraction * chg.initialNumNodes()
            && round < _context.preprocessing.community_detection.max_pass_iterations; ++round) {

            if (!_deactivate_random) {
                utils::Randomize::instance().parallelShuffleVector(nodes, 0UL, nodes.size());
            }

            tbb::enumerable_thread_specific<size_t> local_nr_nodes_moved(0);
            tbb::parallel_for(0UL, nodes.size(), [&](const size_t i) {
                const CommunityMove cm = calculateBestMove(chg, communities, nodes[i]);
                if (makeMove(chg, communities, cm)) {
                    ++local_nr_nodes_moved.local();
                }
                });
            nr_nodes_moved = local_nr_nodes_moved.combine(std::plus<>());
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

    //! Iterator for the Community volumes
    IteratorRange<CommunityVolumeIterator> communityVolumes(const ds::CommunityHypergraph& chg) const {
        return IteratorRange<CommunityVolumeIterator>(_community_volumes.cbegin(), _community_volumes.cbegin() + chg.initialNumNodes());
    }

    size_t overall_checks = 0;
    size_t pruned_by_old = 0;

private:

    // ! only for testing
    bool communityEdgeContributionisEmpty() {
        bool result = true;
        for (const HyperedgeWeight hw : _community_edge_contribution.local()) {
            result &= hw == 0;
        }
        return result;
    }

    // ! contains Hyperparameters for the algorithms
    const Context& _context;

    // ! reciprocal total volume
    Volume _reciprocal_vol_total = 0.0L;

    // ! for clearlists
    tbb::enumerable_thread_specific<parallel::scalable_vector<HyperedgeWeight>> _community_edge_contribution;
    tbb::enumerable_thread_specific<parallel::scalable_vector<bool>>_pins_in_community;

    // ! contains (vol_V - vol(C)+vol(v))^d - (vol(V)-vol(C))^d for all valid edgesizes d
    // ! Note the values in here are not cleared after each call to calculateBestMove
    tbb::enumerable_thread_specific<parallel::scalable_vector<Volume>> _powers_of_source_community;

    // ! used in clearlist for calculating the edge contribution
    tbb::enumerable_thread_specific<parallel::scalable_vector<PartitionID>> _community_neighbours_of_node;

    // used in clearlist for calculating the edge contribution
    tbb::enumerable_thread_specific<parallel::scalable_vector<PartitionID>> _community_neighbours_of_edge;

    // ! volumes of each community
    parallel::scalable_vector<AtomicHyperedgeWeight> _community_volumes;

    // ! deactivates random node order in local moving
    const bool _deactivate_random;

};
}

namespace mt_kahypar::metrics {
Volume hyp_modularity(const ds::CommunityHypergraph& hypergraph, const parallel::scalable_vector<HypernodeID>& communities, const community_detection::HypergraphLocalMovingModularity& hlmm);
}