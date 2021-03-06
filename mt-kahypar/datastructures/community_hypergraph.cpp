#include "community_hypergraph.h"

#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/datastructures/concurrent_bucket_map.h"
#include "mt-kahypar/utils/timer.h"

#include <tbb/parallel_reduce.h>

namespace mt_kahypar::ds {

//TODO:  "single pin edges" save them for dbg purposes
CommunityHypergraph CommunityHypergraph::contract(StaticHypergraph& hypergraph, parallel::scalable_vector<HypernodeID>& communities) {
    utils::Timer::instance().start_timer("community_hypergraph_contract", "CommunityHypergaph Contraction");
    hypergraph = _hg->contract(communities, 0, false);
    utils::Timer::instance().start_timer("community_specific", "Community Specific");
    CommunityHypergraph chg(_context);
    chg._hg = &hypergraph;
    chg._vol_v = _vol_v;
    chg._valid_edge_sizes = std::move(_valid_edge_sizes);
    chg._d_edge_weights = std::move(_d_edge_weights);
    chg._total_edge_weight = _total_edge_weight;
    chg._community_count_locks = std::move(_community_count_locks);

    Array<parallel::AtomicWrapper<HyperedgeWeight>>& tmp_node_volumes = _tmp_community_hypergraph_buffer->tmp_node_volumes;
    Array<HypernodeID>& multipin_mapping = _tmp_community_hypergraph_buffer->multipin_mapping;
    const Array<HypernodeID>& incidence_array = hypergraph._incidence_array;


    tbb::parallel_for(0U, initialNumNodes(), [&](const HypernodeID i) {
        tmp_node_volumes[i].store(0);
        });
    // map according to hypergraph mapping
    tbb::parallel_for(0U, initialNumNodes(), [&](const HypernodeID i) {
        const HypernodeID coarse_i = communities[i];
        tmp_node_volumes[coarse_i] += nodeVolume(i);
        });

    chg._node_volumes.resize(chg.initialNumNodes());
    // actually writing the volumes to the graph
    tbb::parallel_for(0U, chg.initialNumNodes(), [&](const HypernodeID i) {
        chg._node_volumes[i] = std::move(tmp_node_volumes[i]);
        });

    chg._community_counts.resize(chg.initialNumEdges());
    tbb::parallel_for(0U, chg.initialNumEdges(), [&](const HyperedgeID he) {
        chg._community_counts[he] = chg.edgeSize(he) > _context.preprocessing.community_detection.hyperedge_size_caching_threshold
            ? std::make_unique<CommunityCount<Map>>(chg.edgeSize(he), chg.pins(he), true) : std::unique_ptr<CommunityCount<Map>>(nullptr);
        });


    const size_t num_pins = hypergraph.initialNumPins();
    const size_t num_hyperedges = hypergraph.initialNumEdges();
    multipin_mapping.assign(num_pins, 0U);
    tbb::parallel_for(0UL, num_hyperedges, [&](const HyperedgeID he) {
        HypernodeID previous = std::numeric_limits<HypernodeID>::max();
        const size_t incidence_array_start = hypergraph.hyperedge(he).firstEntry();
        const size_t incidence_array_end = hypergraph.hyperedge(he).firstInvalidEntry();
        for (size_t i = incidence_array_start; i < incidence_array_end; ++i) {
            if (previous != incidence_array[i]) {
                previous = incidence_array[i];
                multipin_mapping[i] = 1;
            }
        }
        });

    parallel::TBBPrefixSum<HypernodeID, Array> multipin_mapping_prefix_sum(multipin_mapping);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, num_pins), multipin_mapping_prefix_sum);
    tbb::parallel_for(0UL, num_hyperedges, [&](const HyperedgeID he) {
        const size_t incidence_array_start = hypergraph.hyperedge(he).firstEntry();
        const size_t incidence_array_end = hypergraph.hyperedge(he).firstInvalidEntry();
        HypernodeID previous = incidence_array[incidence_array_start];
        size_t write_pos = multipin_mapping_prefix_sum[incidence_array_start];
        //ASSERT(he < _multipin_indexes.size() && write_pos < _multipin_incidence_array.size());
        _multipin_indexes[he] = write_pos;
        HypernodeID count = 0;
        for (size_t i = incidence_array_start; i < incidence_array_end; ++i) {
            if (previous != incidence_array[i]) {
                //ASSERT(_multipin_incidence_array.begin() + write_pos < _multipin_incidence_array.end());
                ASSERT(previous < hypergraph.initialNumNodes());
                //_multipin_incidence_array[write_pos].id = previous;
                //_multipin_incidence_array[write_pos].multiplicity = count;
                _multipin_id[write_pos] = previous;
                _multipin_multiplicity[write_pos] = count;
                previous = incidence_array[i];
                count = 1;
                ++write_pos;
            } else {
                ++count;
            }
        }
        //_multipin_incidence_array[write_pos].id = previous;
        //_multipin_incidence_array[write_pos].multiplicity = count;
        _multipin_id[write_pos] = previous;
        _multipin_multiplicity[write_pos] = count;
        });
    _multipin_indexes[num_hyperedges] = multipin_mapping_prefix_sum[num_pins];
    chg._multipin_indexes = std::move(_multipin_indexes);
    //chg._multipin_incidence_array = std::move(_multipin_incidence_array);
    chg._multipin_id = std::move(_multipin_id);
    chg._multipin_multiplicity = std::move(_multipin_multiplicity);
    chg._tmp_community_hypergraph_buffer = _tmp_community_hypergraph_buffer;
    _tmp_community_hypergraph_buffer = nullptr;
    utils::Timer::instance().stop_timer("community_specific");
    utils::Timer::instance().stop_timer("community_hypergraph_contract");
    return chg;
}
}