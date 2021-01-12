#include "community_hypergraph.h"

#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/datastructures/concurrent_bucket_map.h"
#include "mt-kahypar/utils/timer.h"

#include <tbb/parallel_reduce.h>


namespace mt_kahypar::ds {

CommunityHypergraph CommunityHypergraph::contract(/*parallel::scalable_vector<HypernodeID>& communities*/) {
    Clustering communities = _hg->_community_ids;
    ASSERT(communities.size() == _hg->_num_hypernodes);

    if (!_tmp_contraction_buffer) {
        allocateTmpContractionBuffer();
    }

    // AUXILLIARY BUFFERS - Reused during multilevel hierarchy to prevent expensive allocations
    Array<size_t>& mapping = _tmp_contraction_buffer->mapping;
    Array<Hypernode>& tmp_hypernodes = _tmp_contraction_buffer->tmp_hypernodes;
    // IncidentNets& tmp_incident_nets = _tmp_contraction_buffer->tmp_incident_nets;
    Array<parallel::IntegralAtomicWrapper<size_t>>& tmp_num_incident_nets =
        _tmp_contraction_buffer->tmp_num_incident_nets;
    Array<parallel::IntegralAtomicWrapper<HypernodeWeight>>& hn_weights =
        _tmp_contraction_buffer->hn_weights;
    // Array<Hyperedge>& tmp_hyperedges = _tmp_contraction_buffer->tmp_hyperedges;
    // IncidenceArray& tmp_incidence_array = _tmp_contraction_buffer->tmp_incidence_array;
    // Array<size_t>& he_sizes = _tmp_contraction_buffer->he_sizes;
    // Array<size_t>& valid_hyperedges = _tmp_contraction_buffer->valid_hyperedges;

    ASSERT(static_cast<size_t>(_hg->_num_hypernodes) <= mapping.size());
    ASSERT(static_cast<size_t>(_hg->_num_hypernodes) <= tmp_hypernodes.size());
    //ASSERT(static_cast<size_t>(_hg->_total_degree) <= tmp_incident_nets.size());
    ASSERT(static_cast<size_t>(_hg->_num_hypernodes) <= tmp_num_incident_nets.size());
    ASSERT(static_cast<size_t>(_hg->_num_hypernodes) <= hn_weights.size());
    // ASSERT(static_cast<size_t>(_num_hyperedges) <= tmp_hyperedges.size());
    // ASSERT(static_cast<size_t>(_num_pins) <= tmp_incidence_array.size());
    // ASSERT(static_cast<size_t>(_num_hyperedges) <= he_sizes.size());
    // ASSERT(static_cast<size_t>(_num_hyperedges) <= valid_hyperedges.size());


    // #################### STAGE 1 ####################
    // Compute vertex ids of coarse hypergraph with a parallel prefix sum
    utils::Timer::instance().start_timer("preprocess_contractions", "Preprocess Contractions");
    mapping.assign(_hg->_num_hypernodes, 0);

    _hg->doParallelForAllNodes([&](const HypernodeID& hn) {
        ASSERT(static_cast<size_t>(communities[hn]) < mapping.size());
        mapping[communities[hn]] = 1UL;
        });

    // Prefix sum determines vertex ids in coarse hypergraph
    parallel::TBBPrefixSum<size_t, Array> mapping_prefix_sum(mapping);
    tbb::parallel_scan(tbb::blocked_range<size_t>(0UL, _hg->_num_hypernodes), mapping_prefix_sum);
    HypernodeID num_hypernodes = mapping_prefix_sum.total_sum();

    // Remap community ids
    tbb::parallel_for(ID(0), _hg->_num_hypernodes, [&](const HypernodeID& hn) {
        if (_hg->nodeIsEnabled(hn)) {
            communities[hn] = mapping_prefix_sum[communities[hn]];
        } else {
            communities[hn] = kInvalidHypernode;
        }

        // Reset tmp contraction buffer
        if ( hn < num_hypernodes ) {
          hn_weights[hn] = 0;
          tmp_hypernodes[hn] = Hypernode(true);
          tmp_num_incident_nets[hn] = 0;
        }
        });

    // Mapping from a vertex id of the current hypergraph to its
    // id in the coarse hypergraph
    auto map_to_coarse_hypergraph = [&](const HypernodeID hn) {
        ASSERT(hn < communities.size());
        return communities[hn];
    };


    _hg->doParallelForAllNodes([&](const HypernodeID& hn) {
        const HypernodeID coarse_hn = map_to_coarse_hypergraph(hn);
        ASSERT(coarse_hn < num_hypernodes, V(coarse_hn) << V(num_hypernodes));
        // Weight vector is atomic => thread-safe
        hn_weights[coarse_hn] += _hg->nodeWeight(hn);
        // In case community detection is enabled all vertices matched to one vertex
        // in the contracted hypergraph belong to same community. Otherwise, all communities
        // are default assigned to community 0
        // Aggregate upper bound for number of incident nets of the contracted vertex
        tmp_num_incident_nets[coarse_hn] += _hg->nodeDegree(hn);
        });
    utils::Timer::instance().stop_timer("preprocess_contractions");

    return std::move(*this);
}

}