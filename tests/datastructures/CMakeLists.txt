target_sources(mt_kahypar_fast_tests PRIVATE
        static_hypergraph_test.cc
        partitioned_hypergraph_test.cc
        delta_partitioned_hypergraph_test.cc
        partitioned_hypergraph_smoke_test.cc
        gain_updates_test.cc
        graph_test.cc
        connectivity_set_test.cc
        priority_queue_test.cc
        sparsifier_hypergraph_test.cc
        array_test.cc
        sparse_map_test.cc
        pin_count_in_part_test.cc
)

target_sources(mt_kahypar_strong_tests PRIVATE
        dynamic_hypergraph_test.cc
        partitioned_hypergraph_test.cc
        delta_partitioned_hypergraph_test.cc
        partitioned_hypergraph_smoke_test.cc
        contraction_tree_test.cc
        gain_updates_test.cc
        incident_net_array_test.cc
        nlevel_smoke_test.cc
)
