#include <memory>
#include <vector>
#include <iostream>

#include <libkahypar.h>

// Install library interface via 'sudo make install.library' in build folder
// Compile with: g++ -std=c++14 -DNDEBUG -O3 partition_hypergraph_lib_2.cc -o example_2 -lkahypar
int main(int argc, char* argv[]) {

  // Initialize thread pool with 8 threads and NUMA allocation policy INTERLEAVED
  mt_kahypar_initialize_thread_pool(8, true /* activate interleaved NUMA allocation policy */ );

  // Load context from file
  mt_kahypar_context_t* context = mt_kahypar_context_new();
  mt_kahypar_configure_context_from_file(context, "../config/quality_preset.ini");

  // Setup Hypergraph
  mt_kahypar_hypernode_id_t num_vertices = 0;
  mt_kahypar_hyperedge_id_t num_hyperedges = 0;
  size_t* hyperedge_indices(nullptr);
  mt_kahypar_hyperedge_id_t* hyperedges(nullptr);
  mt_kahypar_hypernode_weight_t* hypernode_weights(nullptr);
  mt_kahypar_hyperedge_weight_t* hyperedge_weights(nullptr);
  mt_kahypar_read_hypergraph_from_file("contracted_ibm01.hgr", &num_vertices, &num_hyperedges,
    &hyperedge_indices, &hyperedges, &hyperedge_weights, &hypernode_weights);

  const double imbalance = 0.03;
  const mt_kahypar_partition_id_t k = 2;

  mt_kahypar_hyperedge_weight_t objective = 0;

  std::vector<mt_kahypar_partition_id_t> partition(num_vertices, -1);

  // Partition Hypergraph
  mt_kahypar_partition(num_vertices, num_hyperedges,
       	               imbalance, k, 0 /* seed */,
               	       nullptr /* unit vertex_weights */, hyperedge_weights,
               	       hyperedge_indices, hyperedges,
       	               &objective, context, partition.data(),
                       false /* verbose output */ );

  // Print objective
  std::cout << "Objective: " << objective << std::endl;

  delete(hyperedge_indices);
  delete(hyperedges);
  delete(hypernode_weights);
  delete(hyperedge_weights);
  mt_kahypar_context_free(context);
}