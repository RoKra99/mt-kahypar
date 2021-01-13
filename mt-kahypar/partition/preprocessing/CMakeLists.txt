set(PreprocessingSources
        community_detection/parallel_louvain.cpp
        community_detection/local_moving_modularity.cpp
        community_detection/hypergraph_local_moving_modularity.cpp
        community_detection/hypergraph_louvain.cpp)

foreach(modtarget IN LISTS TARGETS_WANTING_ALL_SOURCES)
    target_sources(${modtarget} PRIVATE ${PreprocessingSources})
endforeach()