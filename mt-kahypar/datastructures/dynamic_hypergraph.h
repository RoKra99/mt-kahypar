/*******************************************************************************
 * This file is part of KaHyPar.
 *
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * KaHyPar is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * KaHyPar is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with KaHyPar.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#pragma once

#include <mutex>

#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_scan.h"

#include "kahypar/meta/mandatory.h"

#include "mt-kahypar/datastructures/community_support.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/incident_net_vector.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
namespace ds {

// Forward
class DynamicHypergraphFactory;

class DynamicHypergraph {

  static constexpr bool debug = false;
  static constexpr bool enable_heavy_assert = false;

  static_assert(std::is_unsigned<HypernodeID>::value, "Hypernode ID must be unsigned");
  static_assert(std::is_unsigned<HyperedgeID>::value, "Hyperedge ID must be unsigned");

  /**
   * Represents a hypernode of the hypergraph and contains all information
   * associated with a vertex.
   */
  class Hypernode {
   public:
    using IDType = HypernodeID;

    Hypernode() :
      _weight(1),
      _community_id(0),
      _valid(false) { }

    Hypernode(const bool valid) :
      _weight(1),
      _community_id(0),
      _valid(valid) { }

    bool isDisabled() const {
      return _valid == false;
    }

    void enable() {
      ASSERT(isDisabled());
      _valid = true;
    }

    void disable() {
      ASSERT(!isDisabled());
      _valid = false;
    }

    HyperedgeWeight weight() const {
      return _weight;
    }

    void setWeight(HyperedgeWeight weight) {
      ASSERT(!isDisabled());
      _weight = weight;
    }

    PartitionID communityID() const {
      ASSERT(!isDisabled());
      return _community_id;
    }

    void setCommunityID(const PartitionID community_id) {
      ASSERT(!isDisabled());
      _community_id = community_id;
    }

   private:
    // ! Hypernode weight
    HyperedgeWeight _weight;
    // ! Community id
    PartitionID _community_id;
    // ! Flag indicating whether or not the element is active.
    bool _valid;
  };

  /**
   * Represents a hyperedge of the hypergraph and contains all information
   * associated with a net (except connectivity information).
   */
  class Hyperedge {
   public:
    using IDType = HyperedgeID;

    Hyperedge() :
      _begin(0),
      _size(0),
      _weight(1),
      _hash(kEdgeHashSeed),
      _valid(false) { }

    // Sentinel Constructor
    Hyperedge(const size_t begin) :
      _begin(begin),
      _size(0),
      _weight(1),
      _hash(kEdgeHashSeed),
      _valid(false) { }

    // ! Disables the hypernode/hyperedge. Disable hypernodes/hyperedges will be skipped
    // ! when iterating over the set of all nodes/edges.
    void disable() {
      ASSERT(!isDisabled());
      _valid = false;
    }

    void enable() {
      ASSERT(isDisabled());
      _valid = true;
    }

    bool isDisabled() const {
      return _valid == false;
    }

    // ! Returns the index of the first element in _incidence_array
    size_t firstEntry() const {
      return _begin;
    }

    // ! Sets the index of the first element in _incidence_array to begin
    void setFirstEntry(size_t begin) {
      ASSERT(!isDisabled());
      _begin = begin;
    }

    // ! Returns the index of the first element in _incidence_array
    size_t firstInvalidEntry() const {
      return _begin + _size;
    }

    size_t size() const {
      ASSERT(!isDisabled());
      return _size;
    }

    void setSize(size_t size) {
      ASSERT(!isDisabled());
      _size = size;
    }

    void decrementSize() {
      ASSERT(!isDisabled());
      --_size;
    }

    HyperedgeWeight weight() const {
      ASSERT(!isDisabled());
      return _weight;
    }

    void setWeight(HyperedgeWeight weight) {
      ASSERT(!isDisabled());
      _weight = weight;
    }

    size_t& hash() {
      return _hash;
    }

    size_t hash() const {
      return _hash;
    }

    bool operator== (const Hyperedge& rhs) const {
      return _begin == rhs._begin && _size == rhs._size && _weight == rhs._weight;
    }

    bool operator!= (const Hyperedge& rhs) const {
      return _begin != rhs._begin || _size != rhs._size || _weight != rhs._weight;
    }

   private:
    // ! Index of the first element in _incidence_array
    size_t _begin;
    // ! Number of pins
    size_t _size;
    // ! hyperedge weight
    HyperedgeWeight _weight;
    // ! Hash of pins
    size_t _hash;
    // ! Flag indicating whether or not the element is active.
    bool _valid;
  };

  /*!
   * Iterator for HypergraphElements (Hypernodes/Hyperedges)
   *
   * The iterator is used in for-each loops over all hypernodes/hyperedges.
   * In order to support iteration over coarsened hypergraphs, this iterator
   * skips over HypergraphElements marked as invalid.
   * Iterating over the set of vertices \f$V\f$ therefore is linear in the
   * size \f$|V|\f$ of the original hypergraph - even if it has been coarsened
   * to much smaller size. The same also holds true for for-each loops over
   * the set of hyperedges.
   *
   * In order to be as generic as possible, the iterator does not expose the
   * internal Hypernode/Hyperedge representations. Instead only handles to
   * the respective elements are returned, i.e. the IDs of the corresponding
   * hypernodes/hyperedges.
   *
   */
  template <typename ElementType>
  class HypergraphElementIterator :
    public std::iterator<std::forward_iterator_tag,    // iterator_category
                         typename ElementType::IDType,   // value_type
                         std::ptrdiff_t,   // difference_type
                         const typename ElementType::IDType*,   // pointer
                         typename ElementType::IDType> {   // reference
   public:
    using IDType = typename ElementType::IDType;

    /*!
     * Construct a HypergraphElementIterator
     * See GenericHypergraph::nodes() or GenericHypergraph::edges() for usage.
     *
     * If start_element is invalid, the iterator advances to the first valid
     * element.
     *
     * \param start_element A pointer to the starting position
     * \param id The index of the element the pointer points to
     * \param max_id The maximum index allowed
     */
    HypergraphElementIterator(const ElementType* start_element, IDType id, IDType max_id) :
      _id(id),
      _max_id(max_id),
      _element(start_element) {
      if (_id != _max_id && _element->isDisabled()) {
        operator++ ();
      }
    }

    // ! Returns the id of the element the iterator currently points to.
    IDType operator* () const {
      return _id;
    }

    // ! Prefix increment. The iterator advances to the next valid element.
    HypergraphElementIterator & operator++ () {
      ASSERT(_id < _max_id);
      do {
        ++_id;
        ++_element;
      } while (_id < _max_id && _element->isDisabled());
      return *this;
    }

    // ! Postfix increment. The iterator advances to the next valid element.
    HypergraphElementIterator operator++ (int) {
      HypergraphElementIterator copy = *this;
      operator++ ();
      return copy;
    }

    bool operator!= (const HypergraphElementIterator& rhs) {
      return _id != rhs._id;
    }

    bool operator== (const HypergraphElementIterator& rhs) {
      return _id == rhs._id;
    }

   private:
    // Handle to the HypergraphElement the iterator currently points to
    IDType _id = 0;
    // Maximum allowed index
    IDType _max_id = 0;
    // HypergraphElement the iterator currently points to
    const ElementType* _element = nullptr;
  };

  static_assert(std::is_trivially_copyable<Hypernode>::value, "Hypernode is not trivially copyable");
  static_assert(std::is_trivially_copyable<Hyperedge>::value, "Hyperedge is not trivially copyable");

  enum class ContractionResult : uint8_t {
    CONTRACTED = 0,
    PENDING_CONTRACTIONS = 1,
    WEIGHT_LIMIT_REACHED = 2
  };

  using IncidenceArray = Array<HypernodeID>;
  using IncidentNets = parallel::scalable_vector<IncidentNetVector<HyperedgeID>>;
  using OwnershipVector = parallel::scalable_vector<parallel::IntegralAtomicWrapper<bool>>;
  using ReferenceCountVector = parallel::scalable_vector<HypernodeID>;
  using MementoVector = parallel::scalable_vector<Memento>;
  using ThreadLocalHyperedgeVector = tbb::enumerable_thread_specific<parallel::scalable_vector<HyperedgeID>>;

 public:
  static constexpr bool is_static_hypergraph = true;
  static constexpr bool is_partitioned = false;
  static constexpr size_t SIZE_OF_HYPERNODE = sizeof(Hypernode);
  static constexpr size_t SIZE_OF_HYPEREDGE = sizeof(Hyperedge);

  // ! Iterator to iterate over the hypernodes
  using HypernodeIterator = HypergraphElementIterator<const Hypernode>;
  // ! Iterator to iterate over the hyperedges
  using HyperedgeIterator = HypergraphElementIterator<const Hyperedge>;
  // ! Iterator to iterate over the pins of a hyperedge
  using IncidenceIterator = typename IncidenceArray::const_iterator;
  // ! Iterator to iterate over the incident nets of a hypernode
  using IncidentNetsIterator = typename IncidentNetVector<HyperedgeID>::const_iterator;
  // ! Iterator to iterate over the set of communities contained in a hyperedge
  using CommunityIterator = typename CommunitySupport<StaticHypergraph>::CommunityIterator;


  explicit DynamicHypergraph() :
    _num_hypernodes(0),
    _num_removed_hypernodes(0),
    _num_hyperedges(0),
    _num_removed_hyperedges(0),
    _max_edge_size(0),
    _num_pins(0),
    _total_degree(0),
    _total_weight(0),
    _hypernodes(),
    _incident_nets(),
    _acquired_hns(),
    _hn_ref_count(),
    _contraction_tree(),
    _hyperedges(),
    _incidence_array(),
    _acquired_hes(),
    _tmp_incident_nets(),
    _failed_hyperedge_contractions(),
    _community_support() { }

  DynamicHypergraph(const DynamicHypergraph&) = delete;
  DynamicHypergraph & operator= (const DynamicHypergraph &) = delete;

  DynamicHypergraph(DynamicHypergraph&& other) :
    _num_hypernodes(other._num_hypernodes),
    _num_removed_hypernodes(other._num_removed_hypernodes),
    _num_hyperedges(other._num_hyperedges),
    _num_removed_hyperedges(other._num_removed_hyperedges),
    _max_edge_size(other._max_edge_size),
    _num_pins(other._num_pins),
    _total_degree(other._total_degree),
    _total_weight(other._total_weight),
    _hypernodes(std::move(other._hypernodes)),
    _incident_nets(std::move(other._incident_nets)),
    _acquired_hns(std::move(other._acquired_hns)),
    _hn_ref_count(std::move(other._hn_ref_count)),
    _contraction_tree(std::move(other._contraction_tree)),
    _hyperedges(std::move(other._hyperedges)),
    _incidence_array(std::move(other._incidence_array)),
    _acquired_hes(std::move(other._acquired_hes)),
    _tmp_incident_nets(std::move(other._tmp_incident_nets)),
    _failed_hyperedge_contractions(std::move(other._failed_hyperedge_contractions)),
    _community_support(std::move(other._community_support)) { }

  DynamicHypergraph & operator= (DynamicHypergraph&& other) {
    _num_hypernodes = other._num_hypernodes;
    _num_removed_hypernodes = other._num_removed_hypernodes;
    _num_hyperedges = other._num_hyperedges;
    _num_removed_hyperedges = other._num_removed_hyperedges;
    _max_edge_size = other._max_edge_size;
    _num_pins = other._num_pins;
    _total_degree = other._total_degree;
    _total_weight = other._total_weight;
    _hypernodes = std::move(other._hypernodes);
    _incident_nets = std::move(other._incident_nets);
    _acquired_hns = std::move(other._acquired_hns);
    _hn_ref_count = std::move(other._hn_ref_count);
    _contraction_tree = std::move(other._contraction_tree);
    _hyperedges = std::move(other._hyperedges);
    _incidence_array = std::move(other._incidence_array);
    _acquired_hes = std::move(other._acquired_hes);
    _tmp_incident_nets = std::move(other._tmp_incident_nets);
    _failed_hyperedge_contractions = std::move(other._failed_hyperedge_contractions);
    _community_support = std::move(other._community_support);
    return *this;
  }

  ~DynamicHypergraph() {
    freeInternalData();
  }

  // ####################### General Hypergraph Stats #######################

  // ! Initial number of hypernodes
  HypernodeID initialNumNodes() const {
    return _num_hypernodes;
  }

  // ! Number of removed hypernodes
  HypernodeID numRemovedHypernodes() const {
    return _num_removed_hypernodes;
  }

  // ! Initial number of hyperedges
  HyperedgeID initialNumEdges() const {
    return _num_hyperedges;
  }

  HyperedgeID numGraphEdges() const {
    ERROR("numGraphEdges() is not supported in dynamic hypergraph");
    return 0;
  }

  HyperedgeID numNonGraphEdges() const {
    ERROR("numNonGraphEdges() is not supported in dynamic hypergraph");
    return initialNumEdges();
  }

  // ! Number of removed hyperedges
  HyperedgeID numRemovedHyperedges() const {
    return _num_removed_hyperedges;
  }

  // ! Set the number of removed hyperedges
  void setNumRemovedHyperedges(const HyperedgeID num_removed_hyperedges) {
    _num_removed_hyperedges = num_removed_hyperedges;
  }

  // ! Initial number of pins
  HypernodeID initialNumPins() const {
    return _num_pins;
  }

  // ! Initial sum of the degree of all vertices
  HypernodeID initialTotalVertexDegree() const {
    return _total_degree;
  }

  // ! Total weight of hypergraph
  HypernodeWeight totalWeight() const {
    return _total_weight;
  }

  // ! Recomputes the total weight of the hypergraph (parallel)
  void updateTotalWeight(const TaskGroupID) {
    _total_weight = tbb::parallel_reduce(tbb::blocked_range<HypernodeID>(ID(0), _num_hypernodes), 0,
      [this](const tbb::blocked_range<HypernodeID>& range, HypernodeWeight init) {
        HypernodeWeight weight = init;
        for (HypernodeID hn = range.begin(); hn < range.end(); ++hn) {
          weight += this->_hypernodes[hn].weight();
        }
        return weight;
      }, std::plus<HypernodeWeight>());
  }

  // ! Recomputes the total weight of the hypergraph (sequential)
  void updateTotalWeight() {
    _total_weight = 0;
    for ( const HypernodeID& hn : nodes() ) {
      _total_weight += nodeWeight(hn);
    }
  }

  // ####################### Iterators #######################

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const F& f) {
    static_cast<const DynamicHypergraph&>(*this).doParallelForAllNodes(f);
  }

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const F& f) const {
    tbb::parallel_for(ID(0), _num_hypernodes, [&](const HypernodeID& hn) {
      if ( nodeIsEnabled(hn) ) {
        f(hn);
      }
    });
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const F& f) {
    static_cast<const DynamicHypergraph&>(*this).doParallelForAllEdges(f);
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const F& f) const {
    tbb::parallel_for(ID(0), _num_hyperedges, [&](const HyperedgeID& he) {
      if ( edgeIsEnabled(he) ) {
        f(he);
      }
    });
  }

  // ! Returns a range of the active nodes of the hypergraph
  IteratorRange<HypernodeIterator> nodes() const {
    return IteratorRange<HypernodeIterator>(
      HypernodeIterator(_hypernodes.data(), ID(0), _num_hypernodes),
      HypernodeIterator(_hypernodes.data() + _num_hypernodes, _num_hypernodes, _num_hypernodes));
  }

  // ! Returns a range of the active edges of the hypergraph
  IteratorRange<HyperedgeIterator> edges() const {
    return IteratorRange<HyperedgeIterator>(
      HyperedgeIterator(_hyperedges.data(), ID(0), _num_hyperedges),
      HyperedgeIterator(_hyperedges.data() + _num_hyperedges, _num_hyperedges, _num_hyperedges));
  }

  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetsIterator> incidentEdges(const HypernodeID u, const bool is_disabled = false) const {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    ASSERT(is_disabled || !hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return IteratorRange<IncidentNetsIterator>(
      _incident_nets[u].cbegin(), _incident_nets[u].cend());
  }

  // ! Returns a range to loop over the pins of hyperedge e.
  IteratorRange<IncidenceIterator> pins(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    const Hyperedge& he = hyperedge(e);
    return IteratorRange<IncidenceIterator>(
      _incidence_array.cbegin() + he.firstEntry(),
      _incidence_array.cbegin() + he.firstInvalidEntry());
  }

  // ! Returns a range to loop over the pins of hyperedge e that belong to a certain community.
  // ! Note, this function fails if community hyperedges are not initialized.
  IteratorRange<IncidenceIterator> pins(const HyperedgeID e, const PartitionID community_id) const {
    return _community_support.pins(*this, e, community_id);
  }

  // ! Returns a range to loop over the set of communities contained in hyperedge e.
  IteratorRange<CommunityIterator> communities(const HyperedgeID e) const {
    return _community_support.communities(e);
  }

  // ####################### Hypernode Information #######################

  // ! Weight of a vertex
  HypernodeWeight nodeWeight(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).weight();
  }

  // ! Sets the weight of a vertex
  void setNodeWeight(const HypernodeID u, const HypernodeWeight weight) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).setWeight(weight);
  }

  // ! Degree of a hypernode
  HyperedgeID nodeDegree(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return _incident_nets[u].size();
  }

  // ! Returns, whether a hypernode is enabled or not
  bool nodeIsEnabled(const HypernodeID u) const {
    return !hypernode(u).isDisabled();
  }

  // ! Enables a hypernode (must be disabled before)
  void enableHypernode(const HypernodeID u) {
    hypernode(u).enable();
  }

  // ! Disables a hypernode (must be enabled before)
  void disableHypernode(const HypernodeID u) {
    hypernode(u).disable();
  }

  // ! Removes a hypernode (must be enabled before)
  void removeHypernode(const HypernodeID u) {
    hypernode(u).disable();
    ++_num_removed_hypernodes;
  }

  // ####################### Hyperedge Information #######################

  // ! Weight of a hyperedge
  HypernodeWeight edgeWeight(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    return hyperedge(e).weight();
  }

  // ! Sets the weight of a hyperedge
  void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    return hyperedge(e).setWeight(weight);
  }

  // ! Number of pins of a hyperedge
  HypernodeID edgeSize(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    return hyperedge(e).size();
  }

  // ! Maximum size of a hyperedge
  HypernodeID maxEdgeSize() const {
    return _max_edge_size;
  }

  // ! Hash value defined over the pins of a hyperedge
  size_t edgeHash(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    return hyperedge(e).hash();
  }

  // ! Returns, whether a hyperedge is enabled or not
  bool edgeIsEnabled(const HyperedgeID e) const {
    return !hyperedge(e).isDisabled();
  }

  // ! Enables a hyperedge (must be disabled before)
  void enableHyperedge(const HyperedgeID e) {
    hyperedge(e).enable();
  }

  // ! Disabled a hyperedge (must be enabled before)
  void disableHyperedge(const HyperedgeID e) {
    hyperedge(e).disable();
  }

  HyperedgeID graphEdgeID(const HyperedgeID) const {
    ERROR("graphEdgeID(e) is not supported in dynamic hypergraph");
    return kInvalidHyperedge;
  }

  HyperedgeID nonGraphEdgeID(const HyperedgeID) const {
    ERROR("nonGraphEdgeID(e) is not supported in dynamic hypergraph");
    return kInvalidHyperedge;
  }

  HypernodeID graphEdgeHead(const HyperedgeID, const HypernodeID) const {
    ERROR("nonGraphEdgeID(e) is not supported in dynamic hypergraph");
    return kInvalidHyperedge;
  }

  // ####################### Community Hyperedge Information #######################

  // ! Weight of a community hyperedge
  HypernodeWeight edgeWeight(const HyperedgeID e, const PartitionID community_id) const {
    return _community_support.edgeWeight(e, community_id);
  }

  // ! Sets the weight of a community hyperedge
  void setEdgeWeight(const HyperedgeID e, const PartitionID community_id, const HyperedgeWeight weight) {
    _community_support.setEdgeWeight(e, community_id, weight);
  }

  // ! Number of pins of a hyperedge that are assigned to a community
  HypernodeID edgeSize(const HyperedgeID e, const PartitionID community_id) const {
    return _community_support.edgeSize(e, community_id);
  }

  // ! Hash value defined over the pins of a hyperedge that belongs to a community
  size_t edgeHash(const HyperedgeID e, const PartitionID community_id) const {
    return _community_support.edgeHash(e, community_id);
  }

  // ####################### Community Information #######################

  // ! Number of communities
  PartitionID numCommunities() const {
    return _community_support.numCommunities();
  }

  // ! Community id which hypernode u is assigned to
  PartitionID communityID(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).communityID();
  }

  // ! Assign a community to a hypernode
  // ! Note, in order to use all community-related functions, initializeCommunities()
  // ! have to be called after assigning to each vertex a community id
  void setCommunityID(const HypernodeID u, const PartitionID community_id) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).setCommunityID(community_id);
  }

  // ! Consider hypernode u is part of community C = {v_1, ..., v_n},
  // ! than this function returns a unique id for hypernode u in the
  // ! range [0,n).
  HypernodeID communityNodeId(const HypernodeID u) const {
    return _community_support.communityNodeId(u);
  }

  // ! Number of hypernodes in community
  HypernodeID numCommunityHypernodes(const PartitionID community) const {
    return _community_support.numCommunityHypernodes(community);
  }

  // ! Number of pins in community
  HypernodeID numCommunityPins(const PartitionID community) const {
    return _community_support.numCommunityPins(community);
  }

  // ! Total degree of community
  HyperedgeID communityDegree(const PartitionID community) const {
    return _community_support.communityDegree(community);
  }

  // ! Number of communities which pins of hyperedge belongs to
  size_t numCommunitiesInHyperedge(const HyperedgeID e) const {
    return _community_support.numCommunitiesInHyperedge(e);
  }

  // ####################### Contract / Uncontract #######################

  DynamicHypergraph contract(parallel::scalable_vector<HypernodeID>&,
                             const TaskGroupID) {
    ERROR("contract(c, id) is not supported in dynamic hypergraph");
    return DynamicHypergraph();
  }

  /**!
   * Registers a contraction in the hypergraph whereas vertex u is the representative
   * of the contraction and v its contraction partner. Several threads can call this function
   * in parallel. The function adds the contraction of u and v to a contraction tree that determines
   * a parallel execution order and synchronization points for all running contractions.
   * The contraction can be executed by calling function contract(v, max_node_weight).
   */
  bool registerContraction(const HypernodeID u, const HypernodeID v) {
    // Acquires ownership of vertex v that gives the calling thread exclusive rights
    // to modify the contraction tree entry of v
    acquireHypernode(v);

    // If there is no other contraction registered for vertex v
    // we try to determine its representative in the contraction tree
    if ( _contraction_tree[v] == v ) {

      HypernodeID w = u;
      bool cycle_detected = false;
      while ( true ) {
        // Search for representative of u in the contraction tree.
        // It is either a root of the contraction tree or a vertex
        // with a reference count greater than zero, which indicates
        // that there are still ongoing contractions on this node that
        // have to be processed.
        while ( _contraction_tree[w] != w && _hn_ref_count[w] == 0 ) {
          w = _contraction_tree[w];
          if ( w == v ) {
            cycle_detected = true;
            break;
          }
        }

        if ( !cycle_detected ) {
          // In case contraction of u and v does not induce any
          // cycle in the contraction tree we try to acquire vertex w
          if ( w < v ) {
            // Acquire ownership in correct order to prevent deadlocks
            releaseHypernode(v);
            acquireHypernode(w);
            acquireHypernode(v);
            if ( _contraction_tree[v] != v ) {
              releaseHypernode(v);
              return false;
            }
          } else {
            acquireHypernode(w);
          }

          // Double-check condition of while loop above after acquiring
          // ownership of w
          if ( _contraction_tree[w] != w && _hn_ref_count[w] == 0 ) {
            // In case something changed, we release ownership of w and
            // search again for the representative of u.
            releaseHypernode(w);
          } else {
            // Otherwise we perform final cycle check to verify that
            // contraction of u and v will not introduce any new cycle.
            HypernodeID x = w;
            do {
              x = _contraction_tree[x];
              if ( x == v ) {
                cycle_detected = true;
                break;
              }
            } while ( _contraction_tree[x] != x );

            if ( cycle_detected ) {
              releaseHypernode(w);
              releaseHypernode(v);
              return false;
            }

            // All checks succeded, we can safely increment the
            // reference count of w and update the contraction tree
            break;
          }
        } else {
          releaseHypernode(v);
          return false;
        }
      }

      // Increment reference count of w indicating that there pending
      // contraction at vertex w and update contraction tree.
      ++_hn_ref_count[w];
      _contraction_tree[v] = w;

      releaseHypernode(w);
      releaseHypernode(v);
      return true;
    } else {
      releaseHypernode(v);
      return false;
    }
  }

  /**!
   * Contracts a previously registered contraction. Representative u of vertex v is looked up
   * in the contraction tree and performed if there are no pending contractions in the subtree
   * of v and the contractions respects the maximum allowed node weight. If (u,v) is the last
   * pending contraction in the subtree of u then the function recursively contracts also
   * u (if any contraction is registered). Therefore, function can return several contractions
   * or also return an empty contraction vector.
   */
  MementoVector contract(const HypernodeID v,
                         const HypernodeWeight max_node_weight = std::numeric_limits<HypernodeWeight>::max()) {
    ASSERT(_contraction_tree[v] != v, "No contraction registered for hypernode" << v);

    MementoVector mementos;
    HypernodeID x = _contraction_tree[v];
    HypernodeID y = v;
    ContractionResult res = ContractionResult::CONTRACTED;
    // We perform all contractions registered in the contraction tree
    // as long as there are no pending contractions (_hn_ref_count[y] == 0
    // is equivalent with no pending contractions)
    while ( x != y && res != ContractionResult::PENDING_CONTRACTIONS) {
      // Perform Contraction
      res = contract(x, y, max_node_weight);
      if ( res == ContractionResult::CONTRACTED ) {
        mementos.emplace_back(Memento { x, y });
      }
      y = x;
      x = _contraction_tree[y];
    }
    return mementos;
  }

  // ! Only for testing
  HypernodeID contractionTree(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return _contraction_tree[u];
  }

  // ! Only for testing
  HypernodeID referenceCount(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return _hn_ref_count[u];
  }

  // ! Only for testing
  void decrementReferenceCount(const HypernodeID u) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    --_hn_ref_count[u];
  }

  // ####################### Remove / Restore Hyperedges #######################

  /*!
  * Removes a hyperedge from the hypergraph. This includes the removal of he from all
  * of its pins and to disable the hyperedge.
  *
  * NOTE, this function is not thread-safe and should only be called in a single-threaded
  * setting.
  */
  void removeEdge(const HyperedgeID he) {
    ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");
    for ( const HypernodeID& pin : pins(he) ) {
      removeIncidentEdgeFromHypernode(he, pin);
    }
    ++_num_removed_hyperedges;
    disableHyperedge(he);
  }

  /*!
  * Removes a hyperedge from the hypergraph. This includes the removal of he from all
  * of its pins and to disable the hyperedge. Note, in contrast to removeEdge, this function
  * removes hyperedge from all its pins in parallel.
  *
  * NOTE, this function is not thread-safe and should only be called in a single-threaded
  * setting.
  */
  void removeLargeEdge(const HyperedgeID he) {
    ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");
    const size_t incidence_array_start = hyperedge(he).firstEntry();
    const size_t incidence_array_end = hyperedge(he).firstInvalidEntry();
    tbb::parallel_for(incidence_array_start, incidence_array_end, [&](const size_t pos) {
      const HypernodeID pin = _incidence_array[pos];
      removeIncidentEdgeFromHypernode(he, pin);
    });
    disableHyperedge(he);
  }

  /*!
   * Restores a large hyperedge previously removed from the hypergraph.
   */
  void restoreLargeEdge(const HyperedgeID& he) {
    ASSERT(!edgeIsEnabled(he), "Hyperedge" << he << "is enabled");
    enableHyperedge(he);
    const size_t incidence_array_start = hyperedge(he).firstEntry();
    const size_t incidence_array_end = hyperedge(he).firstInvalidEntry();
    tbb::parallel_for(incidence_array_start, incidence_array_end, [&](const size_t pos) {
      const HypernodeID pin = _incidence_array[pos];
      connectHypernodeWithIncidentEdge(he, pin);
    });
  }

  // ####################### Initialization / Reset Functions #######################

  /*!
   * Initializes community-related information after all vertices are assigned to a community.
   * This includes:
   *  1.) Number of Communities
   *  2.) Number of Vertices per Community
   *  3.) Number of Pins per Community
   *  4.) For each hypernode v of community C, we compute a unique id within
   *      that community in the range [0, |C|)
   */
  void initializeCommunities() {
    _community_support.initialize(*this);
  }

  /*!
  * Initializes community hyperedges.
  * This includes:
  *   1.) Sort the pins of each hyperedge in increasing order of their community id
  *   2.) Introduce for each community id contained in a hyperedge a seperate
  *       community hyperedge pointing to a range of consecutive pins with
  *       same community in that hyperedge
  */
  void initializeCommunityHyperedges(const TaskGroupID) {
    _community_support.initializeCommunityHyperedges(*this);
  }

  /*!
   * Removes all community hyperedges from the hypergraph after parallel community
   * coarsening terminates.
   */
  void removeCommunityHyperedges(const TaskGroupID,
                                 const parallel::scalable_vector<HypernodeID>& contraction_index = {}) {
    _community_support.removeCommunityHyperedges(contraction_index);
  }

  // ! Reset internal community information
  void setCommunityIDs(const parallel::scalable_vector<PartitionID>& community_ids) {
    if ( _community_support.isInitialized() ) {
      _community_support.freeInternalData();
    }

    ASSERT(community_ids.size() == UI64(_num_hypernodes));
    doParallelForAllNodes([&](const HypernodeID& hn) {
      hypernode(hn).setCommunityID(community_ids[hn]);
    });

    initializeCommunities();
  }

  // ####################### Copy #######################

  // ! Copy dynamic hypergraph in parallel
  DynamicHypergraph copy(const TaskGroupID task_group_id) {
    DynamicHypergraph hypergraph;

    hypergraph._num_hypernodes = _num_hypernodes;
    hypergraph._num_removed_hypernodes = _num_removed_hypernodes;
    hypergraph._num_hyperedges = _num_hyperedges;
    hypergraph._num_removed_hyperedges = _num_removed_hyperedges;
    hypergraph._max_edge_size = _max_edge_size;
    hypergraph._num_pins = _num_pins;
    hypergraph._total_degree = _total_degree;
    hypergraph._total_weight = _total_weight;

    tbb::parallel_invoke([&] {
      hypergraph._hypernodes.resize(_hypernodes.size());
      memcpy(hypergraph._hypernodes.data(), _hypernodes.data(),
        sizeof(Hypernode) * _hypernodes.size());
    }, [&] {
      tbb::parallel_invoke([&] {
        hypergraph._incident_nets.resize(_incident_nets.size());
      }, [&] {
        hypergraph._acquired_hns.resize(_acquired_hns.size());
      });
      tbb::parallel_for(ID(0), _num_hypernodes, [&](const HypernodeID& hn) {
        hypergraph._incident_nets[hn].resize(_incident_nets[hn].size());
        hypergraph._acquired_hns[hn] = _acquired_hns[hn];
        memcpy(hypergraph._incident_nets[hn].data(), _incident_nets[hn].data(),
          sizeof(HyperedgeID) * _incident_nets[hn].size());
      });
    }, [&] {
      hypergraph._hn_ref_count.resize(_hn_ref_count.size());
      memcpy(hypergraph._hn_ref_count.data(), _hn_ref_count.data(),
        sizeof(HypernodeID) * _hn_ref_count.size());
    }, [&] {
      hypergraph._contraction_tree.resize(_contraction_tree.size());
      memcpy(hypergraph._contraction_tree.data(), _contraction_tree.data(),
        sizeof(HypernodeID) * _contraction_tree.size());
    }, [&] {
      hypergraph._hyperedges.resize(_hyperedges.size());
      memcpy(hypergraph._hyperedges.data(), _hyperedges.data(),
        sizeof(Hyperedge) * _hyperedges.size());
    }, [&] {
      hypergraph._incidence_array.resize(_incidence_array.size());
      memcpy(hypergraph._incidence_array.data(), _incidence_array.data(),
        sizeof(HypernodeID) * _incidence_array.size());
    }, [&] {
      hypergraph._acquired_hes.resize(_num_hyperedges);
      tbb::parallel_for(ID(0), _num_hyperedges, [&](const HyperedgeID& he) {
        hypergraph._acquired_hes[he] = _acquired_hes[he];
      });
    }, [&] {
      hypergraph._community_support = _community_support.copy(task_group_id);
    });
    return hypergraph;
  }

  // ! Copy dynamic hypergraph sequential
  DynamicHypergraph copy() {
    DynamicHypergraph hypergraph;

    hypergraph._num_hypernodes = _num_hypernodes;
    hypergraph._num_removed_hypernodes = _num_removed_hypernodes;
    hypergraph._num_hyperedges = _num_hyperedges;
    hypergraph._num_removed_hyperedges = _num_removed_hyperedges;
    hypergraph._max_edge_size = _max_edge_size;
    hypergraph._num_pins = _num_pins;
    hypergraph._total_degree = _total_degree;
    hypergraph._total_weight = _total_weight;

    hypergraph._hypernodes.resize(_hypernodes.size());
    memcpy(hypergraph._hypernodes.data(), _hypernodes.data(),
      sizeof(Hypernode) * _hypernodes.size());
    hypergraph._incident_nets.resize(_incident_nets.size());
    hypergraph._acquired_hns.resize(_num_hypernodes);
    for ( HypernodeID hn = 0; hn < _num_hypernodes; ++hn ) {
      hypergraph._incident_nets[hn].resize(_incident_nets[hn].size());
      hypergraph._acquired_hns[hn] = _acquired_hns[hn];
      memcpy(hypergraph._incident_nets[hn].data(), _incident_nets[hn].data(),
        sizeof(HyperedgeID) * _incident_nets[hn].size());
    }
    hypergraph._hn_ref_count.resize(_hn_ref_count.size());
    memcpy(hypergraph._hn_ref_count.data(), _hn_ref_count.data(),
      sizeof(HypernodeID) * _hn_ref_count.size());
    hypergraph._contraction_tree.resize(_contraction_tree.size());
    memcpy(hypergraph._contraction_tree.data(), _contraction_tree.data(),
      sizeof(HypernodeID) * _contraction_tree.size());
    hypergraph._hyperedges.resize(_hyperedges.size());
    memcpy(hypergraph._hyperedges.data(), _hyperedges.data(),
      sizeof(Hyperedge) * _hyperedges.size());
    hypergraph._incidence_array.resize(_incidence_array.size());
    memcpy(hypergraph._incidence_array.data(), _incidence_array.data(),
      sizeof(HypernodeID) * _incidence_array.size());
    hypergraph._acquired_hes.resize(_num_hyperedges);
    for ( HyperedgeID he = 0; he < _num_hyperedges; ++he ) {
      hypergraph._acquired_hes[he] = _acquired_hes[he];
    }

    hypergraph._community_support = _community_support.copy();

    return hypergraph;
  }

  // ! Free internal data in parallel
  void freeInternalData() {
    if ( _num_hypernodes > 0 || _num_hyperedges > 0 ) {
      _community_support.freeInternalData();
    }
    _num_hypernodes = 0;
    _num_hyperedges = 0;
  }

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    parent->addChild("Hypernodes", sizeof(Hypernode) * _hypernodes.size());
    parent->addChild("Incident Nets", sizeof(HyperedgeID) * _incidence_array.size());
    parent->addChild("Hypernode Ownership Vector", sizeof(bool) * _acquired_hns.size());
    parent->addChild("Hypernode Reference Counts", sizeof(HypernodeID) * _hn_ref_count.size());
    parent->addChild("Contraction Tree", sizeof(HypernodeID) * _contraction_tree.size());
    parent->addChild("Hyperedges", sizeof(Hyperedge) * _hyperedges.size());
    parent->addChild("Incidence Array", sizeof(HypernodeID) * _incidence_array.size());
    parent->addChild("Hyperedge Ownership Vector", sizeof(bool) * _acquired_hes.size());

    utils::MemoryTreeNode* community_support_node = parent->addChild("Community Support");
    _community_support.memoryConsumption(community_support_node);
  }

 private:
  friend class DynamicHypergraphFactory;
  template<typename Hypergraph>
  friend class CommunitySupport;

  // ####################### Acquiring / Releasing Ownership #######################

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void acquireHypernode(const HypernodeID u) {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    bool expected = false;
    bool desired = true;
    while ( !_acquired_hns[u].compare_exchange_strong(expected, desired) ) {
      expected = false;
    }
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool tryAcquireHypernode(const HypernodeID u) {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    bool expected = false;
    bool desired = true;
    return _acquired_hns[u].compare_exchange_strong(expected, desired);
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void releaseHypernode(const HypernodeID u) {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    ASSERT(_acquired_hns[u], "Hypernode" << u << "is not acquired!");
    _acquired_hns[u] = false;
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void acquireHyperedge(const HyperedgeID e) {
    ASSERT(e < _num_hyperedges, "Hyperedge" << e << "does not exist");
    bool expected = false;
    bool desired = true;
    while ( !_acquired_hes[e].compare_exchange_strong(expected, desired) ) {
      expected = false;
    }
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE bool tryAcquireHyperedge(const HyperedgeID e) {
    ASSERT(e < _num_hyperedges, "Hyperedge" << e << "does not exist");
    bool expected = false;
    bool desired = true;
    return _acquired_hes[e].compare_exchange_strong(expected, desired);
  }

  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void releaseHyperedge(const HyperedgeID e) {
    ASSERT(e < _num_hyperedges, "Hyperedge" << e << "does not exist");
    ASSERT(_acquired_hes[e], "Hyperedge" << e << "is not acquired!");
    _acquired_hes[e] = false;
  }

  // ####################### Hypernode Information #######################

  // ! Accessor for hypernode-related information
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Hypernode& hypernode(const HypernodeID u) const {
    ASSERT(u <= _num_hypernodes, "Hypernode" << u << "does not exist");
    return _hypernodes[u];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hypernode& hypernode(const HypernodeID u) {
    return const_cast<Hypernode&>(static_cast<const DynamicHypergraph&>(*this).hypernode(u));
  }

  // ####################### Hyperedge Information #######################

  // ! Accessor for hyperedge-related information
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Hyperedge& hyperedge(const HyperedgeID e) const {
    ASSERT(e <= _num_hyperedges, "Hyperedge" << e << "does not exist");
    return _hyperedges[e];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hyperedge& hyperedge(const HyperedgeID e) {
    return const_cast<Hyperedge&>(static_cast<const DynamicHypergraph&>(*this).hyperedge(e));
  }

  // ####################### Contract / Uncontract #######################

  /**!
   * Contracts a previously registered contraction. The contraction of u and v is
   * performed if there are no pending contractions in the subtree of v and the
   * contractions respects the maximum allowed node weight. In case the contraction
   * was performed successfully, enum type CONTRACTED is returned. If contraction
   * was not performed either WEIGHT_LIMIT_REACHED (in case sum of both vertices is
   * greater than the maximum allowed node weight) or PENDING_CONTRACTIONS (in case
   * there are some unfinished contractions in the subtree of v) is returned.
   */
  ContractionResult contract(const HypernodeID u,
                             const HypernodeID v,
                             const HypernodeWeight max_node_weight) {

    // Acquire ownership in correct order to prevent deadlocks
    if ( u < v ) {
      acquireHypernode(u);
      acquireHypernode(v);
    } else {
      acquireHypernode(v);
      acquireHypernode(u);
    }

    // Contraction is valid if
    //  1.) Contraction partner v is enabled
    //  2.) There are no pending contractions on v
    //  3.) Resulting node weight is less or equal than a predefined upper bound
    const bool contraction_partner_valid = nodeIsEnabled(v) && _hn_ref_count[v] == 0;
    const bool less_or_equal_than_max_node_weight =
      hypernode(u).weight() + hypernode(v).weight() <= max_node_weight;
    if ( contraction_partner_valid && less_or_equal_than_max_node_weight ) {
      ASSERT(nodeIsEnabled(u), "Hypernode" << u << "is disabled!");
      hypernode(u).setWeight(nodeWeight(u) + nodeWeight(v));
      hypernode(v).disable();
      releaseHypernode(u);
      releaseHypernode(v);

      parallel::scalable_vector<HyperedgeID>& tmp_incident_nets = _tmp_incident_nets.local();
      parallel::scalable_vector<HyperedgeID>& failed_hyperedge_contractions = _failed_hyperedge_contractions.local();
      for ( const HyperedgeID& he : incidentEdges(v, true) ) {
        // Try to acquire ownership of hyperedge. In case of success, we perform the
        // contraction and otherwise, we remember the hyperedge and try later again.
        if ( tryAcquireHyperedge(he) ) {
          contractHyperedge(u, v, he, tmp_incident_nets);
          releaseHyperedge(he);
        } else {
          failed_hyperedge_contractions.push_back(he);
        }
      }

      // Perform contraction on which we failed to acquire ownership on the first try
      for ( const HyperedgeID& he : failed_hyperedge_contractions ) {
        acquireHyperedge(he);
        contractHyperedge(u, v, he, tmp_incident_nets);
        releaseHyperedge(he);
      }

      // tmp_incident_nets contains all hyperedges to which vertex u is
      // adjacent after the contraction, we use a special insert function to make
      // sure that iterators are not invalidated while inserting into the vector.
      if ( tmp_incident_nets.size() > 0 ) {
        _incident_nets[u].bulk_insert(tmp_incident_nets);
      }
      tmp_incident_nets.clear();
      failed_hyperedge_contractions.clear();

      acquireHypernode(u);
      --_hn_ref_count[u];
      releaseHypernode(u);
      return ContractionResult::CONTRACTED;
    } else {
      ContractionResult res = ContractionResult::PENDING_CONTRACTIONS;
      if ( !less_or_equal_than_max_node_weight ) {
        --_hn_ref_count[u];
        _contraction_tree[v] = v;
        res = ContractionResult::WEIGHT_LIMIT_REACHED;
      }
      releaseHypernode(u);
      releaseHypernode(v);
      return res;
    }
  }

  // ! Performs the contraction of (u,v) inside hyperedge he
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void contractHyperedge(const HypernodeID u, const HypernodeID v, const HyperedgeID he,
                                                         parallel::scalable_vector<HyperedgeID>& tmp_incident_nets) {
    Hyperedge& e = hyperedge(he);
    const HypernodeID pins_begin = e.firstEntry();
    const HypernodeID pins_end = e.firstInvalidEntry();
    HypernodeID slot_of_u = pins_end - 1;
    HypernodeID last_pin_slot = pins_end - 1;

    for (HypernodeID idx = pins_begin; idx != last_pin_slot; ++idx) {
      const HypernodeID pin = _incidence_array[idx];
      if (pin == v) {
        std::swap(_incidence_array[idx], _incidence_array[last_pin_slot]);
        --idx;
      } else if (pin == u) {
        slot_of_u = idx;
      }
    }

    ASSERT(_incidence_array[last_pin_slot] == v, "v is not last entry in adjacency array!");

    if (slot_of_u != last_pin_slot) {
      // Case 1:
      // Hyperedge e contains both u and v. Thus we don't need to connect u to e and
      // can just cut off the last entry in the edge array of e that now contains v.
      DBG << V(he) << ": Case 1";
      e.hash() -= kahypar::math::hash(v);
      e.decrementSize();
    } else {
      DBG << V(he) << ": Case 2";
      // Case 2:
      // Hyperedge e does not contain u. Therefore we  have to connect e to the representative u.
      // This reuses the pin slot of v in e's incidence array (i.e. last_pin_slot!)
      e.hash() -= kahypar::math::hash(v);
      e.hash() += kahypar::math::hash(u);
      _incidence_array[last_pin_slot] = u;
      tmp_incident_nets.push_back(he);
    }
  }

  // ####################### Remove / Restore Hyperedges #######################

  // ! Removes hyperedge e from the incident nets of vertex hn
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void removeIncidentEdgeFromHypernode(const HyperedgeID e,
                                                                       const HypernodeID u) {
    using std::swap;
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");

    IncidentNetVector<HyperedgeID>& incident_nets_of_u = _incident_nets[u];
    size_t pos = 0;
    for ( ; pos < incident_nets_of_u.size(); ++pos ) {
      if ( incident_nets_of_u[pos] == e ) {
        break;
      }
    }
    ASSERT(pos < incident_nets_of_u.size());
    swap(incident_nets_of_u[pos], incident_nets_of_u.back());
    incident_nets_of_u.pop_back();
  }

  // ! Inserts hyperedge he to incident nets array of vertex hn
  KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void connectHypernodeWithIncidentEdge(const HyperedgeID e,
                                                                        const HypernodeID u) {
    using std::swap;
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    IncidentNetVector<HyperedgeID>& incident_nets_of_u = _incident_nets[u];
    HEAVY_REFINEMENT_ASSERT(std::count(incident_nets_of_u.cbegin(), incident_nets_of_u.cend(), e) == 0,
                        "HN" << u << "is already connected to HE" << e);
    incident_nets_of_u.push_back(e);
  }

  // ! Number of hypernodes
  HypernodeID _num_hypernodes;
  // ! Number of removed hypernodes
  HypernodeID _num_removed_hypernodes;
  // ! Number of hyperedges
  HyperedgeID _num_hyperedges;
  // ! Number of removed hyperedges
  HyperedgeID _num_removed_hyperedges;
  // ! Maximum size of a hyperedge
  HypernodeID _max_edge_size;
  // ! Number of pins
  HypernodeID _num_pins;
  // ! Total degree of all vertices
  HypernodeID _total_degree;
  // ! Total weight of hypergraph
  HypernodeWeight _total_weight;

  // ! Hypernodes
  Array<Hypernode> _hypernodes;
  // ! Pins of hyperedges
  IncidentNets _incident_nets;
  // ! Atomic bool vector used to acquire unique ownership of hypernodes
  OwnershipVector _acquired_hns;
  // ! Indicates how many contractions are currently registered on a hypernode
  ReferenceCountVector _hn_ref_count;
  // ! Tracks the contraction tree of the hypergraph
  parallel::scalable_vector<HypernodeID> _contraction_tree;


  // ! Hyperedges
  Array<Hyperedge> _hyperedges;
  // ! Incident nets of hypernodes
  IncidenceArray _incidence_array;
  // ! Atomic bool vector used to acquire unique ownership of hyperedges
  OwnershipVector _acquired_hes;
  // ! Collects hyperedes that will be adjacent to a vertex after a contraction
  ThreadLocalHyperedgeVector _tmp_incident_nets;
  // ! Collects hyperedge contractions that failed due to failed acquired ownership
  ThreadLocalHyperedgeVector _failed_hyperedge_contractions;

  // ! Community Information and Stats
  CommunitySupport<DynamicHypergraph> _community_support;

};

} // namespace ds
} // namespace mt_kahypar