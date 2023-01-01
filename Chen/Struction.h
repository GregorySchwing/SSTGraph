#ifndef STRUCTION_H
#define STRUCTION_H

template <typename T, typename SM> struct SET_ANTI_EDGES_F {
  T *numberAntiEdges;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  SET_ANTI_EDGES_F(T* _numberAntiEdges, SM &_G) : 
  numberAntiEdges(_numberAntiEdges),
  G(_G)  {}
  // Assumes undirected graph
  inline bool update(uint32_t s, uint32_t d) { // Update
    //  Self-edges aren't considered in common neighbors, 
    //  therefore, the source isn't included in the common neighbors
    //  of the destination, so subtract 1.
    numberAntiEdges[d] += (G.getDegree(d) - G.common_neighbors(s,d) - 1);
    //printf("s %u d %u nAE %u degree(%u) : %u common neighbors (%u, %u) : %u\n", s , d, G.getDegree(d) - G.common_neighbors(s,d) - 1,
    //d, G.getDegree(d), s, d, G.common_neighbors(s,d));
    return true;
  }
  // Assumes undirected graph
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    // Assuming I am not a neighbor to myself,
    // thus G.common_neighbors(s,d) > 0 indicates a triangle.
    uint32_t edgeIndex = __sync_fetch_and_add(&numberAntiEdges[d], (G.getDegree(d) - G.common_neighbors(s,d) - 1));
    return true;
  }
  
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};


template <typename T, typename SM> 
struct SET_STRUCTION_F {
  T *numberAntiEdges;
  T *performStruction;
  const SM &G;

  explicit SET_STRUCTION_F(T *_numberAntiEdges, T *_performStruction, const SM &_G) : numberAntiEdges(_numberAntiEdges), performStruction(_performStruction), G(_G) {}
  inline bool operator()(uintE i) {
    return performStruction[i] = (G.getDegree(i) > numberAntiEdges[i]) && numberAntiEdges[i];
  }
};

template <typename T, typename SM> 
struct GET_STRUCTION_SET_F {
  T *performStruction;
  T *numStructionNeighbors;
  const SM &G;

  explicit GET_STRUCTION_SET_F(T *_performStruction, T *_numStructionNeighbors, const SM &_G) : 
  performStruction(_performStruction), 
  numStructionNeighbors(_numStructionNeighbors),
  G(_G) {}
  inline bool operator()(uintE i) {
    return performStruction[i];
  }
};

template <typename T, typename SM> 
struct GET_STRUCTION_SET_AND_NEIGHBORS_F {
  T *performStruction;
  T *numStructionNeighbors;
  const SM &G;

  explicit GET_STRUCTION_SET_AND_NEIGHBORS_F(T *_performStruction, T *_numStructionNeighbors, const SM &_G) : 
  performStruction(_performStruction), 
  numStructionNeighbors(_numStructionNeighbors),
  G(_G) {}
  inline bool operator()(uintE i) {
    return performStruction[i] || numStructionNeighbors[i];
  }
};


template <typename T, typename SM> struct DELETE_NEIGHBORHOOD_OF_STRUCTION_VERTEX_F {
  T *performStruction;
  T *maxVertex;
  T *numStructionNeighbors;
  T *numToEdgesRemove;
  T *maxNumToEdgesRemove;
  std::tuple<el_t, el_t> *edgesToRemove;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  DELETE_NEIGHBORHOOD_OF_STRUCTION_VERTEX_F(T *_performStruction, T *_maxVertex, T *_numStructionNeighbors, std::tuple<el_t, el_t> *_edgesToRemove, 
  T *_numToEdgesRemove, T *_maxNumToEdgesRemove, SM &_G) : 
  performStruction(_performStruction),
  maxVertex(_maxVertex),
  numStructionNeighbors(_numStructionNeighbors),
  maxNumToEdgesRemove(_maxNumToEdgesRemove),
  numToEdgesRemove(_numToEdgesRemove),
  edgesToRemove(_edgesToRemove),
  G(_G)  {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    bool deletedEdge = false;
    if (performStruction[d]){
      if (G.has(s,d)){
        uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
        if (edgeIndex < *maxNumToEdgesRemove) 
          edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{s, d};
        deletedEdge = true;
      }
      if (G.has(d,s)){
        uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
        if (edgeIndex < *maxNumToEdgesRemove) 
          edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{d, s};
        deletedEdge = true;
      }
    } else if (numStructionNeighbors[d]){
        bool deletedEdge = false;
        // source has an edge to the struction v0
        if (G.has(s,maxVertex[d])){
          uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
          if (edgeIndex < *maxNumToEdgesRemove) 
            edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{s, d};
          deletedEdge = true;
        }
        if (G.has(maxVertex[d],s)){
          uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
          if (edgeIndex < *maxNumToEdgesRemove) 
            edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{d, s};
          deletedEdge = true;
        }
    }
    return deletedEdge;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    bool deletedEdge = false;
    if (performStruction[s]){
      if (G.has(s,d)){
        uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
        if (edgeIndex < *maxNumToEdgesRemove) 
          edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{s, d};
        deletedEdge = true;
      }
      if (G.has(d,s)){
        uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
        if (edgeIndex < *maxNumToEdgesRemove) 
          edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{d, s};
        deletedEdge = true;
      }
    } else if (numStructionNeighbors[s]){
        bool deletedEdge = false;
        // source has an edge to the struction v0
        if (G.has(d,maxVertex[s])){
          uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
          if (edgeIndex < *maxNumToEdgesRemove) 
            edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{s, d};
          deletedEdge = true;
        }
        if (G.has(maxVertex[s],d)){
          uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
          if (edgeIndex < *maxNumToEdgesRemove) 
            edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{d, s};
          deletedEdge = true;
        }
    }
    return deletedEdge;
  }
  
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};


template <typename T, typename SM> struct DELETE_ALL_VERTICES_IN_VERTEX_SUBSET_F {
  T *numToEdgesRemove;
  T *maxNumToEdgesRemove;
  std::tuple<el_t, el_t> *edgesToRemove;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  DELETE_ALL_VERTICES_IN_VERTEX_SUBSET_F(std::tuple<el_t, el_t> *_edgesToRemove, 
  T *_numToEdgesRemove, T *_maxNumToEdgesRemove, SM &_G) : 
  maxNumToEdgesRemove(_maxNumToEdgesRemove),
  numToEdgesRemove(_numToEdgesRemove),
  edgesToRemove(_edgesToRemove),
  G(_G)  {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    bool deletedEdge = false;
    if (G.has(s,d)){
      uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
      if (edgeIndex < *maxNumToEdgesRemove) 
        edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{s, d};
      deletedEdge = true;
    }
    if (G.has(d,s)){
      uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
      if (edgeIndex < *maxNumToEdgesRemove) 
        edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{d, s};
      deletedEdge = true;
    }
    return deletedEdge;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    bool deletedEdge = false;
    if (G.has(s,d)){
      uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
      if (edgeIndex < *maxNumToEdgesRemove) 
        edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{s, d};
      deletedEdge = true;
    }
    if (G.has(d,s)){
      uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
      if (edgeIndex < *maxNumToEdgesRemove) 
        edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{d, s};
      deletedEdge = true;
    }
    return deletedEdge;
  }
  
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};


template <typename T, typename SM> struct STRUCTION_OP_F {
  T *performStruction;
  T *maxVertex;
  T *numStructionNeighbors;
  T *numToEdgesRemove;
  T *maxNumToEdgesRemove;
  std::tuple<el_t, el_t> *edgesToRemove;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  STRUCTION_OP_F(T *_performStruction, T *_maxVertex, T *_numStructionNeighbors, std::tuple<el_t, el_t> *_edgesToRemove, 
  T *_numToEdgesRemove, T *_maxNumToEdgesRemove, SM &_G) : 
  performStruction(_performStruction),
  maxVertex(_maxVertex),
  numStructionNeighbors(_numStructionNeighbors),
  maxNumToEdgesRemove(_maxNumToEdgesRemove),
  numToEdgesRemove(_numToEdgesRemove),
  edgesToRemove(_edgesToRemove),
  G(_G)  {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    bool deletedEdge = false;
    if (performStruction[d]){
      if (G.has(s,d)){
        uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
        if (edgeIndex < *maxNumToEdgesRemove) 
          edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{s, d};
        deletedEdge = true;
      }
      if (G.has(d,s)){
        uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
        if (edgeIndex < *maxNumToEdgesRemove) 
          edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{d, s};
        deletedEdge = true;
      }
    } else if (numStructionNeighbors[d]){
        bool deletedEdge = false;
        // source has an edge to the struction v0
        if (G.has(s,maxVertex[d])){
          uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
          if (edgeIndex < *maxNumToEdgesRemove) 
            edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{s, d};
          deletedEdge = true;
        }
        if (G.has(maxVertex[d],s)){
          uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
          if (edgeIndex < *maxNumToEdgesRemove) 
            edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{d, s};
          deletedEdge = true;
        }
    }
    return deletedEdge;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    bool deletedEdge = false;
    if (performStruction[s]){
      if (G.has(s,d)){
        uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
        if (edgeIndex < *maxNumToEdgesRemove) 
          edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{s, d};
        deletedEdge = true;
      }
      if (G.has(d,s)){
        uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
        if (edgeIndex < *maxNumToEdgesRemove) 
          edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{d, s};
        deletedEdge = true;
      }
    } else if (numStructionNeighbors[s]){
        bool deletedEdge = false;
        // source has an edge to the struction v0
        if (G.has(d,maxVertex[s])){
          uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
          if (edgeIndex < *maxNumToEdgesRemove) 
            edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{s, d};
          deletedEdge = true;
        }
        if (G.has(maxVertex[s],d)){
          uint32_t edgeIndex = __sync_fetch_and_add(numToEdgesRemove, 1);
          if (edgeIndex < *maxNumToEdgesRemove) 
            edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{d, s};
          deletedEdge = true;
        }
    }
    return deletedEdge;
  }
  
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};

template <typename T, typename SM> struct SOLVE_MIS_F {
  T *performStruction;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  const SM &G;
  SOLVE_MIS_F(T *_performStruction, const SM &_G) : performStruction(_performStruction), G(_G)  {}

  /*
    In dense mode, EdgeMap loops over all vertices and
    looks at incoming edges to see if the source is part of the
    vertex set. This does not require locking because each vertex
    only updates itself and is preferred when the vertex set is large.

    I am the destination.  I only update inCover[d].
    Degrees are const in this class.

  */

  inline bool update(uint32_t s, uint32_t d) { // Update
    if (performStruction[s] > performStruction[d]) {
      performStruction[d] = 0;
      return false;
    } else if (performStruction[s] == performStruction[d]){
      // Arbitrary tiebreaker, flipped to get struction_a example to work.
      //performStruction[d] &= h(s) < h(d);
      performStruction[d] &= h(s) > h(d);
      return false;
    } else {
      // Start at 1, no need to set
      //inCover[d] &= 1;
      return false;
    }
  }
  /*
    In sparse mode, EdgeMap
    iterates over the outgoing edges of each vertex in the subset
    and updates the destination vertex for each edge. Because it_i is
    run in parallel, synchronization must be used when accessing
    the destination vertex data.

    I am the source.  I update inCover[d] using synchronization.
    Degrees are const in this class.

  */
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    if (performStruction[s] > performStruction[d]) {
      //inCover[d] = 0;
      __sync_fetch_and_and(&performStruction[d], 0);
      return false;
    } else if (performStruction[s] == performStruction[d]){
      __sync_fetch_and_and(&performStruction[d], h(s) < h(d));
      return false;
    } else {
      // Start at 1, no need to set
      //inCover[d] &= 1;
      return false;
    }
  }
  // cond function checks if vertex in remaining vertices set has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
    //return G.getDegree(d) > 0;
  }
};


template <typename T, typename SM> struct SET_NUM_STRUCTION_NEIGHBORS_F {
  T *performStruction;
  T *numStructionNeighbors;

  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  const SM &G;
  SET_NUM_STRUCTION_NEIGHBORS_F(T *_performStruction, T *_numStructionNeighbors, const SM &_G) : 
  performStruction(_performStruction), 
  numStructionNeighbors(_numStructionNeighbors),
  G(_G)  {}

  /*
    In dense mode, EdgeMap loops over all vertices and
    looks at incoming edges to see if the source is part of the
    vertex set. This does not require locking because each vertex
    only updates itself and is preferred when the vertex set is large.

    I am the destination.  I only update inCover[d].
    Degrees are const in this class.

  */

  inline bool update(uint32_t s, uint32_t d) { // Update
    if (performStruction[s])
      numStructionNeighbors[d] += performStruction[s];
    return false;
  }
  /*
    In sparse mode, EdgeMap
    iterates over the outgoing edges of each vertex in the subset
    and updates the destination vertex for each edge. Because it_i is
    run in parallel, synchronization must be used when accessing
    the destination vertex data.

    I am the source.  I update inCover[d] using synchronization.
    Degrees are const in this class.

  */
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    if (performStruction[s]){
      __sync_fetch_and_and(&numStructionNeighbors[d], performStruction[s]);
    }
    return false;
  }
  // cond function checks if vertex in remaining vertices set has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
    //return G.getDegree(d) > 0;
  }
};


template <typename T, typename SM> struct SET_LARGEST_VERTEX_STRUCT_F {
  T *maxVertex;
  T *numStructionNeighbors;
  T *performStruction;
   // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  SET_LARGEST_VERTEX_STRUCT_F(T* _maxVertex, T* _numStructionNeighbors, T* _performStruction, SM &_G) : 
  maxVertex(_maxVertex),
  numStructionNeighbors(_numStructionNeighbors),
  performStruction(_performStruction),
  G(_G)  {}
  // Only set this for non-struction v_0, looking for examples 
  // where non-struction v_0 is neighbors with two struction candidates.
  inline bool update(uint32_t s, uint32_t d) { // Update
    if (numStructionNeighbors[d] > 1){
      if (performStruction[s] && h(s) > h(maxVertex[d])){
        maxVertex[d] = s;
      }
    }
    return true;
  }
  // May need to use the write_min and cmp and swap formalism from Bellman-Ford
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    if (numStructionNeighbors[d] > 1){
      if (performStruction[s] && h(s) > h(maxVertex[d])){
        maxVertex[d] = s;
      }
    }
  }
  
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};


template <typename T, typename SM> struct RESOLVE_CONFLICTS_STRUCT_F {
  T *maxVertex;
  T *numStructionNeighbors;
  T *performStruction;

   // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  RESOLVE_CONFLICTS_STRUCT_F(T* _maxVertex, T* _numStructionNeighbors, T* _performStruction, SM &_G) : 
  maxVertex(_maxVertex),
  numStructionNeighbors(_numStructionNeighbors),
  performStruction(_performStruction),

  G(_G)  {}
  // Only set this for non-struction v_0, looking for examples 
  // where non-struction v_0 is neighbors with two struction candidates.
  inline bool update(uint32_t s, uint32_t d) { // Update
    if (numStructionNeighbors[s] > 1 && maxVertex[s] != d){
      performStruction[d] = 0;
    }
    return false;
  }
  // Assumes undirected graph
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    if (numStructionNeighbors[s] > 1 && maxVertex[s] != d){
      performStruction[d] = 0;
    }
    return false;
  }
  
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};

#endif