
template <typename T, typename SM> struct SET_DOMINATED_F {
  T *dominates;
  SM &G;
  SET_DOMINATED_F( T *_dominates, SM &_G) : 
  dominates(_dominates),
  G(_G)  {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    // d dominates s; 
    // && (G.getDegree(d)!=1)) - avoid leaves, since they technically dominate
    dominates[d] |= ((G.common_neighbors(s,d) - G.getDegree(s) - 1) == 0) && (G.getDegree(d)!=1);
    return dominates[d];
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    // && (G.getDegree(d)!=1)) - avoid leaves, since they technically dominate
    uint32_t vertexDominates = __sync_or_and_fetch(&dominates[d], ((G.common_neighbors(s,d) - G.getDegree(s) - 1) == 0) && (G.getDegree(d)!=1));
    return true;
  }
  
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};
