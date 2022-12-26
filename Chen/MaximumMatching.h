
#pragma once
#include "../SparseMatrix.hpp"
#include "Match.h"
#include <limits>

template <typename T, typename SM> 
struct GET_LEVEL_I_F {
  int64_t i;
  T *evenlevel;
  T *oddlevel;
  const SM &G;

  explicit GET_LEVEL_I_F(int64_t _i, T *_evenlevel, T *_oddlevel, const SM &_G) : i(_i), evenlevel(_evenlevel), oddlevel(_oddlevel), G(_G) {}
  inline bool operator()(uintE v) {
    return min(evenlevel[v], oddlevel[v]) == i;
  }
};



template <typename T, typename SM> struct EVEN_LEVEL_F {
  T * evenlevel;
  T * oddlevel;
  T * bridges;
  T * predecessors;
  T * anomalies;
  bool *unvisited_v;
  bool *unused_e;
  bool *unvisited_e;
  int64_t i;
  const SM &G;
  EVEN_LEVEL_F(T * _evenlevel,
        T * _oddlevel,
        T * _bridges,
        T * _predecessors,
        T * _anomalies,
        bool *_unvisited_v,
        bool *_unused_e,
        bool *_unvisited_e,
        int64_t _i,
        const SM &_G) :
  evenlevel(_evenlevel),
  oddlevel(_oddlevel),
  bridges(_bridges),
  predecessors(_predecessors),
  anomalies(_anomalies),
  unvisited_v(_unvisited_v),
  unused_e(_unused_e),
  unvisited_e(_unvisited_e),
  i(_i),
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
    if (evenlevel[s] < std::numeric_limits<int32_t>::max()){
      int64_t temp = (evenlevel[s] + evenlevel[d])/2;
      //bridges[temp] = ???
    } else {
      if (oddlevel[s] == std::numeric_limits<int32_t>::max()){
        oddlevel[s] = i + 1;
      }
      if (oddlevel[s] == i + 1){
        //predecessors[s] = ???
      }
      if (oddlevel[s] < i){
        //anomalies[s] = ???
      }
    } 
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
    //if (performStruction[s]){
    //  __sync_fetch_and_and(&numStructionNeighbors[d], performStruction[s]);
    //}
    return false;
  }
  // cond function checks if vertex in remaining vertices set has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
    //return G.getDegree(d) > 0;
  }
};


template <typename SM> 
class MaximumMatcher {
  public:
    MaximumMatcher(
                    SM &_G,
                    VertexSubset & _remainingVertices) : 
                    G(_G),
                    n(_G.get_rows()),
                    e(_G.get_cols()),
                    remainingVertices(_remainingVertices)
{
  evenlevel = (int32_t *)malloc(n * sizeof(int32_t));
  oddlevel = (int32_t *)malloc(n * sizeof(int32_t));
  blossom = (int32_t *)malloc(n * sizeof(int32_t));
  predecessors = (int32_t *)malloc(n * sizeof(int32_t));
  anomalies = (int32_t *)malloc(n * sizeof(int32_t));
  unvisited_v = (bool *)malloc(n * sizeof(bool));

  unused_e = (bool *)malloc(e * sizeof(bool));
  unvisited_e = (bool *)malloc(e * sizeof(bool));

  bridges = (int32_t *)malloc(n * sizeof(int32_t));

  request = (int32_t *)malloc(n * sizeof(int32_t));
  match = (int32_t *)malloc(n * sizeof(int32_t));

  m.Match(G, remainingVertices, request, match);

  // 0
  for(int64_t i = 0; i < n; i++) { 
    evenlevel[i] = std::numeric_limits<int32_t>::max();
    oddlevel[i] = std::numeric_limits<int32_t>::max();
    blossom[i] = -1;
    predecessors[i] = -1;
    anomalies[i] = -1;
    unvisited_v[i] = true;
  }
  for(int64_t i = 0; i < e; i++) { 
    unused_e[i] = true;
    unvisited_e[i] = true;
  }
  i = -1;

}

    int Search() {
      // 1
      for(int64_t i = 0; i < n; i++)
        if (match[i] < 4) evenlevel[i] = 0;

      // 2 
      VertexSubset level_i_vertices = G.vertexMap(remainingVertices, GET_LEVEL_I_F(i, evenlevel, oddlevel, G), true); // mark visited
      if (!level_i_vertices.non_empty())
        return 0;

      // 3 
      if (i % 2 == 0){

      }
    }

    private:
        SM &G;
        VertexSubset & remainingVertices;
        Matcher m;

        // |V| = n
        int64_t n;
        int64_t e;
        int64_t i;

        // Length V
        int32_t *evenlevel;
        int32_t *oddlevel;
        int32_t *blossom;
        int32_t *predecessors;
        int32_t *anomalies;
        bool *unvisited_v;

        // Length E
        bool *unused_e;
        bool *unvisited_e;

        // Dynamic length
        int32_t *bridges;

        // Matcher variables
        int32_t *request;
        int32_t *match;

};
