
#pragma once
#include "../SparseMatrix.hpp"
#include "Match.h"
#include <limits>

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
      for(int64_t i = 0; i < n; i++) { 
        if (match[i] < 4) evenlevel[i] = 0;
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
