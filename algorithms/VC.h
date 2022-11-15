#pragma once
#include "../SparseMatrix.hpp"
// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
template <typename T, typename SM> struct VC_ROUND_1_F {
  T *inCover;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  const SM &G;
  VC_ROUND_1_F(T *_inCover, const SM &_G) : inCover(_inCover), G(_G)  {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    if (G.getDegree(s) > G.getDegree(d)) {
      inCover[d] = 0;
      return false;
    } else if (G.getDegree(s) == G.getDegree(d)){
      inCover[d] &= s < d;
      return false;
    } else {
      // Start at 1, no need to set
      //inCover[d] &= 1;
      return false;
    }
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    return __sync_bool_compare_and_swap(&inCover[d], -1, s);
  }
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    if (G.getDegree(d) > 0){
      return true;
    }
    return false;
  }
};
  
template <typename T, typename SM> struct VC_ROUND_2_F {
  T *inCover;
  T *solution;

  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  VC_ROUND_2_F(T *_inCover, T *_solution, SM &_G) : inCover(_inCover), solution(_solution), G(_G)  {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    if (inCover[d]) {
      solution[d] = 1;
      G.remove(d,s);
      return false;
    } else if (inCover[s]){
      G.remove(s,d);
    } 
    return true;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    return __sync_bool_compare_and_swap(&inCover[d], -1, s);
  }
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};


template <typename SM> int32_t *VC_with_edge_map(SM &G) {
  int64_t n = G.get_rows();
  if (n == 0) {
    return nullptr;
  }
  // creates inCover array, initialized to all -1, except for start
  int32_t *inCover = (int32_t *)malloc(n * sizeof(int32_t));
  int32_t *solution = (int32_t *)malloc(n * sizeof(int32_t));
  // MEMSET all to 1
  // The plan is to have all nondegree 0 vertices start in the cover, and then exclude most of them.
  parallel_for(int64_t i = 0; i < n; i++) { inCover[i] = 1; }
  parallel_for(int64_t i = 0; i < n; i++) { solution[i] = 0; }

  if (n == 0) {
    return solution;
  }
  VertexSubset remaining_vertices =
      VertexSubset(0, n, true); // initial set contains all vertices
  while (remaining_vertices.non_empty()) { // loop until set of remaining vertices is empty
    parallel_for(int64_t i = 0; i < n; i++) { inCover[i] = G.getDegree(i) > 0; }
    G.edgeMap(remaining_vertices, VC_ROUND_1_F(inCover, G), true, 20);
    remaining_vertices = G.edgeMap(remaining_vertices, VC_ROUND_2_F(inCover, solution, G), true, 20);
    remaining_vertices.print();
  }
  remaining_vertices.del();
  free(inCover);
  return solution;
}
