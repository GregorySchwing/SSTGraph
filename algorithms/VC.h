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
// vertex map function to mark visited vertexSubset
template <typename T, typename SM> 
struct VC_Vertex_F {
  T *inCover;
  const SM &G;

  explicit VC_Vertex_F(T *_inCover, const SM &_G) : inCover(_inCover), G(_G) {}
  inline bool operator()(uintE i) {
    inCover[i] = G.getDegree(i) > 0;
    return G.getDegree(i) > 0;
  }
};

template <typename T, typename SM> 
struct VC_Should_Be_Deleted_F {
  T *inCover;
  const SM &G;

  explicit VC_Should_Be_Deleted_F(T *_inCover, const SM &_G) : inCover(_inCover), G(_G) {}
  inline bool operator()(uintE i) {
    inCover[i] = G.getDegree(i) > 0;
    return G.getDegree(i) > 0;
  }
};

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
    if (G.getDegree(s) > G.getDegree(d)) {
      inCover[d] = 0;
      __sync_fetch_and_and(&inCover[d], 0);
      return false;
    } else if (G.getDegree(s) == G.getDegree(d)){
      __sync_fetch_and_and(&inCover[d], s < d);
      return false;
    } else {
      // Start at 1, no need to set
      //inCover[d] &= 1;
      return false;
    }
  }
  // cond function checks if vertex in remaining vertices set has non-zero degree
  inline bool cond(uint32_t d) {
    return G.getDegree(d) > 0;
  }
};
  
template <typename T, typename SM> struct VC_ROUND_2_F {
  T *inCover;
  T *solution;
  T *numToEdgesRemove;
  T *maxNumToEdgesRemove;
  std::tuple<el_t, el_t> *edgesToRemove;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  VC_ROUND_2_F(T *_inCover, T *_solution, std::tuple<el_t, el_t> *_edgesToRemove, 
  T *_numToEdgesRemove, T *_maxNumToEdgesRemove, SM &_G) : 
  inCover(_inCover), solution(_solution), 
  maxNumToEdgesRemove(_maxNumToEdgesRemove),
  numToEdgesRemove(_numToEdgesRemove),
  edgesToRemove(_edgesToRemove),
  G(_G)  {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    if (inCover[d] && G.getDegree(d)>0) {
      solution[d] = 1;
      if (*numToEdgesRemove + 2 < *maxNumToEdgesRemove){
        edgesToRemove[*numToEdgesRemove] = std::tuple<el_t, el_t>{s, d};
        ++(*numToEdgesRemove);
        edgesToRemove[*numToEdgesRemove] = std::tuple<el_t, el_t>{d, s};
        ++(*numToEdgesRemove);
      } else {
        //printf("Failed to remove edge.  Batch size exceeded!");
        //exit(1);
      }
      return true;
      //return false;
    }
    return false;
    //return true;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    if (inCover[s] && G.getDegree(s)>0) {
      solution[s] = 1;
      if (*numToEdgesRemove + 2 < *maxNumToEdgesRemove){
        uint32_t edgeIndex = __sync_add_and_fetch(numToEdgesRemove, 1);
        edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{s, d};
        edgeIndex = __sync_add_and_fetch(numToEdgesRemove, 1);
        edgesToRemove[edgeIndex] = std::tuple<el_t, el_t>{d, s};
      } else {
        //printf("Failed to remove edge.  Batch size exceeded!");
        //exit(1);
      }
      return true;
      //return false;
    }
    return false;
    //return true;
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
  int32_t b_size = 10000; 
  int32_t b_used = 0; 
  std::tuple<el_t, el_t> *edgesToRemove = (std::tuple<el_t, el_t> *)malloc(b_size * sizeof(std::tuple<el_t, el_t>));

  int32_t *inCover = (int32_t *)malloc(n * sizeof(int32_t));
  int32_t *solution = (int32_t *)malloc(n * sizeof(int32_t));
  
  // Need to do batch removes, by virtue of how edges need to exist for 
  // edge map to find the edges to remove.
  //std::vector<std::tuple<el_t, el_t>> es(b_size);
  // MEMSET all to 1
  // The plan is to have all nondegree 0 vertices start in the cover, and then exclude most of them.
  parallel_for(int64_t i = 0; i < n; i++) { inCover[i] = 1; }
  parallel_for(int64_t i = 0; i < n; i++) { solution[i] = 0; }

  if (n == 0) {
    return solution;
  }
  // Dense
  //VertexSubset remaining_vertices =
  //    VertexSubset(0, n, true); // initial set contains all vertices
  // Sparse
  VertexSubset remaining_vertices =
      VertexSubset(0, n, true); // initial set contains all vertices
  VertexSubset vertices_to_delete;
  while (remaining_vertices.non_empty()) { // loop until set of remaining vertices is empty
    // Read phase
    // Set inCover array (I)
    G.edgeMap(remaining_vertices, VC_ROUND_1_F(inCover, G), false, 20);
    //remaining_vertices.print();
    bool continueOn = false;
    do {
      b_used = 0;
      // returns vertices to delete
      vertices_to_delete = G.edgeMap(remaining_vertices, VC_ROUND_2_F(inCover, solution, edgesToRemove, &b_used, &b_size, G), true, 20);
      printf("VERTS 2 DEL\n");
      vertices_to_delete.print();
      printf("non_empty : %d\n", vertices_to_delete.non_empty());
      printf("batch removes : %d\n", b_used);
      printf("inCover\n");
      for(int64_t i = 0; i < n; i++) { printf("%u ", inCover[i]); }
      printf("\n");
      // Write phase
      G.remove_batch(edgesToRemove, b_used);
    } while(vertices_to_delete.non_empty());
    //printf("B4\n");
    //for(int64_t i = 0; i < n; i++) { printf("%u ", G.getDegree(i)); }
    //printf("\n");
    // remove degree zero vertices
    //printf("after\n");
    //for(int64_t i = 0; i < n; i++) { printf("%u ", G.getDegree(i)); }
    //printf("\n");
    VertexSubset nonzero_degree_remaining_vertices = G.vertexMap(remaining_vertices, VC_Vertex_F(inCover, G), true); // mark visited
    remaining_vertices = nonzero_degree_remaining_vertices;
    printf("REM VERTS\n");
    remaining_vertices.print();
    printf("after\n");
    for(int64_t i = 0; i < n; i++) { printf("%u ", G.getDegree(i)); }
    printf("\n");
  }
  remaining_vertices.del();
  free(inCover);
  free(edgesToRemove);
  return solution;
}
