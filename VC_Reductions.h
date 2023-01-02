#pragma once
#include "SparseMatrix.hpp"
#include "ProgressBar.h"
#include <map>
#include <iterator>
#include <algorithm> // std::set_union, std::sort
#include "Chen/Struction.h"
#include "Chen/Dominated.h"
#include "Chen/GeneralFold.h"
#include "Chen/Match.h"
#include "Chen/MaximumMatching.h"
#include "Chen/CrownReduction.h"

#include "utils/CrownReduction.h"
#include "utils/Struction.h"
#include "utils/Dominated.h"
#include "utils/MaxDegree.h"

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
struct VC_Red_Vertex_F {
  T *inCover;
  const SM &G;

  explicit VC_Red_Vertex_F(T *_inCover, const SM &_G) : inCover(_inCover), G(_G) {}
  inline bool operator()(uintE i) {
    inCover[i] = G.getDegree(i) > 0;
    return G.getDegree(i) > 0;
  }
};



template <typename T, typename SM> 
struct SET_LEAVES_F {
  T *inCover;
  const SM &G;

  explicit SET_LEAVES_F(T *_isLeaf, const SM &_G) : inCover(_isLeaf), G(_G) {}
  inline bool operator()(uintE i) {
    inCover[i] = G.getDegree(i) == 1;
    return inCover[i];
  }
};

template <typename T, typename SM> 
struct SET_POSSIBLE_TRIANGLES_F {
  T *isTriangle;
  const SM &G;

  explicit SET_POSSIBLE_TRIANGLES_F(T *_isTriangle, const SM &_G) : isTriangle(_isTriangle), G(_G) {}
  inline bool operator()(uintE i) {
    isTriangle[i] = G.getDegree(i) == 2;
    return isTriangle[i];
  }
};

template <typename T, typename SM> 
struct SET_LEAVES_1_F {
  T *isLeaf;
  const SM &G;

  explicit SET_LEAVES_1_F(T *_isLeaf, const SM &_G) : isLeaf(_isLeaf), G(_G) {}
  inline bool operator()(uintE i) {
    isLeaf[i] = G.getDegree(i) == 1;
    return isLeaf[i];
  }
};

struct BF_LEAVES_F {
  int *isLeaf;
  int *inCover;
  BF_LEAVES_F(intE *_isLeaf, int *_inCover)
      : isLeaf(_isLeaf), inCover(_inCover) {}
  // Update ShortestPathLen if found a shorter path
  inline bool update(uintE s, uintE d) {
    if (isLeaf[s] > isLeaf[d]) {
      if (inCover[d] == 0) {
        //printf("DELETE LEAF %d %d\n",s,d);
        inCover[d] = 1;
        return 1;
      }
    } else if (isLeaf[s] && isLeaf[d] && h(s) < h(d)){
        if (inCover[d] == 0) {
        //printf("DELETE LEAF %d %d\n",s,d);
        inCover[d] = 1;
        return 1;
      }
    }
    return 0;
  }
  inline bool updateAtomic(uintE s, uintE d) { // atomic Update
    //intE newDist = isLeaf[s] + edgeLen;
    return (CAS(&inCover[d], 0, 1));
  }
  inline bool cond([[maybe_unused]] uintE s) { return true; }
};

struct BF_TRIANGLES_F {
  int *isTri;
  int *inCover;
  BF_TRIANGLES_F(intE *_isTri, int *_inCover)
      : isTri(_isTri), inCover(_inCover) {}
  // Update ShortestPathLen if found a shorter path
  inline bool update(uintE s, uintE d) {
    if (isTri[s] > isTri[d]) {
      if (inCover[d] == 0) {
        //printf("DELETE LEAF %d %d\n",s,d);
        inCover[d] = 1;
        return 1;
      }
    } else if (isTri[s] && isTri[d] && h(s) < h(d)){
        if (inCover[d] == 0) {
        //printf("DELETE LEAF %d %d\n",s,d);
        inCover[d] = 1;
        return 1;
      }
    }
    return 0;
  }
  inline bool updateAtomic(uintE s, uintE d) { // atomic Update
    //intE newDist = isLeaf[s] + edgeLen;
    return (CAS(&inCover[d], 0, 1));
  }
  inline bool cond([[maybe_unused]] uintE s) { return true; }
};


template <typename T, typename SM> struct SET_LEAVES_2_F {
  T *inCover;
  const T *isLeaf;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  const SM &G;
  SET_LEAVES_2_F(const T *_isLeaf, T *_inCover, const SM &_G) : isLeaf(_isLeaf), inCover(_inCover), G(_G)  {}

  /*
    In dense mode, EdgeMap loops over all vertices and
    looks at incoming edges to see if the source is part of the
    vertex set. This does not require locking because each vertex
    only updates itself and is preferred when the vertex set is large.

    I am the destination.  I only update inCover[d].
    Degrees are const in this class.

  */

  inline bool update(uint32_t s, uint32_t d) { // Update
    int bothLeaves = isLeaf[s] && isLeaf[d];
    if (bothLeaves)
      printf("isLeaf[%d] %d isLeaf[%d] %d\n",s,isLeaf[s],d,isLeaf[d]);

    if (isLeaf[s] > isLeaf[d]) {
      printf("isLeaf[%d] %d isLeaf[%d] %d\n",s,isLeaf[s],d,isLeaf[d]);
      inCover[d] |= 1;
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
    if (G.getDegree(s) == 1){
      if (G.getDegree(s) == G.getDegree(d)){
        __sync_fetch_and_and(&inCover[d], h(s) < h(d));
        return false;
      } else {
        __sync_fetch_and_or(&inCover[d], 1);
        return false;
      }
    }
  }
  // cond function checks if vertex in remaining vertices set has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
    //return G.getDegree(d) > 0;
  }
};


template <typename T, typename SM> struct LEAF_REDUCTION_RULE_F {
  T *inCover;
  T *solution;
  T *numToEdgesRemove;
  T *maxNumToEdgesRemove;
  std::tuple<el_t, el_t> *edgesToRemove;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  LEAF_REDUCTION_RULE_F( T *_isLeaf, T *_solution, std::tuple<el_t, el_t> *_edgesToRemove, 
  T *_numToEdgesRemove, T *_maxNumToEdgesRemove, SM &_G) : 
  inCover(_isLeaf),
  solution(_solution),
  maxNumToEdgesRemove(_maxNumToEdgesRemove),
  numToEdgesRemove(_numToEdgesRemove),
  edgesToRemove(_edgesToRemove),
  G(_G)  {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    bool deletedEdge = false;
    // You're a leaf
    if (inCover[s]){
      //printf("LEAF %u %"PRIu32"\n", d, G.getDegree(d));
      // and I'm a leaf
      if (inCover[d]){
        //printf("DOUBLE LEAF %u %u (%"PRIu32") (%"PRIu32")\n", s, d, G.getDegree(s), G.getDegree(d));
        // tiebreak
        if (h(s) < h(d)){
          //printf("SINCE DOUBLE LEAF %u %u (%"PRIu32") (%"PRIu32") only add %u\n", s, d, G.getDegree(s), G.getDegree(d), d);
          //printf("DELETE LEAF %d %d (%d) (%d)\n",d, s, G.getDegree(d), G.getDegree(s));

          solution[s] = 1;
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
        }          
      } else {
        solution[d] = 1;
        //printf("DELETE LEAF %d %d (%d) (%d)\n",s,d, G.getDegree(s), G.getDegree(d));
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
      }
    }

    return deletedEdge;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    printf("USED ATOMIC\n");
    bool deletedEdge = false;
    if (inCover[s]){
      solution[d] = 1;
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
    }
    return deletedEdge;
  }
  
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};


template <typename T, typename SM> 
struct SET_MAX_VERTEX_F {
  T *inCover;
  T maxV;
  const SM &G;

  explicit SET_MAX_VERTEX_F(T *_inCover, T _maxV, const SM &_G) : inCover(_inCover), maxV(_maxV), G(_G) {}
  inline bool operator()(uintE i) {
    inCover[i] = i == maxV;
    //if (i == maxV)
    //  printf("DELETE MAX %d (%d)\n",i, G.getDegree(i));
    return inCover[i];
  }
};

template <typename T, typename SM> struct DELETE_VERTEX_F {
  T *inCover;
  T *solution;
  T *numToEdgesRemove;
  T *maxNumToEdgesRemove;
  std::tuple<el_t, el_t> *edgesToRemove;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  DELETE_VERTEX_F( T *_isLeaf, T *_solution, std::tuple<el_t, el_t> *_edgesToRemove, 
  T *_numToEdgesRemove, T *_maxNumToEdgesRemove, SM &_G) : 
  inCover(_isLeaf),
  solution(_solution),
  maxNumToEdgesRemove(_maxNumToEdgesRemove),
  numToEdgesRemove(_numToEdgesRemove),
  edgesToRemove(_edgesToRemove),
  G(_G)  {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    bool deletedEdge = false;
    if (inCover[d]){
      solution[d] = 1;
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
    }
    return deletedEdge;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    bool deletedEdge = false;
    if (inCover[s]){
      solution[s] = 1;
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
    }
    return deletedEdge;
  }
  
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};

template <typename T, typename SM> struct SET_TRIANGLES_F {
  T *isTriangle;
  T *solution;
  T *numToEdgesRemove;
  T *maxNumToEdgesRemove;
  std::tuple<el_t, el_t> *edgesToRemove;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  SET_TRIANGLES_F( T *_isTriangle, T* _solution, std::tuple<el_t, el_t> *_edgesToRemove, 
  T *_numToEdgesRemove, T *_maxNumToEdgesRemove, SM &_G) : 
  isTriangle(_isTriangle),
  solution(_solution),
  maxNumToEdgesRemove(_maxNumToEdgesRemove),
  numToEdgesRemove(_numToEdgesRemove),
  edgesToRemove(_edgesToRemove),
  G(_G)  {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    // Assuming I am not a neighbor to myself,
    // thus G.common_neighbors(s,d) > 0 indicates a triangle.
    isTriangle[d] = G.getDegree(d) == 2 && G.common_neighbors(s,d);
    return isTriangle[d];
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    // Assuming I am not a neighbor to myself,
    // thus G.common_neighbors(s,d) > 0 indicates a triangle.
    isTriangle[s] = G.getDegree(s) == 2 && G.common_neighbors(s,d);
    return isTriangle[s];
  }
  
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};


template <typename T, typename SM> struct TRIANGLE_REDUCTION_RULE_F {
  T *isTriangle;
  T *solution;
  T *numToEdgesRemove;
  T *maxNumToEdgesRemove;
  std::tuple<el_t, el_t> *edgesToRemove;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  TRIANGLE_REDUCTION_RULE_F( T *_isTriangle, T* _solution, std::tuple<el_t, el_t> *_edgesToRemove, 
  T *_numToEdgesRemove, T *_maxNumToEdgesRemove, SM &_G) : 
  isTriangle(_isTriangle),
  solution(_solution),
  maxNumToEdgesRemove(_maxNumToEdgesRemove),
  numToEdgesRemove(_numToEdgesRemove),
  edgesToRemove(_edgesToRemove),
  G(_G)  {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    bool deletedEdge = false;
    // Assuming I am not a neighbor to myself,
    // thus G.common_neighbors(s,d) > 0 indicates a triangle.
    if (isTriangle[d]){
      solution[s] = 1;
      printf("DELETE TRI %d\n",s);
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
    }
    return deletedEdge;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    bool deletedEdge = false;
    if (isTriangle[s]){
      printf("DELETE TRI %d\n",d);
      solution[d] = 1;
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
    }
    return deletedEdge;
  }
  
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};

class VC_Reductions {
  public:
    template <typename SM> int32_t* RemoveMaxApproximateMVC(SM &G);
    template <typename SM> int32_t* ChenRemoveMaxApproximateMVC(SM &G);
    //template <typename SM> int32_t* Struction(SM &G);
      /*

    template <typename SM> bool Dominated(SM &approxGraph,
                                          VertexSubset &remaining_vertices,
                                          int32_t *vertexDominates,
                                          int32_t *solution,
                                          int32_t b_size,
                                          std::tuple<el_t, el_t> *edgesToRemove,
                                          int32_t &removeCounter);
    template <typename SM> bool Struction(SM &approxGraph,
                                          VertexSubset &remaining_vertices,
                                          int32_t *numberAntiEdges,
                                          int32_t *performStruction,
                                          int32_t *maxVertex,
                                          int32_t *numStructionNeighbors,
                                          int32_t b_size,
                                          std::tuple<el_t, el_t> *edgesToRemove,
                                          std::tuple<el_t, el_t> *edgesToInsert,
                                          int32_t &removeCounter,
                                          int32_t &insertCounter);
    */

    template <typename SM> bool GeneralFold(SM &approxGraph,
                                          VertexSubset &remaining_vertices,
                                          int32_t *request,
                                          int32_t *match,
                                          int32_t *auxMatch,
                                          int32_t *numStructionNeighbors,
                                          int32_t b_size,
                                          std::tuple<el_t, el_t> *edgesToRemove,
                                          std::tuple<el_t, el_t> *edgesToInsert,
                                          int32_t &removeCounter,
                                          int32_t &insertCounter);

template <typename SM> bool Match(SM &approxGraph,
                                  VertexSubset &remaining_vertices,
                                  int32_t *request,
                                  int32_t *match,
                                  int32_t *maxDegree);

template <typename SM> bool AuxilliaryMatch(SM &approxGraph,
                                  VertexSubset &remaining_vertices,
                                  int32_t *request,
                                  int32_t *match,
                                  int32_t *auxMatch,
                                  int32_t *maxDegree);

template <typename SM> bool FindCrown(SM &approxGraph,
                                      VertexSubset &remaining_vertices,
                                      int32_t *request,
                                      int32_t *match,
                                      int32_t *auxMatch,
                                      int32_t *H_n);
};




//template <typename SM> int32_t VC_Reductions::RemoveMaxApproximateMVC(SM &G){
template <typename SM> int32_t* VC_Reductions::ChenRemoveMaxApproximateMVC(SM &G){
  SparseMatrixV<true, bool> approxGraph(G);
  int64_t n = approxGraph.get_rows(); 
  int64_t e = approxGraph.get_cols(); 

  int32_t *solution = (int32_t *)malloc(n * sizeof(int32_t));

  parallel_for(int64_t i = 0; i < n; i++) { solution[i] = 0; }

  int32_t b_size = 10000; 
  std::tuple<el_t, el_t> *edgesToRemove = (std::tuple<el_t, el_t> *)malloc(b_size * sizeof(std::tuple<el_t, el_t>));
  std::tuple<el_t, el_t> *edgesToInsert = (std::tuple<el_t, el_t> *)malloc(b_size * sizeof(std::tuple<el_t, el_t>));
  int32_t removeCounter = 0;
  int32_t insertCounter = 0;
  bool vertexChanged = false;
  bool foundCrown = false;
  bool foundDominating = false;
  bool foundStruction = false;
	bool hasEdges = true;
  VertexSubset remaining_vertices = VertexSubset(0, n, true); // initial set contains all vertices

  MaximumMatcherBlossom mmb(approxGraph,
                      remaining_vertices);
  CrownReduction cr(approxGraph,
                    remaining_vertices,
                    mmb.get_match(),
                    solution);
  Struction struction(approxGraph,
                    remaining_vertices,
                    b_size,
                    edgesToRemove,
                    edgesToInsert);
  Dominated dom(approxGraph,
                    remaining_vertices,
                    b_size,
                    edgesToRemove,
                    solution);
  MaxDegree md(approxGraph,
                    remaining_vertices,
                    b_size,
                    edgesToRemove,
                    solution,
                    hasEdges);

	while (hasEdges)
	{
    // Reduce as much as possible.
    do {
      vertexChanged = false;
      // Must be called before each FindCrown.
      // It makes sense to make CrownReduction own MMB.
      mmb.edmonds();
      foundCrown = cr.FindCrown();
      foundStruction = struction.FindStruction();
      foundDominating = dom.FindDominated();

      printf("foundCrown %s\n", foundCrown ? "true" : "false");
      printf("foundDominating %s\n", foundDominating ? "true" : "false");
      printf("foundStruction %s\n", foundStruction ? "true" : "false");
      vertexChanged = foundCrown || foundDominating || foundStruction;
    } while (vertexChanged);
  
    md.FindMaxDegree();
  }
  // This is the AbuKhzam CR with a shoddy matching.
  /*
  GeneralFold(approxGraph,
            remaining_vertices,
            numberAntiEdges,
            performStruction,
            maxVertex,
            numStructionNeighbors,
            b_size,
            edgesToRemove,
            edgesToInsert,
            removeCounter,
            insertCounter);
  */

  printf("Chen: ");
  for (uint32_t j = 0; j < approxGraph.get_rows(); j++) {
    if (solution[j])
      printf("%lu ", j);
  }
  printf("\n");
  int32_t vc_count = 0;
  for (int j = 0; j < approxGraph.get_rows(); j++) {
    if(solution[j])
      vc_count += 1;
  }
  printf("Chen solution size: %u\n", vc_count);

  free(edgesToRemove);
  free(edgesToInsert);

  return solution;
}

/*
template <typename SM> bool VC_Reductions::Dominated(SM &approxGraph,
                                                    VertexSubset &remaining_vertices,
                                                    int32_t *vertexDominates,
                                                    int32_t *solution,
                                                    int32_t b_size,
                                                    std::tuple<el_t, el_t> *edgesToRemove,
                                                    int32_t &removeCounter){
  int64_t n = approxGraph.get_rows(); 
  parallel_for(int64_t i = 0; i < n; i++) { vertexDominates[i] = 0; }
  VertexSubset dominates = approxGraph.edgeMap(remaining_vertices, SET_DOMINATED_F(vertexDominates, approxGraph), true, 20);
  bool vertexChanged = false;
  VertexSubset vertices_to_delete;
  while (dominates.non_empty()) { // loop until no dominates remain
    printf("DOMINATING VERTICES\n");
    dominates.print();
    vertexChanged = true;
    removeCounter = 0;
    //__sync_fetch_and_and(&b_used, 0);
    // returns vertices to delete
    vertices_to_delete = approxGraph.edgeMap(remaining_vertices, DELETE_VERTEX_F(vertexDominates, solution, edgesToRemove, &removeCounter, &b_size, approxGraph), true, 20);
    // Write phase
    approxGraph.remove_batch(edgesToRemove, min(removeCounter, b_size));
    dominates = approxGraph.edgeMap(remaining_vertices, SET_DOMINATED_F(vertexDominates, approxGraph), true, 20);
  }
  return vertexChanged;
}
*/
//template <typename SM> int32_t VC_Reductions::RemoveMaxApproximateMVC(SM &G){
/*
template <typename SM> bool VC_Reductions::Struction(SM &approxGraph,
                                                    VertexSubset &remaining_vertices,
                                                    int32_t *numberAntiEdges,
                                                    int32_t *performStruction,
                                                    int32_t *maxVertex,
                                                    int32_t *numStructionNeighbors,
                                                    int32_t b_size,
                                                    std::tuple<el_t, el_t> *edgesToRemove,
                                                    std::tuple<el_t, el_t> *edgesToInsert,
                                                    int32_t &removeCounter,
                                                    int32_t &insertCounter){
  bool structionPerformed = false;
  int64_t n = approxGraph.get_rows(); 

  parallel_for(int64_t i = 0; i < n; i++) { numberAntiEdges[i] = 0; }
  parallel_for(int64_t i = 0; i < n; i++) { performStruction[i] = 0; }
  parallel_for(int64_t i = 0; i < n; i++) { numStructionNeighbors[i] = 0; }

  parallel_for(int64_t i = 0; i < n; i++) { maxVertex[i] = 0; }

  VertexSubset struction = approxGraph.edgeMap(remaining_vertices, SET_ANTI_EDGES_F(numberAntiEdges, approxGraph), true, 20);
  parallel_for(int64_t i = 0; i < n; i++) { numberAntiEdges[i] /= 2; }
  VertexSubset firstStructionSet = approxGraph.vertexMap(remaining_vertices, SET_STRUCTION_F(numberAntiEdges, performStruction, approxGraph), true); // mark visited
  #ifdef NDEBUG
  printf("Vertices\n");
  for (uint32_t j = 0; j < n; j++) {
    printf("%lu ", j);
  }
  printf("\nStruction flags\n");
  for (uint32_t j = 0; j < n; j++) {
    printf("%lu ", performStruction[j]);
  }
  printf("\n");
  #endif
  VertexSubset structionMIS = approxGraph.edgeMap(remaining_vertices, SOLVE_MIS_F(performStruction, approxGraph), true, 20);
  #ifdef NDEBUG
  printf("MIS\n");
  printf("Vertices\n");
  for (uint32_t j = 0; j < n; j++) {
    printf("%lu ", j);
  }
  printf("\nStruction flags\n");
  for (uint32_t j = 0; j < n; j++) {
    printf("%lu ", performStruction[j]);
  }
  printf("\n");
  #endif
  VertexSubset structionDeg = approxGraph.edgeMap(remaining_vertices, SET_NUM_STRUCTION_NEIGHBORS_F(performStruction, numStructionNeighbors, approxGraph), true, 20);
  #ifdef NDEBUG
  printf("Degree of struct\n");
  printf("Vertices\n");
  for (uint32_t j = 0; j < n; j++) {
    printf("%lu ", j);
  }
  printf("\nNum neighbors with set flags\n");
  for (uint32_t j = 0; j < n; j++) {
    printf("%lu ", numStructionNeighbors[j]);
  }
  printf("\n");
  #endif
  VertexSubset maxDegree = approxGraph.edgeMap(remaining_vertices, SET_LARGEST_VERTEX_STRUCT_F(maxVertex, numStructionNeighbors, performStruction, approxGraph), true, 20);
  #ifdef NDEBUG
  printf("\nmaxDegree\n");
  for (uint32_t j = 0; j < n; j++) {
    printf("%lu ", maxVertex[j]);
  }
  printf("\n");
  #endif
  VertexSubset fin = approxGraph.edgeMap(remaining_vertices, RESOLVE_CONFLICTS_STRUCT_F(maxVertex, numStructionNeighbors, performStruction, approxGraph), true, 20);
  #ifdef NDEBUG
  printf("\nResolve conflicts\n");
  for (uint32_t j = 0; j < n; j++) {
    printf("%lu ", performStruction[j]);
  }
  printf("\n");  
  #endif
  VertexSubset structionSetAndNeighbors = approxGraph.vertexMap(remaining_vertices, GET_STRUCTION_SET_AND_NEIGHBORS_F(performStruction, numStructionNeighbors, approxGraph), true); // mark visited
  //VertexSubset structionSetAndNeighborsDeleted = approxGraph.edgeMap(remaining_vertices, DELETE_NEIGHBORHOOD_OF_STRUCTION_VERTEX_F(performStruction, maxVertex, numStructionNeighbors, edgesToRemove, &b_used, &b_size, approxGraph), true); // mark visited
  VertexSubset structionSetAndNeighborsDeleted = approxGraph.edgeMap(structionSetAndNeighbors, DELETE_ALL_VERTICES_IN_VERTEX_SUBSET_F(edgesToRemove, &removeCounter, &b_size, approxGraph), true); // mark visited

  VertexSubset structionSetOnly = approxGraph.vertexMap(remaining_vertices, GET_STRUCTION_SET_F(performStruction, numStructionNeighbors, approxGraph), true); // mark visited

  // Assuming all the insertions fit in one batch, this should complete the struction op.
  while (structionSetOnly.non_empty()){
    structionPerformed = true;
    uint32_t v0 = structionSetOnly.pop();
    //printf("Perform struct operation on %lu\n", v0);
    //approxGraph.print_neighbors(v0);
    std::vector<el_t> v0_neighs = approxGraph.get_neighbors(v0);
    #ifdef NDEBUG
    printf("Neighbors of %lu\n", v0);
    for (int i = 0; i < v0_neighs.size(); ++i)
      printf("%u \n", v0_neighs[i]);
    printf("\n");
    #endif
    int usedVertexCounter = 0;
    std::map<std::tuple<el_t,el_t>,el_t> antiEdgeToNodeMap;
    for (int i = 0; i < v0_neighs.size(); ++i)
      for (int j = i+1; j < v0_neighs.size(); ++j){
        if (!approxGraph.has(v0_neighs[i], v0_neighs[j]))
          antiEdgeToNodeMap[std::tuple<el_t, el_t>{v0_neighs[i], v0_neighs[j]}] = v0_neighs[usedVertexCounter++];
      }
    for (auto& t : antiEdgeToNodeMap)
        std::cout << std::get<0>(t.first) << " " 
                  << std::get<1>(t.first) << " " 
                  << t.second << "\n";


    std::map<std::tuple<el_t,el_t>,el_t>::iterator it_i = antiEdgeToNodeMap.begin();
    std::map<std::tuple<el_t,el_t>,el_t>::iterator it_j;

    insertCounter = 0;

    while (it_i != antiEdgeToNodeMap.end())
    {
      // Condition 1 - 
      // remove the vertices {v0; v1; ... ; vp} from G and
      // introduce a new node vij for every anti-edge {vi; vj} in G
      // where 0 < i < j <= p;

      it_j = std::next(it_i, 1);
      while (it_j != antiEdgeToNodeMap.end())
      {
        std::cout << std::get<0>(it_i->first) << " " 
          << std::get<1>(it_i->first) << " " 
          << it_i->second << "\n"; 
        std::cout << std::get<0>(it_j->first) << " " 
          << std::get<1>(it_j->first) << " " 
          << it_j->second << "\n"; 

        // Condition 2 - 
        // add an edge (vir; vjs) if i = j and
        // (vr; vs) is an edge in G;
        if ((std::get<0>(it_i->first) == std::get<0>(it_j->first))
              && approxGraph.has(std::get<1>(it_i->first), std::get<1>(it_j->first))){
          edgesToInsert[insertCounter++] = std::tuple<el_t, el_t>{it_i->second, it_j->second};
          edgesToInsert[insertCounter++] = std::tuple<el_t, el_t>{it_j->second, it_i->second};
          printf("Adding edge (%lu, %lu)-(%lu, %lu) in new neighborhood\n",std::get<0>(it_i->first), std::get<1>(it_i->first),
          std::get<0>(it_j->first), std::get<1>(it_j->first));
        }

        // Condition 3 - 
        // if i != j, add an edge (vir; vjs)
        if (std::get<0>(it_i->first) != std::get<0>(it_j->first)){
          edgesToInsert[insertCounter++] = std::tuple<el_t, el_t>{it_i->second, it_j->second};
          edgesToInsert[insertCounter++] = std::tuple<el_t, el_t>{it_j->second, it_i->second};
          printf("Adding edge (%lu, %lu)-(%lu, %lu) in new neighborhood\n",std::get<0>(it_i->first), std::get<1>(it_i->first),
          std::get<0>(it_j->first), std::get<1>(it_j->first));
        }
        ++it_j;
      }
      // Condition 4 - 
      // for every u not in {v0; ... ; vp}, 
      // add the edge (vij ; u) if (vi; u)
      // or (vj ; u) is an edge in approxGraph.

      // exnoi
      std::cout << "Find external neighbors of i " << std::get<0>(it_i->first) << std::endl;

      std::vector<el_t> externalNeighborsOf_i = approxGraph.get_neighbors(std::get<0>(it_i->first));

      std::cout << "External neighbors of " << std::get<0>(it_i->first) << " before removal" << std::endl;
      for (auto element : externalNeighborsOf_i) {
        std::cout << element << " ";
      }
      std::cout << std::endl;

      externalNeighborsOf_i.erase( remove_if( begin(externalNeighborsOf_i),end(externalNeighborsOf_i),
          [&](auto x){return find(begin(v0_neighs),end(v0_neighs),x)!=end(v0_neighs);}), end(externalNeighborsOf_i) );

      std::cout << "External neighbors of " << std::get<0>(it_i->first) << " after N(v0) removal" << std::endl;
      for (auto element : externalNeighborsOf_i) {
        std::cout << element << " ";
      }
      std::cout << std::endl;

      externalNeighborsOf_i.erase(std::remove(externalNeighborsOf_i.begin(), externalNeighborsOf_i.end(), v0), externalNeighborsOf_i.end());

      std::cout << "External neighbors of " << std::get<0>(it_i->first) << " after v0 removal" << std::endl;
      for (auto element : externalNeighborsOf_i) {
        std::cout << element << " ";
      }
      std::cout << std::endl;

      // exnoj
      std::vector<el_t> externalNeighborsOf_j = approxGraph.get_neighbors(std::get<1>(it_i->first));
      std::cout << "Find external neighbors of j " << std::get<1>(it_i->first) << std::endl;
      std::cout << "External neighbors of " << std::get<1>(it_i->first) << " before removal" << std::endl;
      for (auto element : externalNeighborsOf_j) {
        std::cout << element << " ";
      }
      std::cout << std::endl;

      externalNeighborsOf_j.erase( remove_if( begin(externalNeighborsOf_j),end(externalNeighborsOf_j),
          [&](auto x){return find(begin(v0_neighs),end(v0_neighs),x)!=end(v0_neighs);}), end(externalNeighborsOf_j) );

      std::cout << "External neighbors of " << std::get<1>(it_i->first) << " after N(v0) removal" << std::endl;
      for (auto element : externalNeighborsOf_j) {
        std::cout << element << " ";
      }
      std::cout << std::endl;

      externalNeighborsOf_j.erase(std::remove(externalNeighborsOf_j.begin(), externalNeighborsOf_j.end(), v0), externalNeighborsOf_j.end());

      std::cout << "External neighbors of " << std::get<1>(it_i->first) << " after v0 removal" << std::endl;
      for (auto element : externalNeighborsOf_j) {
        std::cout << element << " ";
      }
      std::cout << std::endl;

      // Using default function
      std::vector<el_t> unionOfExternalNeighborsOfIAndJ;
      std::set_union(externalNeighborsOf_i.begin(), externalNeighborsOf_i.end(), externalNeighborsOf_j.begin(), externalNeighborsOf_j.end(), std::back_inserter(unionOfExternalNeighborsOfIAndJ));
      std::cout << "External neighbors of I and J are :\n";
      for (const auto &i : unionOfExternalNeighborsOfIAndJ) {
          std::cout << i << ' ';
          edgesToInsert[insertCounter++] = std::tuple<el_t, el_t>{it_i->second, i};
          edgesToInsert[insertCounter++] = std::tuple<el_t, el_t>{i, it_i->second};
      }   
      std::cout << '\n';

      ++it_i;
    }
  }

  //printf("\nBefore batch changes\n");
  //approxGraph.print_arrays();
  approxGraph.remove_batch(edgesToRemove, min(removeCounter, b_size));
  approxGraph.insert_batch(edgesToInsert, min(insertCounter, b_size));
  //printf("\nAfter batch changes\n");
  //approxGraph.print_arrays();

  return structionPerformed;
}
*/



template <typename SM> bool VC_Reductions::GeneralFold(SM &approxGraph,
                                                    VertexSubset &remaining_vertices,
                                                    int32_t *request,
                                                    int32_t *match,
                                                    int32_t *auxMatch,
                                                    int32_t *H_n,
                                                    int32_t b_size,
                                                    std::tuple<el_t, el_t> *edgesToRemove,
                                                    std::tuple<el_t, el_t> *edgesToInsert,
                                                    int32_t &removeCounter,
                                                    int32_t &insertCounter){
    
  int64_t n = approxGraph.get_rows(); 
  int64_t count = 0; 

  parallel_for(int64_t i = 0; i < n; i++) { auxMatch[i] = 0; }
  parallel_for(int64_t i = 0; i < n; i++) { H_n[i] = 0; }

  Match(approxGraph,
        remaining_vertices,
        request,
        match,
        H_n);
  printf("maximal matching M1\n");
  for(int64_t i = 0; i < n; i++) { if(match[i] >= 4) printf("%lu ", i); }
  printf("\n");
  printf("the set O of outsiders\n");
  for(int64_t i = 0; i < n; i++) { if(match[i] < 4) printf("%lu ", i); }
  printf("\n");
  bool step3Crown = AuxilliaryMatch(approxGraph,
                  remaining_vertices,
                  request,
                  match,
                  auxMatch,
                  H_n);

  printf("maximal matching M2\n");
  for(int64_t i = 0; i < n; i++) { if(auxMatch[i] >= 4 && auxMatch[i] != n) printf("%lu ", i); }
  printf("\n");
  printf("the set of vertices unmatched by M2\n");
  for(int64_t i = 0; i < n; i++) { if(auxMatch[i] < 4) printf("%lu ", i); }
  printf("\n");
  bool c = FindCrown(approxGraph,
                    remaining_vertices,
                    request,
                    match,
                    auxMatch,
                    H_n);
  printf("Exited do while loop\n");
  printf("H\n");
  for(int64_t i = 0; i < n; i++) { if(H_n[i]) printf("%lu ", i); }
  printf("\n");
  printf("C\n");
  for(int64_t i = 0; i < n; i++) { if(auxMatch[i] == 1) printf("%lu ", i); }
  printf("\n");
}


template <typename SM> bool VC_Reductions::Match(SM &approxGraph,
                                                    VertexSubset &remaining_vertices,
                                                    int32_t *request,
                                                    int32_t *match,
                                                    int32_t *maxDegree){
  int count = 0;
  int UL = 256;

  int64_t n = approxGraph.get_rows(); 

  parallel_for(int64_t i = 0; i < n; i++) { match[i] = 0; }

  //parallel_for(int64_t i = 0; i < n; i++) { numStructionNeighbors[i] = 0; }

  VertexSubset unmatchedVertices;
  do {
    uint randomNumber = rand();
    parallel_for(int64_t i = 0; i < n; i++) { request[i] = n; }
    //parallel_for(int64_t i = 0; i < n; i++) { maxDegree[i] = INT32_MAX; }

    //printf("match round %d\n", count);
    unmatchedVertices = approxGraph.vertexMap(remaining_vertices, SELECT_COLOR_F(match, approxGraph, randomNumber), true); // mark visited
    /*
    printf("Unmatched verts\n");
    unmatchedVertices.print();
    printf("vert colors\n");
    for(int64_t i = 0; i < n; i++) { printf("%lu %u\n", i, match[i]); }
    printf("\n");

    */

    //approxGraph.edgeMap(remaining_vertices, GET_MIN_DEGREE_F(match, maxDegree, approxGraph), false, 20);
    //approxGraph.edgeMap(remaining_vertices, REQUEST_2_F(match, request, maxDegree, approxGraph), false, 20);
    approxGraph.edgeMap(remaining_vertices, REQUEST_3_F(match, request, maxDegree, approxGraph), false, 20);

    /*
    printf("vert REQUEST after request\n");
    for(int64_t i = 0; i < n; i++) { printf("%lu %u\n", i, request[i]); }
    printf("\n");
    */
    //VertexSubset struction2 = approxGraph.edgeMap(remaining_vertices, RESPOND_2_F(match, request, approxGraph), true, 20);
    //approxGraph.edgeMap(remaining_vertices, RESPOND_2_F(match, request, approxGraph), false, 20);
    approxGraph.edgeMap(remaining_vertices, RESPOND_2_F(match, request, approxGraph), false, 20);

    /*
    printf("vert REQUEST after respond\n");
    for(int64_t i = 0; i < n; i++) { printf("%lu %u\n", i, request[i]); }
    printf("\n");
    */
    unmatchedVertices = approxGraph.vertexMap(remaining_vertices, MATCH_F(match, request, approxGraph), true); // mark visited
  } while (unmatchedVertices.non_empty() && ++count < UL);
  /*
  printf("Unmatched verts\n");
  unmatchedVertices.print();  
  printf("vert requests\n");
  for(int64_t i = 0; i < n; i++) { printf("%lu %u\n", i, request[i]); }
  printf("\n");
  printf("vert colors\n");
  for(int64_t i = 0; i < n; i++) { printf("%lu %u\n", i, match[i]); }
  printf("\n");
  */
  return true;
}


template <typename SM> bool VC_Reductions::AuxilliaryMatch(SM &approxGraph,
                                                    VertexSubset &remaining_vertices,
                                                    int32_t *request,
                                                    int32_t *match,
                                                    int32_t *auxMatch,
                                                    int32_t *maxDegree){
  int count = 0;
  int UL = 256;

  int64_t n = approxGraph.get_rows(); 

  parallel_for(int64_t i = 0; i < n; i++) { auxMatch[i] = n; }

  VertexSubset unmatchedVertices;

  //printf("match round %d\n", count);
  approxGraph.edgeMap(remaining_vertices, SELECT_COLOR_AUX_N_O_F(match, auxMatch, approxGraph), false, 20); // mark visited

  // UMV only consists of O and N(O)
  unmatchedVertices = approxGraph.vertexMap(remaining_vertices, SELECT_COLOR_AUX_F(match, auxMatch, approxGraph), true); // mark visited
  // Just for printing later.
  VertexSubset N_O = approxGraph.vertexMap(unmatchedVertices, GET_UNMATCHED_N_O(auxMatch, approxGraph), true); // mark visited

  do {
    uint randomNumber = rand();

    parallel_for(int64_t i = 0; i < n; i++) { request[i] = n; }
    //parallel_for(int64_t i = 0; i < n; i++) { maxDegree[i] = INT32_MAX; }

    /*
    printf("Unmatched verts\n");
    unmatchedVertices.print();
    printf("vert colors\n");
    for(int64_t i = 0; i < n; i++) { printf("%lu %u\n", i, match[i]); }
    printf("\n");
    */
    //approxGraph.edgeMap(remaining_vertices, GET_MIN_DEGREE_F(auxMatch, maxDegree, approxGraph), false, 20);
    //approxGraph.edgeMap(remaining_vertices, REQUEST_2_F(auxMatch, request, maxDegree, approxGraph), false, 20);
    approxGraph.edgeMap(remaining_vertices, REQUEST_3_F(match, request, maxDegree, approxGraph), false, 20);

    /*
    printf("vert REQUEST after request\n");
    for(int64_t i = 0; i < n; i++) { printf("%lu %u\n", i, request[i]); }
    printf("\n");
    */
    //VertexSubset struction2 = approxGraph.edgeMap(remaining_vertices, RESPOND_2_F(match, request, approxGraph), true, 20);
    approxGraph.edgeMap(remaining_vertices, RESPOND_2_F(auxMatch, request, approxGraph), false, 20);

    /*
    printf("vert REQUEST after respond\n");
    for(int64_t i = 0; i < n; i++) { printf("%lu %u\n", i, request[i]); }
    printf("\n");
    */
    unmatchedVertices = approxGraph.vertexMap(remaining_vertices, MATCH_F(auxMatch, request, approxGraph), true); // mark visited
  } while (unmatchedVertices.non_empty() && ++count < UL);

  VertexSubset unmatched_N_O = approxGraph.vertexMap(unmatchedVertices, GET_UNMATCHED_N_O(auxMatch, approxGraph), true); // mark visited

  printf("Unmatched verts\n");
  unmatchedVertices.print();  
  printf("unmatched_N_O verts\n");
  unmatched_N_O.print(); 
  if (!unmatched_N_O.non_empty()){
    printf("Identified a crown!\n");
    printf("Neighbors of O\n");
    N_O.print();
  }
  /*
  printf("vert requests\n");
  for(int64_t i = 0; i < n; i++) { printf("%lu %u\n", i, request[i]); }
  printf("\n");
  printf("vert colors\n");
  for(int64_t i = 0; i < n; i++) { printf("%lu %u\n", i, auxMatch[i]); }
  printf("\n");
  */
  return !unmatched_N_O.non_empty();
}


template <typename SM> bool VC_Reductions::FindCrown(SM &approxGraph,
                                                    VertexSubset &remaining_vertices,
                                                    int32_t *request,
                                                    int32_t *match,
                                                    int32_t *auxMatch,
                                                    int32_t *H_n){
  int count = 0;
  int UL = 256;

  int64_t n = approxGraph.get_rows(); 


  parallel_for(int64_t i = 0; i < n; i++) { H_n[i] = 0; }
  parallel_for(int64_t i = 0; i < n; i++) { request[i] = 0; }
  //parallel_for(int64_t i = 0; i < n; i++) { auxMatch[i] = auxMatch[i] == 1; }
  //VertexSubset struction = approxGraph.edgeMap(remaining_vertices, REQUEST_2_F(match, request, approxGraph), true, 20);
  printf("I_0\n");
  for(int64_t i = 0; i < n; i++) { if(auxMatch[i] == 1) printf("%lu ", i); }
  printf("\n");
  bool I_N_expanded = false;
  do {
    // Step 5a
    approxGraph.edgeMap(remaining_vertices, SET_NEIGHBORS_OF_UNMATCHED_O_F(auxMatch, H_n, approxGraph), false, 20);
    printf("Hn\n");
    for(int64_t i = 0; i < n; i++) { if(H_n[i]) printf("%lu ", i); }
    printf("\n");
    // Set NM2(Hn).
    approxGraph.edgeMap(remaining_vertices, SET_NEIGHBORS_WITHIN_AUX_MATCHING_F(auxMatch, H_n, request, approxGraph), false, 20);
    printf("NM2(Hn)\n");
    for(int64_t i = 0; i < n; i++) { if(request[i]) printf("%lu ", i); }
    printf("\n");
    // Step 5b
    // Repeat steps 5a and 5b until n = N so that In-1 = In.
    // In-1 = In; means no expansion happened.
    // N iterations should be the upper limit.
    I_N_expanded = false;
    for(int64_t i = 0; i < n; i++) { 
      if (auxMatch[i] != 1 && request[i] == 1){
        I_N_expanded = true;
        auxMatch[i] = 1;
      }
    }
  } while (I_N_expanded && ++count < n);


  return true;
}


//template <typename SM> int32_t VC_Reductions::RemoveMaxApproximateMVC(SM &G){
template <typename SM> int32_t* VC_Reductions::RemoveMaxApproximateMVC(SM &G){

  SparseMatrixV<true, bool> approxGraph(G);
  int64_t n = approxGraph.get_rows(); 
  if (n == 0) {
    return 0;
  }

  // creates inCover array, initialized to all -1, except for start
  int32_t b_size = 10000; 
  int32_t b_used = 0; 
  std::tuple<el_t, el_t> *edgesToRemove = (std::tuple<el_t, el_t> *)malloc(b_size * sizeof(std::tuple<el_t, el_t>));

  int32_t *isLeaf = (int32_t *)malloc(n * sizeof(int32_t));
  int32_t *inCover = (int32_t *)malloc(n * sizeof(int32_t));
  int32_t *solution = (int32_t *)malloc(n * sizeof(int32_t));

  parallel_for(int64_t i = 0; i < n; i++) { isLeaf[i] = 0; }
  parallel_for(int64_t i = 0; i < n; i++) { inCover[i] = 0; }
  parallel_for(int64_t i = 0; i < n; i++) { solution[i] = 0; }

  VertexSubset remaining_vertices = VertexSubset(0, n, true); // initial set contains all vertices
  VertexSubset vertices_to_delete;

  //SM approxGraph(G);
  ProgressBar progressBar;
  progressBar.setAlgorithmStartSize(n);
  //progressBar.printNumberOfRemainingVertices(remaining_vertices.get_n());
  //std::cout << std::endl;
	int32_t minimum = 0;
	bool hasEdges = true;
	while (hasEdges)
	{
		bool leafHasChanged = false, triangleHasChanged = false;
		unsigned int iterationCounter = 0;
		do {
      leafHasChanged = false;
      triangleHasChanged = false;

      //printf("Degrees\n");

      approxGraph.vertexMap(remaining_vertices, SET_LEAVES_1_F(isLeaf, approxGraph), false); // mark visited
      VertexSubset leaves = approxGraph.edgeMap(remaining_vertices, BF_LEAVES_F(isLeaf, inCover), true, 20);
      /*
      if (leaves.non_empty()){
        for (unsigned int i = 0; i < n; i++)
          printf("%d ", approxGraph.getDegree(i));
        printf("\n");
      }
      */
      
      while (leaves.non_empty()) { // loop until no leaves remain
        leafHasChanged = true;
        b_used = 0;
        //__sync_fetch_and_and(&b_used, 0);
        // returns vertices to delete
        vertices_to_delete = approxGraph.edgeMap(remaining_vertices, DELETE_VERTEX_F(inCover, solution, edgesToRemove, &b_used, &b_size, approxGraph), true, 20);
        // Write phase
        approxGraph.remove_batch(edgesToRemove, min(b_used, b_size));
        //if(percentage < 95.0)
        //  progressBar.printIterationBar(vertices_to_delete.get_n());
        approxGraph.vertexMap(remaining_vertices, SET_LEAVES_1_F(isLeaf, approxGraph), false); // mark visited
        leaves = approxGraph.edgeMap(remaining_vertices, BF_LEAVES_F(isLeaf, inCover), true, 20);
      }
      /*
      printf("IsLeaf\n");
      for (unsigned int i = 0; i < n; i++)
        printf("%d ", isLeaf[i]);
      printf("\n");
      printf("In cover\n");
      for (unsigned int i = 0; i < n; i++)
        printf("%d ", inCover[i]);
      printf("\n");
      printf("NumLeaves %u\n",leaves.get_n());
			*/
      // sets all triangles
      parallel_for(int64_t i = 0; i < n; i++) { inCover[i] = 0; }
      VertexSubset triangles = approxGraph.edgeMap(remaining_vertices, SET_TRIANGLES_F(isLeaf, solution, edgesToRemove, &b_used, &b_size, approxGraph), true, 20);
      // Prune any extra triangles - i.e. for a single connected component triangle,  only 2 vertices should be included
      leaves = approxGraph.edgeMap(remaining_vertices, BF_TRIANGLES_F(isLeaf, inCover), true, 20);

      while (triangles.non_empty()) { // loop until no leaves remain
        triangleHasChanged = true;
        b_used = 0;
        // returns vertices to delete
        vertices_to_delete = approxGraph.edgeMap(remaining_vertices, DELETE_VERTEX_F(inCover, solution, edgesToRemove, &b_used, &b_size, approxGraph), true, 20);
        // Write phase
        approxGraph.remove_batch(edgesToRemove, min(b_used, b_size));
        parallel_for(int64_t i = 0; i < n; i++) { inCover[i] = 0; }
        // sets all triangles
        triangles = approxGraph.edgeMap(remaining_vertices, SET_TRIANGLES_F(inCover, solution, edgesToRemove, &b_used, &b_size, approxGraph), true, 20);
        // Prune any extra triangles - i.e. for a single connected component triangle,  only 2 vertices should be included
        leaves = approxGraph.edgeMap(remaining_vertices, BF_TRIANGLES_F(isLeaf, inCover), true, 20);
      }
      
      //parallel_for(int64_t i = 0; i < n; i++) { inCover[i] = 0; }
      /*
      VertexSubset triangles = approxGraph.edgeMap(remaining_vertices, SET_TRIANGLES_F(inCover, solution, edgesToRemove, &b_used, &b_size, approxGraph), true, 20);
			while (triangles.non_empty()) { // loop until no leaves remain
        b_used = 0; 
        for (unsigned int i = 0; i < n; i++)
          printf("%d ", approxGraph.getDegree(i));
        printf("\n");
        vertices_to_delete = approxGraph.edgeMap(triangles, TRIANGLE_REDUCTION_RULE_F(inCover, solution, edgesToRemove, &b_used, &b_size, approxGraph), true, 20);
        if (vertices_to_delete.non_empty()){
          triangleHasChanged = true;
        }
        // Write phase
        approxGraph.remove_batch(edgesToRemove, min(b_used, b_size));
        triangles = approxGraph.edgeMap(remaining_vertices, SET_TRIANGLES_F(inCover, solution, edgesToRemove, &b_used, &b_size, approxGraph), true, 20);
      }
      */
		} while (leafHasChanged || triangleHasChanged);
		//} while (leafHasChanged);

		int32_t maxV;
		int32_t maxD = 0;
		for (unsigned int i = 0; i < n; i++)
		{
			if (approxGraph.getDegree(i) > maxD)
			{
				maxV = i;
				maxD = approxGraph.getDegree(i);
			}
		}
		if (maxD == 0)
			hasEdges = false;
		else
		{
      /*
      for (unsigned int i = 0; i < n; i++)
        printf("%d ", approxGraph.getDegree(i));
      printf("\n");
      */
      VertexSubset maxVSet = approxGraph.vertexMap(remaining_vertices, SET_MAX_VERTEX_F(inCover, maxV, approxGraph), true); // mark visited
			//approxGraph.deleteVertex(maxV);
      bool firstBatchIteration = true; 
      float percentage = progressBar.getAlgorithmPercentage(remaining_vertices.get_n());
      do {
        b_used = 0;
        //__sync_fetch_and_and(&b_used, 0);
        // returns vertices to delete
        vertices_to_delete = approxGraph.edgeMap(remaining_vertices, DELETE_VERTEX_F(inCover, solution, edgesToRemove, &b_used, &b_size, approxGraph), true, 20);
        if (firstBatchIteration && percentage < 95.0){
          progressBar.setIterationStartSize(vertices_to_delete.get_n());
          firstBatchIteration = false;
        }
        // Write phase
        approxGraph.remove_batch(edgesToRemove, min(b_used, b_size));
        //if(percentage < 95.0)
        //  progressBar.printIterationBar(vertices_to_delete.get_n());
      } while(vertices_to_delete.non_empty());
			++minimum;
		}
	}
  // Destructor is automatically called
	//approxGraph.del();

  free(edgesToRemove);
  free(inCover);
	return solution;
}
