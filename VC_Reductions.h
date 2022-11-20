#pragma once
#include "SparseMatrix.hpp"
#include "ProgressBar.h"

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
struct SET_LEAVES_F {
  T *isLeaf;
  const SM &G;

  explicit SET_LEAVES_F(T *_isLeaf, const SM &_G) : isLeaf(_isLeaf), G(_G) {}
  inline bool operator()(uintE i) {
    isLeaf[i] = G.getDegree(i) == 1;
    return isLeaf[i];
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

template <typename T, typename SM> struct LEAF_REDUCTION_RULE_F {
  T *isLeaf;
  T *numToEdgesRemove;
  T *maxNumToEdgesRemove;
  std::tuple<el_t, el_t> *edgesToRemove;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  LEAF_REDUCTION_RULE_F( T *_isLeaf, std::tuple<el_t, el_t> *_edgesToRemove, 
  T *_numToEdgesRemove, T *_maxNumToEdgesRemove, SM &_G) : 
  isLeaf(_isLeaf),
  maxNumToEdgesRemove(_maxNumToEdgesRemove),
  numToEdgesRemove(_numToEdgesRemove),
  edgesToRemove(_edgesToRemove),
  G(_G)  {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    bool deletedEdge = false;
    if (isLeaf[d]){
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
    if (isLeaf[s]){
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


template <typename T, typename SM> struct TRIANGLE_REDUCTION_RULE_F {
  T *isTriangle;
  T *numToEdgesRemove;
  T *maxNumToEdgesRemove;
  std::tuple<el_t, el_t> *edgesToRemove;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  TRIANGLE_REDUCTION_RULE_F( T *_isTriangle, std::tuple<el_t, el_t> *_edgesToRemove, 
  T *_numToEdgesRemove, T *_maxNumToEdgesRemove, SM &_G) : 
  isTriangle(_isTriangle),
  maxNumToEdgesRemove(_maxNumToEdgesRemove),
  numToEdgesRemove(_numToEdgesRemove),
  edgesToRemove(_edgesToRemove),
  G(_G)  {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    bool deletedEdge = false;
    // Assuming I am not a neighbor to myself,
    // thus G.common_neighbors(s,d) > 0 indicates a triangle.
    if (isTriangle[d] && G.common_neighbors(s,d)){
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
    if (isTriangle[s] && G.common_neighbors(s,d)){
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
    template <typename SM> int32_t RemoveMaxApproximateMVC(SM &G);
};

template <typename SM> int32_t VC_Reductions::RemoveMaxApproximateMVC(SM &G){

  int64_t n = G.get_rows();
  if (n == 0) {
    return 0;
  }

  // creates inCover array, initialized to all -1, except for start
  int32_t b_size = 10000; 
  int32_t b_used = 0; 
  std::tuple<el_t, el_t> *edgesToRemove = (std::tuple<el_t, el_t> *)malloc(b_size * sizeof(std::tuple<el_t, el_t>));

  int32_t *isLeaf = (int32_t *)malloc(n * sizeof(int32_t));

  parallel_for(int64_t i = 0; i < n; i++) { isLeaf[i] = 0; }

  VertexSubset remaining_vertices = VertexSubset(0, n, true); // initial set contains all vertices
  VertexSubset vertices_to_delete;

  SM approxGraph(G);

	int32_t minimum = 0;
	bool hasEdges = true;
	while (hasEdges)
	{
		bool leafHasChanged = false, triangleHasChanged = false;
		unsigned int iterationCounter = 0;
		do {
      leafHasChanged = false;
      triangleHasChanged = false;
      VertexSubset leaves = G.vertexMap(remaining_vertices, SET_LEAVES_F(isLeaf, G), true); // mark visited
			while (leaves.non_empty()) { // loop until no leaves remain
        b_used = 0; 
        vertices_to_delete = G.edgeMap(remaining_vertices, LEAF_REDUCTION_RULE_F(isLeaf, edgesToRemove, &b_used, &b_size, G), true, 20);
        if (vertices_to_delete.non_empty())
          leafHasChanged = true;
        // Write phase
        G.remove_batch(edgesToRemove, min(b_used, b_size));
        VertexSubset leaves = G.vertexMap(remaining_vertices, SET_LEAVES_F(isLeaf, G), true); // mark visited
      }

      VertexSubset possibleTriangles = G.vertexMap(remaining_vertices, SET_POSSIBLE_TRIANGLES_F(isLeaf, G), true); // mark visited
			while (possibleTriangles.non_empty()) { // loop until no leaves remain
        b_used = 0; 
        vertices_to_delete = G.edgeMap(possibleTriangles, TRIANGLE_REDUCTION_RULE_F(isLeaf, edgesToRemove, &b_used, &b_size, G), true, 20);
        if (vertices_to_delete.non_empty())
          triangleHasChanged = true;
        // Write phase
        G.remove_batch(edgesToRemove, min(b_used, b_size));
        VertexSubset possibleTriangles = G.vertexMap(remaining_vertices, SET_POSSIBLE_TRIANGLES_F(isLeaf, G), true); // mark visited
      }
      
		} while (leafHasChanged || triangleHasChanged);

		unsigned int maxV;
		int maxD = 0;
		for (unsigned int i = 0; i < approxGraph.get_rows(); i++)
		{
      if (maxD)
			//if (approxGraph.degree[i] > maxD)
			{
				maxV = i;
				//maxD = approxGraph.degree[i];
			}
		}
		if (maxD == 0)
			hasEdges = false;
		else
		{
			//approxGraph.deleteVertex(maxV);
			++minimum;
		}
	}

  // Destructor is automatically called
	//approxGraph.del();
  free(edgesToRemove);
  free(isLeaf);
	return minimum;
}
