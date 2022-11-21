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
    // I'm a leaf
    if (inCover[d]){
      //printf("LEAF %u %"PRIu32"\n", d, G.getDegree(d));
      // and you're a leaf
      if (inCover[s]){
        //printf("DOUBLE LEAF %u %u (%"PRIu32") (%"PRIu32")\n", s, d, G.getDegree(s), G.getDegree(d));
        // tiebreak
        if (h(s) < h(d)){
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
      } else {
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
    }
    return deletedEdge;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    printf("USED ATOMIC\n");
    exit(1);
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
    isTriangle[d] |= G.getDegree(d) == 2 && G.common_neighbors(s,d);
    return isTriangle[d];
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    // Assuming I am not a neighbor to myself,
    // thus G.common_neighbors(s,d) > 0 indicates a triangle.
    isTriangle[s] |= G.getDegree(s) == 2 && G.common_neighbors(s,d);
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
    if (isTriangle[s]){
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



class VC_Reductions {
  public:
    template <typename SM> int32_t RemoveMaxApproximateMVC(SM &G);
};

//template <typename SM> int32_t VC_Reductions::RemoveMaxApproximateMVC(SM &G){
template <typename SM> int32_t VC_Reductions::RemoveMaxApproximateMVC(SM &approxGraph){
  //int64_t n = G.get_rows();

  int64_t n = approxGraph.get_rows();
  if (n == 0) {
    return 0;
  }

  // creates inCover array, initialized to all -1, except for start
  int32_t b_size = 10000; 
  int32_t b_used = 0; 
  std::tuple<el_t, el_t> *edgesToRemove = (std::tuple<el_t, el_t> *)malloc(b_size * sizeof(std::tuple<el_t, el_t>));

  int32_t *inCover = (int32_t *)malloc(n * sizeof(int32_t));
  int32_t *solution = (int32_t *)malloc(n * sizeof(int32_t));

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
      VertexSubset leaves = approxGraph.vertexMap(remaining_vertices, SET_LEAVES_F(inCover, approxGraph), true); // mark visited
      //printf("NumLeaves %u\n",leaves.get_n());
			while (leaves.non_empty()) { // loop until no leaves remain
        b_used = 0; 
        vertices_to_delete = approxGraph.edgeMap(remaining_vertices, LEAF_REDUCTION_RULE_F(inCover, solution, edgesToRemove, &b_used, &b_size, approxGraph), true, 20);
        if (vertices_to_delete.non_empty())
          leafHasChanged = true;
        // Write phase
        approxGraph.remove_batch(edgesToRemove, min(b_used, b_size));
        leaves = approxGraph.vertexMap(remaining_vertices, SET_LEAVES_F(inCover, approxGraph), true); // mark visited
      }
    /*
      VertexSubset triangles = approxGraph.edgeMap(remaining_vertices, SET_TRIANGLES_F(inCover, solution, edgesToRemove, &b_used, &b_size, approxGraph), true, 20);
      //printf("NumPossTri %u\n",possibleTriangles.get_n());
			while (triangles.non_empty()) { // loop until no leaves remain
        b_used = 0; 
        vertices_to_delete = approxGraph.edgeMap(triangles, TRIANGLE_REDUCTION_RULE_F(inCover, solution, edgesToRemove, &b_used, &b_size, approxGraph), true, 20);
        if (vertices_to_delete.non_empty())
          triangleHasChanged = true;
        // Write phase
        approxGraph.remove_batch(edgesToRemove, min(b_used, b_size));
        triangles = approxGraph.edgeMap(remaining_vertices, SET_TRIANGLES_F(inCover, solution, edgesToRemove, &b_used, &b_size, approxGraph), true, 20);
      }
		} while (leafHasChanged || triangleHasChanged);
    */
		} while (leafHasChanged);

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
        if(percentage < 95.0)
          progressBar.printIterationBar(vertices_to_delete.get_n());
      } while(vertices_to_delete.non_empty());
			++minimum;
		}

	}
  // Destructor is automatically called
	//approxGraph.del();
  free(edgesToRemove);
  free(inCover);
  free(solution);
	return minimum;
}
