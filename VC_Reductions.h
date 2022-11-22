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
        inCover[d] = 1;
        return 1;
      }
    } else if (isLeaf[s] && isLeaf[d] && h(s) < h(d)){
        if (inCover[d] == 0) {
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
      exit(1);
      inCover[d] |= 1;
      return false;
    }
  }
  /*
    In sparse mode, EdgeMap
    iterates over the outgoing edges of each vertex in the subset
    and updates the destination vertex for each edge. Because it is
    run in parallel, synchronization must be used when accessing
    the destination vertex data.

    I am the source.  I update inCover[d] using synchronization.
    Degrees are const in this class.

  */
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
  exit(1);
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
          printf("DELETE LEAF %d %d (%d) (%d)\n",d, s, G.getDegree(d), G.getDegree(s));

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
        printf("DELETE LEAF %d %d (%d) (%d)\n",s,d, G.getDegree(s), G.getDegree(d));
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
    if (i == maxV)
      printf("DELETE MAX %d (%d)\n",i, G.getDegree(i));
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

      printf("Degrees\n");
      for (unsigned int i = 0; i < n; i++)
        printf("%d ", approxGraph.getDegree(i));
      printf("\n");
      approxGraph.vertexMap(remaining_vertices, SET_LEAVES_1_F(isLeaf, approxGraph), false); // mark visited
      VertexSubset leaves = approxGraph.edgeMap(remaining_vertices, BF_LEAVES_F(isLeaf, inCover), true, 20);
      printf("IsLeaf\n");
      for (unsigned int i = 0; i < n; i++)
        printf("%d ", isLeaf[i]);
      printf("\n");
      printf("In cover\n");
      for (unsigned int i = 0; i < n; i++)
        printf("%d ", inCover[i]);
      printf("\n");
      printf("NumLeaves %u\n",leaves.get_n());
			while (leaves.non_empty()) { // loop until no leaves remain
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
      
      //parallel_for(int64_t i = 0; i < n; i++) { inCover[i] = 0; }

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
		} while (leafHasChanged || triangleHasChanged);

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
      for (unsigned int i = 0; i < n; i++)
        printf("%d ", approxGraph.getDegree(i));
      printf("\n");
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
  int32_t vc_count = 0;
  for (int j = 0; j < approxGraph.get_rows(); j++) {
    vc_count += solution[j];
  }
  printf("\nVC_Red size : %u\n", vc_count);
  exit(1);
  free(edgesToRemove);
  free(inCover);
  free(solution);
	return minimum;
}
