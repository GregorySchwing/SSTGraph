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
struct H_F {
  int32_t *Parents;
  explicit H_F(int32_t *_Parents) : Parents(_Parents) {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    if (Parents[d] == -1) {
      Parents[d] = s;
      return true;
    }
    return false;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    return __sync_bool_compare_and_swap(&Parents[d], -1, s);
  }
  // cond function checks if vertex has been visited yet
  inline bool cond(uint32_t d) { return (Parents[d] == -1); }
};

struct I_F {
  int32_t *Parents;
  int32_t *match;
  explicit I_F(int32_t *_Parents, int32_t *_match) : Parents(_Parents), match(_match) {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    if (Parents[d] == -1) {
      Parents[d] = s;
      return true;
    }
    return false;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    return __sync_bool_compare_and_swap(&Parents[d], -1, s);
  }
  // cond function checks if vertex has been visited yet
  // also destination should be in M.
  inline bool cond(uint32_t d) { return (Parents[d] == -1 && match[d] != -1); }
};

struct CYCLE_1_F {
  int32_t *Parents;
  int32_t *Pairs;
  explicit CYCLE_1_F(int32_t *_Parents, int32_t *_Pairs) : Parents(_Parents), Pairs(_Pairs) {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    if (Parents[d] == -1) {
      Parents[d] = s;
      return true;
    }
    return false;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    return __sync_bool_compare_and_swap(&Parents[d], -1, s);
  }
  // cond function checks if vertex has been visited yet
  inline bool cond(uint32_t d) { return (Parents[d] == -1); }
};

template <typename SM> int32_t *CR_with_edge_map(const SM &G, int* match, uint32_t src) {
  int64_t start = src;
  int64_t n = G.get_rows();
  int32_t i = 0;
  if (n == 0) {
    return nullptr;
  }
  // creates Parents array, initialized to all -1, except for start
  int32_t *Parents = (int32_t *)malloc(n * sizeof(int32_t));
  // creates Parents array, initialized to all -1, except for start
  int32_t *Parents = (int32_t *)malloc(n * sizeof(int32_t));
  // creates Pairs array, initialized to all -1
  // when an edge is shared between two vertices in H or I,
  // they set each other, so they may backtrack until they converge.
  int32_t *Pairs = (int32_t *)malloc(n * sizeof(int32_t));
  parallel_for(int64_t i = 0; i < n; i++) { Parents[i] = -1; }
  parallel_for(int64_t i = 0; i < n; i++) { Pairs[i] = -1; }

  if (n == 0) {
    return Parents;
  }
  Parents[start] = start;
  VertexSubset frontier = VertexSubset(start, n); // creates initial frontier
  while (frontier.non_empty()) { // loop until frontier is empty
    VertexSubset H = G.edgeMap(frontier, H_F(Parents), true, 20);

    printf("H\n");
    H.print();
    // Check for cycles in H
    VertexSubset I = G.edgeMap(H, I_F(Parents, match), true, 20);
    printf("I\n");
    I.print();
    // Check for cycles in I
    // {M={M\{<xq, NM(xq)>}} ∪ {<NM(xq), xq−1>},
    //    ^ added by me    ^ to indicate we are removing some edges 
    // and adding another (others), and not only removing. 
    // the set difference only applies to the {<xq, NM(xq)>}} term.
    frontier.del();
    // Not sure if its this simple
    // Technically need to ensure no previous H vertices are in the frontier.
    // 5.4 else Hi=N(Ii) \ U_j=0..i−1 Hj; } 
    // possible could do a vertex map on frontier depth. though not currently tracking this.
    frontier = I;
  }
  frontier.del();
  return Parents;
}
