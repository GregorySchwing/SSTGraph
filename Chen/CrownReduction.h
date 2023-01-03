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
template <typename SM> 
struct Update_Remaining_V_F {
  const SM &G;

  explicit Update_Remaining_V_F(const SM &_G) : G(_G) {}
  inline bool operator()(uintE i) {
    return G.getDegree(i) > 0;
  }
};

template <typename T, typename SM> 
struct GET_UNMATCHED_F {
  const SM &G;
  int *match;
  T *Cycles;
  explicit GET_UNMATCHED_F(int *_match, T *_Cycles, const SM &_G) : 
  G(_G),
  match(_match),
  Cycles(_Cycles) {}
  inline bool operator()(uintE i) {
    return match[i] == -1 && !Cycles[i];
  }
};

template <typename T, typename SM> 
struct GET_START_F {
  const SM &G;
  int *match;
  T *Cycles;
  explicit GET_START_F(int *_match, T *_Cycles, const SM &_G) : 
  G(_G),
  match(_match),
  Cycles(_Cycles) {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    if (Cycles[s] || match[s])
      return false;
    // has at least 1 neighbor to a noncycle vertex.
    return !Cycles[d];
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    if (Cycles[s] || match[s])
      return false;
    // has at least 1 neighbor to a noncycle vertex.
    return !Cycles[d];
  }
  // cond function checks if vertex has been visited yet
  inline bool cond(uint32_t d) { return true; }
};


struct H_F {
  int32_t *Parents;
  int32_t *Cycles;
  int32_t *Depth;
  explicit H_F(int32_t *_Parents, int32_t *_Cycles, int32_t *_Depth) : Parents(_Parents), Cycles(_Cycles), Depth(_Depth) {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    if (Cycles[d])
      return false;
    if (Parents[d] == -1) {
      Parents[d] = s;
      Depth[d] = Depth[s]+1;
      return true;
    }
    return false;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    if (Cycles[d])
      return false;
    return __sync_bool_compare_and_swap(&Parents[d], -1, s) && 
    __sync_bool_compare_and_swap(&Depth[d], -1, Depth[s]+1);
  }
  // cond function checks if vertex has been visited yet
  inline bool cond(uint32_t d) { return (Parents[d] == -1); }
};

struct I_F {
  int32_t *Parents;
  int32_t *Cycles;
  int32_t *Depth;
  int32_t *match;
  explicit I_F(int32_t *_Parents, int32_t *_Cycles, int32_t *_Depth, int32_t *_match) : 
  Parents(_Parents), Cycles(_Cycles),
  Depth(_Depth), match(_match) {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    if (Cycles[d])
      return false;
    if (Parents[d] == -1) {
      Parents[d] = s;
      Depth[d] = Depth[s]+1;
      return true;
    }
    return false;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    if (Cycles[d])
      return false;
    return __sync_bool_compare_and_swap(&Parents[d], -1, s) && 
    __sync_bool_compare_and_swap(&Depth[d], -1, Depth[s]+1);
  }
  // cond function checks if vertex has been visited yet
  // also destination should be in M.
  inline bool cond(uint32_t d) { return (Parents[d] == -1 && match[d] != -1); }
};

struct CYCLE_DETECTION_F {
  int32_t *Parents;
  int32_t *lock;
  int32_t *Depth;
  int32_t *CycleEdge_u;
  int32_t *CycleEdge_v;

  explicit CYCLE_DETECTION_F(int32_t *_Parents, int32_t *_lock, 
  int32_t *_Depth, int32_t *_CycleEdge_u, int32_t *_CycleEdge_v) : 
  Parents(_Parents), lock(_lock), Depth(_Depth), CycleEdge_u(_CycleEdge_u), CycleEdge_v(_CycleEdge_v) {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    //printf("%d depth %d %d depth %d\n",s, Depth[s], d, Depth[d]);
    if (Depth[d] == Depth[s] && __sync_bool_compare_and_swap(lock, -1, 0)){
      *CycleEdge_u = s;
      *CycleEdge_v = d;
      return true;
    }
    return false;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    //printf("%d depth %d %d depth %d\n",s, Depth[s], d, Depth[d]);
    if (Depth[d] == Depth[s] && __sync_bool_compare_and_swap(lock, -1, 0)){
      *CycleEdge_u = s;
      *CycleEdge_v = d;
      return true;
    }
    return false;
  }
  // cond function checks if vertex has been visited yet
  // Only check for cycles amongst already visted vertices.
  inline bool cond(uint32_t d) { return (Parents[d] != -1); }
};

// Honestly I'm not sure if this is terminating because 6 is root
// or because 6 splits the cycle.
struct CYCLE_BT_F {
  int32_t *Parents;
  int32_t *Pair;

  explicit CYCLE_BT_F(int32_t *_Parents, int32_t *_Pair) : 
  Parents(_Parents), Pair(_Pair) {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    if (Parents[s] == d) {
      ++Pair[d];
      return true;
    }
    return false;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    if (Parents[s] == d) {
      // The first call to this will return true, but if they are converging,
      // the second will return false.  True & False == False, so the BT with terminate empty.
      uint32_t nc = __sync_fetch_and_add(&Pair[d], 1);
      return true;
    }
    return false;
  }
  // Only BT to parents while the cycle hasn't converged.
  inline bool cond(uint32_t d) { return true; }
};

// Honestly I'm not sure if this is terminating because 6 is root
// or because 6 splits the cycle.
struct CYCLE_BT_2_F {
  int32_t *Parents;
  int32_t *Pairs;

  explicit CYCLE_BT_2_F(int32_t *_Parents, int32_t *_Pairs) : 
  Parents(_Parents), Pairs(_Pairs) {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    if (Parents[s] == d) {
      Pairs[d] = Pairs[s];
      return true;
    }
    return false;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    if (Parents[s] == d) {
      __sync_bool_compare_and_swap(&Pairs[d], -1, Pairs[s]);
        return true;
    }
    return false;
  }
  // Only BT to parents while the cycle hasn't converged.
  inline bool cond(uint32_t d) { return true; }
};


template <typename T> 
struct GET_XQ_F {
  T *Pair;
  explicit GET_XQ_F(T *_Pair) : 
  Pair(_Pair) {}
  inline bool operator()(uintE v) {
    return Pair[v];
  }
};

template <typename T> 
struct GET_XQ_2_F {
  T *Pair;
  T *Parents;
  explicit GET_XQ_2_F(T *_Pair, T *_Parents) : 
  Pair(_Pair),
  Parents(_Parents) {}
  inline bool operator()(uintE v) {
    //printf("Pair[%d] %d Pair[Parents[%d]] %d\n", v, Pair[v], v, Pair[Parents[v]]);
    return Pair[v] != Pair[Parents[v]];
  }
};

template <typename T> 
struct SET_SOLUTION_H_F {
  T *Solution;
  explicit SET_SOLUTION_H_F(T *_Solution) : 
  Solution(_Solution){}
  inline bool operator()(uintE v) {
    //printf("Pair[%d] %d Pair[Parents[%d]] %d\n", v, Pair[v], v, Pair[Parents[v]]);
    return Solution[v] = 1;
  }
};

// Honestly I'm not sure if this is terminating because 6 is root
// or because 6 splits the cycle.
struct H_SET_CYCLE_F {
  int32_t *Parents;
  int32_t *Pair;
  int32_t *Cycles;
  int32_t xq;

  explicit H_SET_CYCLE_F(int32_t *_Parents, int32_t *_Pair, int32_t *_Cycles, int32_t _xq) : 
  Parents(_Parents), Pair(_Pair), Cycles(_Cycles), xq(_xq) {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    if ((Parents[s] == d || Pair[s] == d) && !Cycles[d]) {
      Cycles[d] = 1;
      return xq != d;
    }
    return false;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    if ((Parents[s] == d || Pair[s] == d) && !Cycles[d]) {
      Cycles[d] = 1;
      return xq != d;
    }
    return false;
  }
  // Only BT to parents while the cycle hasn't converged.
  inline bool cond(uint32_t d) { return true; }
};
