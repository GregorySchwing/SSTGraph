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
  int32_t *Depth;
  explicit H_F(int32_t *_Parents, int32_t *_Depth) : Parents(_Parents), Depth(_Depth) {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    if (Parents[d] == -1) {
      Parents[d] = s;
      Depth[d] = Depth[s]+1;
      return true;
    }
    return false;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    return __sync_bool_compare_and_swap(&Parents[d], -1, s) && 
    __sync_bool_compare_and_swap(&Depth[d], -1, Depth[s]+1);
  }
  // cond function checks if vertex has been visited yet
  inline bool cond(uint32_t d) { return (Parents[d] == -1); }
};

struct I_F {
  int32_t *Parents;
  int32_t *Depth;
  int32_t *match;
  explicit I_F(int32_t *_Parents, int32_t *_Depth, int32_t *_match) : Parents(_Parents), Depth(_Depth), match(_match) {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    if (Parents[d] == -1) {
      Parents[d] = s;
      Depth[d] = Depth[s]+1;
      return true;
    }
    return false;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    return __sync_bool_compare_and_swap(&Parents[d], -1, s) && 
    __sync_bool_compare_and_swap(&Depth[d], -1, Depth[s]+1);
  }
  // cond function checks if vertex has been visited yet
  // also destination should be in M.
  inline bool cond(uint32_t d) { return (Parents[d] == -1 && match[d] != -1); }
};

struct CYCLE_DETECTION_F {
  int32_t *Parents;
  int32_t *Pairs;
  int32_t *Depth;

  explicit CYCLE_DETECTION_F(int32_t *_Parents, int32_t *_Pairs, int32_t *_Depth) : 
  Parents(_Parents), Pairs(_Pairs), Depth(_Depth) {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    printf("%d depth %d %d depth %d\n",s, Depth[s], d, Depth[d]);
    if (Depth[d] == Depth[s]) {
      Pairs[d] = s;
      return true;
    }
    return false;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    return __sync_bool_compare_and_swap(&Depth[d], Depth[s], s);
  }
  // cond function checks if vertex has been visited yet
  // Only check for cycles amongst already visted vertices.
  inline bool cond(uint32_t d) { return (Parents[d] != -1); }
};


struct CYCLE_BT_F {
  int32_t *Parents;
  int32_t *Pairs;

  explicit CYCLE_BT_F(int32_t *_Parents, int32_t *_Pairs) : 
  Parents(_Parents), Pairs(_Pairs) {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    if (Parents[s] == d) {
      return true;
    }
    return false;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    if (Parents[s] == d) {
      return true;
    }
    return false;
  }
  // cond function checks if vertex has been visited yet
  // Only check for cycles amongst already visted vertices.
  inline bool cond(uint32_t d) { return (Parents[d] != -1); }
};
