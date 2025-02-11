#include "BitArray.hpp"
#include "PMA.hpp"
#include "SparseMatrix.hpp"
#include "TinySet.hpp"
#include "TinySet_small.h"
#include "TinySet_small.hpp"
#include "algorithms/BC.h"
#include "algorithms/BFS.h"
#include "algorithms/BellmanFord.h"
#include "algorithms/Components.h"
#include "algorithms/PageRank.h"
#include "algorithms/TC.h"
#include "algorithms/Touchall.h"
#include "algorithms/VC.h"
#include "algorithms/VC_BnB.h"
#include "helpers.h"
#include "parallel.h"
#include "rmat_util.h"
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

constexpr int batch_size = 1000;

std::mt19937 rng;
void srand(uint32_t seed) { rng.seed(seed); }
int random_int() { return rng(); }

uint32_t rand_in_range(uint32_t max) { return random_int() % max; }
double random_float(uint32_t max) {
  return static_cast<double>(random_int() % max);
}

int timing_inserts(uint64_t max_size) {
  printf("std::set, b, 32,");
  uint64_t start = get_usecs();
  std::set<uint32_t> s32;
  uint64_t end = get_usecs();
  printf("creation, %lu, ", end - start);
  start = get_usecs();
  for (uint32_t i = 0; i < max_size; i++) {
    s32.insert(i);
  }
  end = get_usecs();
  printf("insertion, %lu,", end - start);
  start = get_usecs();
  uint32_t sum = 0;
  for (auto el : s32) {
    sum += el;
  }
  end = get_usecs();
  printf("sum_time, %lu, sum_total, %u\n", end - start, sum);

  printf("std::set, b, 16,");
  start = get_usecs();
  std::set<uint16_t> s16;
  end = get_usecs();
  printf("creation, %lu, ", end - start);
  start = get_usecs();
  for (uint32_t i = 0; i < max_size; i++) {
    s16.insert(i);
  }
  end = get_usecs();
  printf("insertion, %lu,", end - start);
  start = get_usecs();
  sum = 0;
  for (auto el : s16) {
    sum += el;
  }
  end = get_usecs();
  printf("sum_time, %lu, sum_total, %u\n", end - start, sum);

  printf("std::unordered_set, b, 32, ");
  start = get_usecs();
  std::unordered_set<uint32_t> us32;
  end = get_usecs();
  printf("creation, %lu, ", end - start);
  start = get_usecs();
  for (uint32_t i = 0; i < max_size; i++) {
    us32.insert(i);
  }
  end = get_usecs();
  printf("insertion, %lu,", end - start);
  start = get_usecs();
  sum = 0;
  for (auto el : us32) {
    sum += el;
  }
  end = get_usecs();
  printf("sum_time, %lu, sum_total, %u\n", end - start, sum);

  printf("std::unordered_set, b, 16, ");
  start = get_usecs();
  std::unordered_set<uint16_t> us16;
  end = get_usecs();
  printf("creation, %lu, ", end - start);
  start = get_usecs();
  for (uint32_t i = 0; i < max_size; i++) {
    us16.insert(i);
  }
  end = get_usecs();
  printf("insertion, %lu,", end - start);
  start = get_usecs();
  sum = 0;
  for (auto el : us16) {
    sum += el;
  }
  end = get_usecs();
  printf("sum_time, %lu, sum_total, %u\n", end - start, sum);

  if (max_size <= 100000) {
    printf("std::vector, b , 32,");
    start = get_usecs();
    std::vector<uint32_t> v;
    end = get_usecs();
    printf("creation, %lu, ", end - start);
    start = get_usecs();
    for (uint32_t i = 0; i < max_size; i++) {
      v.insert(v.begin(), i);
    }
    end = get_usecs();
    printf("insertion, %lu,", end - start);
    start = get_usecs();
    sum = 0;
    for (auto el : v) {
      sum += el;
    }
    end = get_usecs();
    printf("sum_time, %lu, sum_total, %u\n", end - start, sum);
  }
  if (max_size <= 100000) {
    printf("std::vector, b , 16,");
    start = get_usecs();
    std::vector<uint16_t> v;
    end = get_usecs();
    printf("creation, %lu, ", end - start);
    start = get_usecs();
    for (uint32_t i = 0; i < max_size; i++) {
      v.insert(v.begin(), i);
    }
    end = get_usecs();
    printf("insertion, %lu,", end - start);
    start = get_usecs();
    sum = 0;
    for (auto el : v) {
      sum += el;
    }
    end = get_usecs();
    printf("sum_time, %lu, sum_total, %u\n", end - start, sum);
  }
  if (max_size <= 100000) {
    printf("std::vector, b , 8,");
    start = get_usecs();
    std::vector<uint8_t> v;
    end = get_usecs();
    printf("creation, %lu, ", end - start);
    start = get_usecs();
    for (uint32_t i = 0; i < max_size; i++) {
      v.insert(v.begin(), i);
    }
    end = get_usecs();
    printf("insertion, %lu,", end - start);
    start = get_usecs();
    sum = 0;
    for (auto el : v) {
      sum += el;
    }
    end = get_usecs();
    printf("sum_time, %lu, sum_total, %u\n", end - start, sum);
  }
  for (uint8_t b = 8; b <= 32; b += 8) {
    if ((1UL << b) < max_size) {
      continue;
    }
    printf("PMA, b, %u, ", b);
    start = get_usecs();
    PMA<4> pma;
    end = get_usecs();
    printf("creation, %lu, ", end - start);
    start = get_usecs();
    for (uint32_t i = 0; i < max_size; i++) {
      pma.insert(i);
    }
    end = get_usecs();
    printf("insertion, %lu,", end - start);
    start = get_usecs();
    sum = pma.sum_keys();
    end = get_usecs();
    printf("sum_time, %lu, sum_total, %u\n", end - start, sum);
  }

  printf("TinySet");
  start = get_usecs();
  TinySetV ts(max_size);
  end = get_usecs();
  printf("creation, %lu, ", end - start);
  start = get_usecs();
  for (uint32_t i = 0; i < max_size; i++) {
    ts.insert(i);
  }
  end = get_usecs();
  printf("insertion, %lu,", end - start);
  start = get_usecs();
  sum = ts.sum_keys();
  end = get_usecs();
  printf("sum_time, %lu, sum_total, %u\n", end - start, sum);
  // for (uint8_t b = 8; b <= 32; b += 8) {
  //   printf("PackedArray, b, %u,", b);
  //   start = get_usecs();
  //   block_t *array = create_PackedArray<0>(32);
  //   end = get_usecs();
  //   printf("creation, %lu, ", end - start);
  //   start = get_usecs();
  //   for (uint32_t i = 0; i < max_size; i++) {
  //     array = PackedArray_insert<0>(array, i, i, 0);
  //   }
  //   end = get_usecs();
  //   printf("insertion, %lu,", end - start);
  //   start = get_usecs();
  //   sum = 0;
  //   for (uint32_t i = 0; i < max_size; i++) {
  //     sum += PackedArray_get<0>(array, i, b);
  //   }
  //   end = get_usecs();
  //   printf("sum_time, %lu, sum_total, %u\n", end - start, sum);
  // }
  return 0;
}

int timing_random_inserts(uint64_t max_size, uint64_t num_inserts) {
  srand(0);
  printf("std::set, b, 32,");
  uint64_t start = get_usecs();
  std::set<uint32_t> s32;
  uint64_t end = get_usecs();
  printf("creation, %lu, ", end - start);
  start = get_usecs();
  for (uint32_t i = 0; i < num_inserts; i++) {
    s32.insert(rand_in_range(max_size));
  }
  end = get_usecs();
  printf("insertion, %lu,", end - start);
  start = get_usecs();
  uint32_t sum = 0;
  for (auto el : s32) {
    sum += el;
  }
  end = get_usecs();
  printf("sum_time, %lu, sum_total, %u\n", end - start, sum);

  srand(0);
  printf("std::unordered_set, b, 32, ");
  start = get_usecs();
  std::unordered_set<uint32_t> us32;
  end = get_usecs();
  printf("creation, %lu, ", end - start);
  start = get_usecs();
  for (uint32_t i = 0; i < num_inserts; i++) {
    us32.insert(rand_in_range(max_size));
  }
  end = get_usecs();
  printf("insertion, %lu,", end - start);
  start = get_usecs();
  sum = 0;
  for (auto el : us32) {
    sum += el;
  }
  end = get_usecs();
  printf("sum_time, %lu, sum_total, %u\n", end - start, sum);

  srand(0);
  printf("TinySet");
  start = get_usecs();
  TinySetV ts(max_size);
  end = get_usecs();
  printf("creation, %lu, ", end - start);
  start = get_usecs();
  for (uint32_t i = 0; i < num_inserts; i++) {
    ts.insert(rand_in_range(max_size));
  }
  end = get_usecs();
  printf("insertion, %lu,", end - start);
  start = get_usecs();
  sum = ts.sum_keys();
  end = get_usecs();
  printf("sum_time, %lu, sum_total, %u\n", end - start, sum);
  srand(0);
  printf("BitArray, ");
  start = get_usecs();
  BitArray bitarray = BitArray(max_size);
  end = get_usecs();
  printf("creation, %lu, ", end - start);
  start = get_usecs();
  for (uint32_t i = 0; i < num_inserts; i++) {
    bitarray.set(rand_in_range(max_size));
  }
  end = get_usecs();
  printf("insertion, %lu,", end - start);
  start = get_usecs();
  sum = 0;
  for (uint32_t i = 0; i < max_size; i++) {
    if (bitarray.get(i)) {
      sum += i;
    }
  }
  end = get_usecs();
  printf("sum_time, %lu, sum_total, %u\n", end - start, sum);
  return 0;
}

void perf_test_tinyset(uint32_t N) {
  printf("testing tinyset N = %u\n", N);
  uint64_t start = get_usecs();
  TinySetV ts(N);
  for (uint32_t i = 0; i < N; i++) {
    ts.insert(i);
    if (!ts.has(i)) {
      printf("don't have element %u, stopping\n", i);
      return;
    }
  }
  for (uint32_t i = 0; i < N; i++) {
    if (!ts.has(i)) {
      printf("don't have element %u, stopping\n", i);
      return;
    }
  }
  uint64_t duration = get_usecs() - start;
  printf("it took %lu seconds\n", duration / 1000000);
  start = get_usecs();
  uint64_t sum = ts.sum_keys();
  uint64_t duration2 = get_usecs() - start;
  printf("it took %lu ms to sum\n", duration2 / 1000);
  printf("%u, %lu, %lu, %lu\n", N, duration, ts.get_size(), sum);
}

template <int b> void PMA_add_test(uint32_t el_count, bool check = false) {
  uint32_t max = 0;
  if constexpr (b == 4) {
    max = UINT32_MAX;
  } else {
    max = 1U << (b * 8);
  }
  PMA<b> pma;
  el_count = std::min(max - 1, el_count);
  std::vector<uint32_t> random_numbers;
  for (uint32_t i = 0; i < el_count; i++) {
    random_numbers.push_back(random_int() % max);
  }
  uint64_t start = get_usecs();
  for (uint32_t i = 0; i < el_count; i++) {
    // printf("trying to insert %u\n", random_numbers[i]);
    pma.insert(random_numbers[i]);
    if (check) {
      if (!pma.has(random_numbers[i])) {
        printf("FAILED: don't have something we inserted while inserting "
               "elements, %u\n",
               random_numbers[i]);
        pma.print_pma();
        exit(-1);
      }
    }
  }
  for (uint32_t i = 0; i < el_count; i++) {
    if (check) {
      if (!pma.has(random_numbers[i])) {
        printf("FAILED: don't have something we inserted after inserting "
               "elements, "
               "index was %u\n",
               i);
        exit(-1);
      }
    }
  }
  uint64_t insert_duration = get_usecs() - start;
  start = get_usecs();
  uint64_t sum = pma.sum_keys();
  uint64_t sum_duration = get_usecs() - start;
  printf("max, %.9u, b , %.2u, el_count , %.9u, average_insert_time, %.10f, "
         "sum_time, %.6f\n",
         max, b, el_count, ((double)insert_duration) / (1000000 * el_count),
         ((float)sum_duration) / 1000000);
  if (check) {
    std::set<uint32_t> checker(random_numbers.begin(), random_numbers.end());
    uint64_t sum_check = 0;
    for (auto s : checker) {
      sum_check += s;
    }
    if (sum_check != sum) {
      printf("FAILED: got bad result in pma sum in PMA_add_test, got %lu, "
             "expect %lu\n",
             sum, sum_check);
      exit(-1);
    }
  }
}

template <int b> void PMA_remove_test(uint32_t el_count, bool check = false) {
  uint32_t max = 0;
  if constexpr (b == 4) {
    max = UINT32_MAX;
  } else {
    max = 1U << (b * 8);
  }
  PMA<b> pma;
  el_count = std::min(max - 1, el_count);
  std::unordered_set<uint32_t> random_numbers_s;
  while (random_numbers_s.size() < el_count) {
    random_numbers_s.insert(random_int() % max);
  }
  std::vector<uint32_t> random_numbers;
  random_numbers.reserve(el_count);
  for (auto el : random_numbers_s) {
    random_numbers.push_back(el);
  }
  uint64_t start = get_usecs();
  for (uint32_t i = 0; i < el_count; i++) {
    pma.insert(random_numbers[i]);
    // pma.print_pma();
    if (check) {
      if (!pma.has(random_numbers[i])) {
        printf("don't have something we inserted while inserting elements\n");
        exit(-1);
      }
    }
  }
  for (uint32_t i = 0; i < el_count; i++) {
    if (check) {
      if (!pma.has(random_numbers[i])) {
        printf("don't have something we inserted after inserting elements, "
               "index was %u\n",
               i);
        exit(-1);
      }
    }
  }
  uint64_t insert_duration = get_usecs() - start;
  start = get_usecs();
  for (uint32_t i = 0; i < el_count; i++) {
    pma.remove(random_numbers[i]);
    // pma.print_pma();
    if (check) {
      if (pma.has(random_numbers[i])) {
        printf("have something we removed while removing elements, tried to "
               "remove %u\n",
               random_numbers[i]);
        exit(-1);
      }
      for (uint32_t j = i + 1; j < el_count; j++) {
        if (!pma.has(random_numbers[j])) {
          printf("we removed %u when we shouldn't have\n", random_numbers[j]);
          exit(-1);
        }
      }
    }
  }
  if (check && pma.get_n() != 0) {
    printf("still have elements when we shouldn't\n");
    exit(-1);
  }
  uint64_t remove_duration = get_usecs() - start;
  printf("insert duration = %lu, remove duration = %lu\n", insert_duration,
         remove_duration);
}

bool TinySet_add_test(uint32_t max, uint32_t el_count, bool check = false) {
  TinySetV ts(max);
  std::vector<uint32_t> random_numbers;
  for (uint32_t i = 0; i < el_count; i++) {
    uint32_t el = random_int() % max;
    random_numbers.push_back(el);
  }
  uint64_t start = get_usecs();
  for (uint32_t i = 0; i < el_count; i++) {
    ts.insert(random_numbers[i]);
    if (check) {
      // printf("inserting %u\n", random_numbers[i]);
      if (!ts.has(random_numbers[i])) {
        printf("don't have something we inserted while inserting elements\n");
        // ts.print();
        exit(-1);
      }
    }
  }
  for (uint32_t i = 0; i < el_count; i++) {
    if (check) {
      if (!ts.has(random_numbers[i])) {
        printf("don't have something we inserted after inserting elements, "
               "index was %u, number was %u\n",
               i, random_numbers[i]);
        exit(-1);
      }
    }
  }
  uint64_t insert_duration = get_usecs() - start;
  start = get_usecs();
  // printf("pma_sum = %lu\n", ts.pmas[0].pma32.sum());
  uint64_t sum = ts.sum_keys();
  uint64_t sum_duration = get_usecs() - start;
  printf("max, %.10u, el_count , %.9u, average_insert_time milles, %.10f, "
         "sum_time milles, %.6f\n",
         max, el_count, ((double)insert_duration) / (double)(1000UL * el_count),
         ((float)sum_duration) / (double)1000UL);
  if (check) {
    std::set<uint32_t> checker(random_numbers.begin(), random_numbers.end());
    uint64_t sum_check = 0;
    for (auto s : checker) {
      sum_check += s;
    }
    if (sum_check != sum) {
      printf(
          "got bad result in TS sum in TinySet_add_test, got %lu, expect %lu\n",
          sum, sum_check);
      ts.print_pmas();
      ts.print();
      for (auto s : checker) {
        printf("%u, ", s);
      }
      printf("\n");
      exit(-1);
      return false;
    }
  }
  return true;
}
bool TinySet_add_test_fast(uint32_t max, uint32_t el_count) {
  TinySetV ts(max);
  std::vector<uint32_t> random_numbers;
  for (uint32_t i = 0; i < el_count; i++) {
    uint32_t el = random_int() % max;
    random_numbers.push_back(el);
  }
  uint64_t start = get_usecs();
  for (uint32_t i = 0; i < el_count; i++) {
    ts.insert(random_numbers[i]);
  }
  uint64_t insert_duration = get_usecs() - start;
  start = get_usecs();
  uint64_t sum = ts.sum_keys();
  uint64_t sum_duration = get_usecs() - start;
  printf("max, %.10u, el_count , %.9u, fill_frac, %f, "
         "average_insert_time, %.10f, sum_time, %.6f, sum, %lu\n",
         max, el_count, ((double)el_count) / max,
         ((float)insert_duration) / (double)(1000000UL * el_count),
         ((float)sum_duration) / (double)1000000, sum);
  return true;
}

void tinyset_remove_test(uint32_t el_count, uint32_t max_val,
                         bool check = false) {
  TinySetV ts(max_val);
  el_count = std::min(max_val - 1, el_count);
  std::unordered_set<uint32_t> random_numbers_s;
  while (random_numbers_s.size() < el_count) {
    random_numbers_s.insert(random_int() % max_val);
  }
  std::vector<uint32_t> random_numbers;
  random_numbers.reserve(el_count);
  for (auto el : random_numbers_s) {
    random_numbers.push_back(el);
  }
  uint64_t start = get_usecs();
  for (uint32_t i = 0; i < el_count; i++) {
    ts.insert(random_numbers[i]);
    if (check) {
      // printf("adding %u\n", random_numbers[i]);
      if (!ts.has(random_numbers[i])) {
        printf("don't have something we inserted while inserting elements\n");
        exit(-1);
      }
    }
  }
  for (uint32_t i = 0; i < el_count; i++) {
    if (check) {
      if (!ts.has(random_numbers[i])) {
        printf("don't have something we inserted after inserting elements, "
               "index was %u\n",
               i);
        exit(-1);
      }
    }
  }
  uint64_t insert_duration = get_usecs() - start;
  start = get_usecs();
  for (uint32_t i = 0; i < el_count; i++) {
    ts.remove(random_numbers[i]);
    // ts.print_pmas();
    if (check) {
      // printf("removing %u\n", random_numbers[i]);
      if (ts.has(random_numbers[i])) {
        printf("have something we removed while removing elements, tried to "
               "remove %u\n",
               random_numbers[i]);
        exit(-1);
      }
      for (uint32_t j = i + 1; j < el_count; j++) {
        if (!ts.has(random_numbers[j])) {
          printf("we removed %u when we shouldn't have\n", random_numbers[j]);
          exit(-1);
        }
      }
    }
  }
  if (check && ts.get_n() != 0) {
    printf("still have elements when we shouldn't\n");
    exit(-1);
  }
  uint64_t remove_duration = get_usecs() - start;
  printf("insert duration = %lu, remove duration = %lu\n", insert_duration,
         remove_duration);
}

bool real_graph(const std::string &filename, [[maybe_unused]] bool symetric,
                int iters = 20, uint32_t start_node = 0,
                uint32_t max_batch_size = 100000) {
  uint32_t num_nodes = 0;
  uint64_t num_edges = 0;
  std::tuple<el_t, el_t> *edges =
      get_edges_from_file(filename, &num_edges, &num_nodes);

  printf("done reading in the file, n = %u, m = %lu\n", num_nodes, num_edges);
  uint64_t start = get_usecs();
  SparseMatrixV<true, bool> g(num_nodes, num_nodes);
  uint64_t end = get_usecs();
  printf("creation took %lums\n", (end - start) / 1000);
  start = get_usecs();

  uint64_t bfs_milles = 0;
  uint64_t pr_milles = 0;
  uint64_t bc_milles = 0;
  uint64_t cc_milles = 0;
  uint64_t tc_milles = 0;
  uint64_t bf_milles = 0;
  uint64_t vc_milles = 0;
  uint64_t add_batch[10] = {0};
  /*
  for (uint32_t i = 0; i < num_edges; i++) {
    //printf("adding edge %u, (%u, %u)\n", i, srcs[i], dests[i]);
    g.insert(edges[i].x, edges[i].y);
  }
  */

  // uint32_t local_batch_size = 1000;
  uint32_t local_batch_size = batch_size;
  if (num_edges > 10000) {
    local_batch_size = 10000;
  }
  if (num_edges > 100000) {
    local_batch_size = 100000;
  }
  if (num_edges > 10000000) {
    local_batch_size = 10000000;
  }
  if (num_edges > 50000000) {
    local_batch_size = 50000000;
  }
  if (num_edges > 100000000) {
    local_batch_size = 100000000;
  }
  uint64_t i = 0;
  if (num_edges > local_batch_size) {
    for (; i < num_edges - local_batch_size; i += local_batch_size) {
      g.insert_batch(edges + i, local_batch_size);
      // fprintf(stderr, "num_edges added = %lu\n", i + local_batch_size);
    }
  }
  g.insert_batch(edges + i, num_edges % local_batch_size);

  end = get_usecs();

  // g.insert_batch(edges, num_edges);
  free(edges);
  printf("inserting the edges took %lums\n", (end - start) / 1000);
  // uint64_t insert_time = end - start;
  uint64_t size = g.get_memory_size();
  printf("size = %lu bytes, number of edges = %lu, number of nodes = %u\n",
         size, g.M(), num_nodes);
  g.print_statistics();

  SparseMatrixV<true, bool> g2(g);
  g2.print_statistics();
#if 1 
  printf("start vc\n");
  start = get_usecs();
  int32_t *parallel_vc_result = VC_with_edge_map(g);
  end = get_usecs();
  printf("VC: ");
  //for (uint32_t j = 0; j < num_nodes; j++) {
  //  if (parallel_vc_result[j])
  //    printf("%lu ", j);
  //}
  printf("\n");
  bool validSolution = Check_VC_with_edge_map(parallel_vc_result, g2);
  printf("VC is valid solution : %s\n", validSolution ? "true" : "false");
  int32_t vc_count = 0;
  for (int j = 0; j < g.get_rows(); j++) {
    vc_count += parallel_vc_result[j];
  }
  printf("VC size : %u\n", vc_count);
  printf("time to vc %lu micros\n",
         end - start);
#endif

#if 1 
  SparseMatrixV<true, bool> g3(g2);
  printf("start VC_BnB\n");
  start = get_usecs();
  int32_t *parallel_vc_BnB_result = VC_BnB_with_edge_map(g2);
  end = get_usecs();

  printf("time to VC_BnB %lu micros\n",
         end - start);
  bool validSolution_BnB = Check_VC_with_edge_map(parallel_vc_BnB_result, g3);
  printf("VC_BnB is valid solution : %s\n", validSolution_BnB ? "true" : "false");
  int32_t vc_count_BnB = 0;
  for (int j = 0; j < g.get_rows(); j++) {
    vc_count_BnB += parallel_vc_BnB_result[j];
  }
  printf("VC_BnB size : %u\n", vc_count_BnB);
  printf("time to vc %lu micros\n",
         end - start);
#endif

#if 1
  start = get_usecs();
  uint64_t sum1 = g.touch_all_sum();
  end = get_usecs();
  printf("sum of all the edges was = %lu time to count %lu micros\n", sum1,
         end - start);
  start = get_usecs();
  uint64_t sum2 = TouchAll(g);
  end = get_usecs();
  printf("sum of all the edges was = %lu time to count %lu micros\n", sum2,
         end - start);

  int32_t parallel_bfs_result2_ = 0;
  uint64_t parallel_bfs_time2 = 0;

  for (int i = 0; i < iters; i++) {
    start = get_usecs();
    int32_t *parallel_bfs_result = BFS_with_edge_map(g, start_node);
    end = get_usecs();
    parallel_bfs_result2_ += parallel_bfs_result[0];
    if (i == 0 && parallel_bfs_result != nullptr) {
      uint64_t reached = 0;
      for (uint32_t j = 0; j < num_nodes; j++) {
        reached += parallel_bfs_result[j] != -1;
      }
      printf("the bfs from source %u, reached %lu vertices\n", start_node,
             reached);
    }
    std::vector<uint32_t> depths(num_nodes, UINT32_MAX);
    parallel_for(uint32_t j = 0; j < num_nodes; j++) {
      uint32_t current_depth = 0;
      int32_t current_parent = j;
      if (parallel_bfs_result[j] < 0) {
        continue;
      }
      while (current_parent != parallel_bfs_result[current_parent]) {
        current_depth += 1;
        current_parent = parallel_bfs_result[current_parent];
      }
      depths[j] = current_depth;
    }
    std::ofstream myfile;
    myfile.open("bfs.out");
    for (unsigned int i = 0; i < num_nodes; i++) {
      myfile << depths[i] << "\n";
    }
    myfile.close();

    free(parallel_bfs_result);
    parallel_bfs_time2 += (end - start);
  }
  // printf("bfs took %lums, parent of 0 = %d\n", (bfs_time)/(1000*iters),
  // bfs_result_/iters);
  printf("parallel_bfs with edge_map took %lums, parent of 0 = %d\n",
         parallel_bfs_time2 / (1000 * iters), parallel_bfs_result2_ / iters);
  bfs_milles = parallel_bfs_time2 / (1000 * iters);

  start = get_usecs();
  auto *values3 = PR_S<double>(g, 10);
  end = get_usecs();
  printf("pagerank with MAPS took %lums, value of 0 = %f, for %d iters\n",
         (end - start) / (1000), values3[0], iters);
  pr_milles = (end - start) / (1000);
  std::ofstream myfile;
  myfile.open("pr.out");
  for (unsigned int i = 0; i < num_nodes; i++) {
    myfile << values3[i] << "\n";
  }
  myfile.close();
  free(values3);

  start = get_usecs();
  double *values4 = nullptr;
  double dep_0 = 0;
  for (int i = 0; i < iters; i++) {
    if (values4 != nullptr) {
      free(values4);
    }
    values4 = BC(g, start_node);
    dep_0 += values4[0];
  }
  end = get_usecs();
  printf("BC took %lums, value of 0 = %f\n", (end - start) / (1000 * iters),
         dep_0 / iters);
  bc_milles = (end - start) / (1000 * iters);
  if (values4 != nullptr) {
    std::ofstream myfile;
    myfile.open("bc.out");
    for (uint32_t i = 0; i < num_nodes; i++) {
      myfile << values4[i] << "\n";
    }
    myfile.close();
    free(values4);
  }

  start = get_usecs();
  uint32_t *values5 = nullptr;
  uint32_t id_0 = 0;
  for (int i = 0; i < iters; i++) {
    if (values5) {
      free(values5);
    }
    values5 = CC(g);
    id_0 += values5[0];
  }
  end = get_usecs();
  printf("CC took %lums, value of 0 = %u\n", (end - start) / (1000 * iters),
         id_0 / iters);
  cc_milles = (end - start) / (1000 * iters);
  if (values5 != nullptr) {
    std::unordered_map<uint32_t, uint32_t> components;
    for (uint32_t i = 0; i < num_nodes; i++) {
      components[values5[i]] += 1;
    }
    printf("there are %zu components\n", components.size());
    uint32_t curent_max = 0;
    uint32_t curent_max_key = 0;
    for (auto p : components) {
      if (p.second > curent_max) {
        curent_max = p.second;
        curent_max_key = p.first;
      }
    }
    printf("the element with the biggest component is %u, it has %u members "
           "to its component\n",
           curent_max_key, curent_max);
    std::ofstream myfile;
    myfile.open("cc.out");
    for (uint32_t i = 0; i < num_nodes; i++) {
      myfile << values5[i] << "\n";
    }
    myfile.close();
  }

  free(values5);
#if 0
    start = get_usecs();
    TC(g);
    end = get_usecs();
    printf("TC took %lums\n", (end - start) / (1000));
    tc_milles = (end - start) / (1000);
#endif
  start = get_usecs();
  int32_t *bf_values = nullptr;
  int32_t val_0 = 0;
  for (int i = 0; i < iters; i++) {
    if (bf_values != nullptr) {
      free(bf_values);
    }
    bf_values = BF(g, start_node);
    val_0 += bf_values[0];
  }
  end = get_usecs();
  printf("BF took %lums, value of 0 = %d\n", (end - start) / (1000 * iters),
         val_0 / iters);
  bf_milles = (end - start) / (1000 * iters);
  if (bf_values != nullptr) {
    std::ofstream myfile;
    myfile.open("bf.out");
    for (uint32_t i = 0; i < num_nodes; i++) {
      myfile << bf_values[i] << "\n";
    }
    myfile.close();
    free(bf_values);
  }
#endif

// batch updates
#if 1
  auto r = random_aspen();
  uint32_t counter = 0;
  for (uint32_t b_size = 10; b_size <= max_batch_size; b_size *= 10) {
    double batch_insert_time = 0;
    double batch_remove_time = 0;
    for (int it = 0; it < iters; it++) {
      // uint64_t size = g.get_memory_size();
      // printf("size start = %lu\n", size);
      double a = 0.5;
      double b = 0.1;
      double c = 0.1;
      size_t nn = 1UL << (log2_up(num_nodes) - 1);
      auto rmat = rMat<uint32_t>(nn, r.ith_rand(0), a, b, c);
      std::vector<std::tuple<el_t, el_t>> es(b_size);
      parallel_for(uint32_t i = 0; i < b_size; i++) {
        std::pair<uint32_t, uint32_t> edge = rmat(i);
        es[i] = {edge.first, edge.second};
      }
      // std::unordered_map<uint32_t, uint32_t> count_per_vertex;
      // std::map<uint32_t, uint32_t> count_per_count;
      // for (auto x : es) {
      //   count_per_vertex[x.x]++;
      // }
      // for (auto x : count_per_vertex) {
      //   count_per_count[x.second]++;
      // }
      // for (auto x : count_per_count) {
      //   std::cout << x.second << " vertices had " << x.first << "
      //   elements"
      //             << std::endl;
      // }
      start = get_usecs();
      g.insert_batch(es.data(), b_size);
      end = get_usecs();
      batch_insert_time += (double)(end - start);
      // size = g.get_memory_size();
      // printf("size end = %lu\n", size);
      start = get_usecs();
      g.remove_batch(es.data(), b_size);
      end = get_usecs();
      batch_remove_time += (double)(end - start);
    }
    batch_insert_time /= (1000000 * iters);
    batch_remove_time /= (1000000 * iters);
    printf("batch_size = %d, time to insert = %f seconds, throughput = %4.2e "
           "updates/second\n",
           b_size, batch_insert_time, b_size / (batch_insert_time));
    printf("batch_size = %d, time to remove = %f seconds, throughput = %4.2e "
           "updates/second\n",
           b_size, batch_remove_time, b_size / (batch_remove_time));
    add_batch[counter] = b_size / (batch_insert_time * 1000000);
    counter += 1;
  }
#endif
  printf("%s, %lu, %lu, %lu, %lu, %lu, %lu, %lu, %lu, %lu, %lu %lu\n",
         filename.c_str(), bfs_milles, pr_milles, bc_milles, cc_milles,
         bf_milles, tc_milles, add_batch[0], add_batch[1], add_batch[2],
         add_batch[3], add_batch[4]);

  return true;
}

template <typename value_type>
bool real_graph_weights(const std::string &filename,
                        [[maybe_unused]] bool symetric, int iters = 20,
                        uint32_t start_node = 0,
                        uint32_t max_batch_size = 100000) {
  std::cout << "value_type is " << TypeName<value_type>::Get() << std::endl;
  using edge_type =
      typename std::conditional<std::is_same<value_type, bool>::value,
                                std::tuple<el_t, el_t>,
                                std::tuple<el_t, el_t, value_type>>::type;
  uint32_t num_nodes = 0;
  uint64_t num_edges = 0;
  edge_type *edges =
      get_edges_from_file<value_type>(filename, &num_edges, &num_nodes);

  printf("done reading in the file, n = %u, m = %lu\n", num_nodes, num_edges);
  uint64_t start = get_usecs();
  SparseMatrixV<true, value_type> g(num_nodes, num_nodes);
  uint64_t end = get_usecs();
  printf("creation took %lums\n", (end - start) / 1000);
  start = get_usecs();

  uint32_t local_batch_size = batch_size;
  if (num_edges > 10000) {
    local_batch_size = 10000;
  }
  if (num_edges > 100000) {
    local_batch_size = 100000;
  }
  if (num_edges > 10000000) {
    local_batch_size = 10000000;
  }
  if (num_edges > 50000000) {
    local_batch_size = 50000000;
  }
  if (num_edges > 100000000) {
    local_batch_size = 100000000;
  }
  if (num_edges > 500000000) {
    local_batch_size = 500000000;
  }
  uint64_t i = 0;
  if (num_edges > local_batch_size) {
    for (; i < num_edges - local_batch_size; i += local_batch_size) {
      g.insert_batch(edges + i, local_batch_size);
      fprintf(stderr, "num_edges added = %lu\n", i + local_batch_size);
    }
  }
  g.insert_batch(edges + i, num_edges % local_batch_size);

  end = get_usecs();
  // g.insert_batch(edges, num_edges);
  free(edges);
  printf("inserting the edges took %lums\n", (end - start) / 1000);
  // uint64_t insert_time = end - start;
  uint64_t size = g.get_memory_size();
  printf("size = %lu bytes, number of edges = %lu, number of nodes = %u\n",
         size, g.M(), num_nodes);
  // g.print_statistics();

#if 1
  start = get_usecs();
  uint64_t sum1 = 0;
  for (int i = 0; i < iters; i++) {
    sum1 += g.touch_all_sum();
  }
  end = get_usecs();
  printf("sum of all the edges was = %lu time to count %lu micros\n",
         sum1 / iters, (end - start) / iters);
  start = get_usecs();
  uint64_t sum2 = 0;
  for (int i = 0; i < iters; i++) {
    sum2 += TouchAll(g);
  }
  end = get_usecs();
  printf("sum of all the edges was = %lu time to count %lu micros\n",
         sum2 / iters, (end - start) / iters);
#endif
#if 1
  start = get_usecs();
  int32_t *bf_values = nullptr;
  int32_t val_0 = 0;
  for (int i = 0; i < iters; i++) {
    if (bf_values != nullptr) {
      free(bf_values);
    }
    bf_values = BF(g, start_node);
    val_0 += bf_values[0];
  }
  end = get_usecs();
  printf("BF took %lums, value of 0 = %d\n", (end - start) / (1000 * iters),
         val_0 / iters);
  if (bf_values != nullptr) {
    std::ofstream myfile;
    myfile.open("bf_" + TypeName<value_type>::Get() + ".out");
    for (uint32_t i = 0; i < num_nodes; i++) {
      myfile << bf_values[i] << "\n";
    }
    myfile.close();
    free(bf_values);
  }
#endif
// batch updates
#if 1
  auto r = random_aspen();
  for (uint32_t b_size = 10; b_size <= max_batch_size; b_size *= 10) {
    double batch_insert_time = 0;
    double batch_remove_time = 0;
    for (int it = 0; it < iters; it++) {
      // uint64_t size = g.get_memory_size();
      // printf("size start = %lu\n", size);
      double a = 0.5;
      double b = 0.1;
      double c = 0.1;
      size_t nn = 1UL << (log2_up(num_nodes) - 1);
      auto rmat = rMat<uint32_t>(nn, r.ith_rand(0), a, b, c);
      std::vector<edge_type> es(b_size);
      parallel_for(uint32_t i = 0; i < b_size; i++) {
        std::pair<uint32_t, uint32_t> edge = rmat(i);
        if constexpr (std::is_same_v<value_type, bool>) {
          es[i] = {edge.first, edge.second};
        } else {
          es[i] = {edge.first, edge.second, static_cast<value_type>(1)};
        }
      }
      start = get_usecs();
      g.insert_batch(es.data(), b_size);
      end = get_usecs();
      batch_insert_time += (double)(end - start);
      // size = g.get_memory_size();
      // printf("size end = %lu\n", size);
      start = get_usecs();
      g.remove_batch(es.data(), b_size);
      end = get_usecs();
      batch_remove_time += (double)(end - start);
    }
    batch_insert_time /= (1000000 * iters);
    batch_remove_time /= (1000000 * iters);
    printf("batch_size = %d, time to insert = %f seconds, throughput = %4.2e "
           "updates/second\n",
           b_size, batch_insert_time, b_size / (batch_insert_time));
    printf("batch_size = %d, time to remove = %f seconds, throughput = %4.2e "
           "updates/second\n",
           b_size, batch_remove_time, b_size / (batch_remove_time));
  }
#endif
  return true;
}

void get_graph_distribution(const std::string &filename) {
  uint32_t num_nodes = 0;
  uint64_t num_edges = 0;
  std::tuple<el_t, el_t> *edges =
      get_edges_from_file(filename, &num_edges, &num_nodes);

  SparseMatrixV<true, bool> g(num_nodes, num_nodes);

  uint32_t local_batch_size = batch_size;
  if (num_edges > 10000000) {
    local_batch_size = 10000000;
  }
  if (num_edges > 100000000) {
    local_batch_size = 100000000;
  }
  if (num_edges > 500000000) {
    local_batch_size = 500000000;
  }
  uint64_t i = 0;
  if (num_edges > local_batch_size) {
    for (; i < num_edges - local_batch_size; i += local_batch_size) {
      g.insert_batch(edges + i, local_batch_size);
      fprintf(stderr, "num_edges added = %lu\n", i + local_batch_size);
    }
  }
  g.insert_batch(edges + i, num_edges % local_batch_size);
  std::map<uint64_t, uint64_t> degrees;
  for (uint64_t i = 0; i < num_nodes; i++) {
    degrees[g.getDegree(i)] += 1;
  }
  for (auto pair : degrees) {
    printf("%lu, %lu\n", pair.first, pair.second);
  }
}

bool real_graph_static_test(const std::string &filename,
                            [[maybe_unused]] bool symetric, int iters = 10,
                            uint32_t start_node = 0,
                            const std::string &run_info = "") {
  uint32_t num_nodes = 0;
  uint64_t num_edges = 0;
  std::tuple<el_t, el_t> *edges =
      get_edges_from_file(filename, &num_edges, &num_nodes);

  if (num_nodes == 0) {
    printf("graphs needs to have non zero number of nodes\n");
    free(edges);
    return false;
  }

  SparseMatrixV<true, bool> g(num_nodes, num_nodes);

  uint32_t local_batch_size = batch_size;
  if (num_edges > 10000000) {
    local_batch_size = 10000000;
  }
  if (num_edges > 100000000) {
    local_batch_size = 100000000;
  }
  if (num_edges > 500000000) {
    local_batch_size = 500000000;
  }
  uint64_t i = 0;
  if (num_edges > local_batch_size) {
    for (; i < num_edges - local_batch_size; i += local_batch_size) {
      g.insert_batch(edges + i, local_batch_size);
      fprintf(stderr, "num_edges added = %lu\n", i + local_batch_size);
    }
  }
  g.insert_batch(edges + i, num_edges % local_batch_size);

  // g.insert_batch(edges, num_edges);
  free(edges);
  uint64_t start = 0;
  uint64_t end = 0;
  if (g.get_rows() == 0) {
    printf("graph has no vertices\n");
    return false;
  }
  
#if 1
  int32_t *parallel_vc_result = VC_with_edge_map(g);
  int32_t vc_count = 0;
  for (int j = 0; j < g.get_rows(); j++) {
    vc_count += parallel_vc_result[j];
  }
  printf("VC size : %u\n", vc_count);
  free(parallel_vc_result);
  start = get_usecs();
  for (int i = 0; i < iters; i++) {
    int32_t *parallel_vc_result = VC_with_edge_map(g);
    int32_t vc_count = 0;
    for (int j = 0; j < g.get_rows(); j++) {
      vc_count += parallel_vc_result[j];
    }
    printf("VC size : %u\n", vc_count);
    free(parallel_vc_result);
  }
  end = get_usecs();
  printf("tinyset, %d, VC, %d, %s, %s, %f\n", iters, start_node,
          filename.c_str(), run_info.c_str(),
          ((double)(end - start)) / (1000000 * iters));
  fprintf(stderr, "VC done\n");
#endif

#if 1
  auto *values3 = PR_S<double>(g, 10);
  free(values3);
  start = get_usecs();
  for (int i = 0; i < iters; i++) {
    auto *values3 = PR_S<double>(g, 10);
    free(values3);
  }
  end = get_usecs();
  printf("tinyset, %d, PageRank, %d, %s, %s, %f\n", iters, start_node,
         filename.c_str(), run_info.c_str(),
         ((double)(end - start)) / (1000000 * iters));
  fprintf(stderr, "PR done\n");

  int32_t *parallel_bfs_result = BFS_with_edge_map(g, start_node);
  free(parallel_bfs_result);
  start = get_usecs();
  for (int i = 0; i < iters; i++) {
    int32_t *parallel_bfs_result = BFS_with_edge_map(g, start_node);
    free(parallel_bfs_result);
  }
  end = get_usecs();
  printf("tinyset, %d, BFS, %d, %s, %s, %f\n", iters, start_node,
         filename.c_str(), run_info.c_str(),
         ((double)(end - start)) / (1000000 * iters));
  fprintf(stderr, "BFS done\n");

  double *values4 = BC(g, start_node);
  free(values4);
  start = get_usecs();
  for (int i = 0; i < iters; i++) {
    values4 = BC(g, start_node);
    free(values4);
  }
  end = get_usecs();
  printf("tinyset, %d, BC, %d, %s, %s, %f\n", iters, start_node,
         filename.c_str(), run_info.c_str(),
         ((double)(end - start)) / (1000000 * iters));
  fprintf(stderr, "BC done\n");
  uint32_t *values5 = CC(g);
  start = get_usecs();
  for (int i = 0; i < iters; i++) {
    if (values5) {
      free(values5);
    }
    values5 = CC(g);
  }
  end = get_usecs();
  printf("tinyset, %d, Components, %d, %s, %s, %f\n", iters, start_node,
         filename.c_str(), run_info.c_str(),
         ((double)(end - start)) / (1000000 * iters));
  fprintf(stderr, "CC done\n");
  free(values5);
#endif
#if 1
  TC(g);
  start = get_usecs();
  for (int i = 0; i < iters; i++) {
    TC(g);
  }
  end = get_usecs();
  printf("tinyset, %d, TC, %d, %s, %s, %f\n", iters, start_node,
         filename.c_str(), run_info.c_str(),
         ((double)(end - start)) / (1000000 * iters));
#endif
  return true;
}

void rewrite_graph(const std::string &filename) {
  printf("rewriting graph %s\n", filename.c_str());
  uint32_t num_nodes = 0;
  uint64_t num_edges = 0;
  std::tuple<el_t, el_t> *edges =
      get_edges_from_file(filename, &num_edges, &num_nodes);
  printf("num_nodes = %u\n", num_nodes);
  std::vector<uint32_t> new_node_ids(num_nodes, 0);
  for (uint32_t i = 0; i < num_nodes; i++) {
    new_node_ids[i] = i;
  }
  std::shuffle(new_node_ids.begin(), new_node_ids.end(), rng);
  printf("node 35 in the old graph is node %u in the new\n", new_node_ids[35]);
  std::string f_name = filename + "el.shuf";
  FILE *fw = fopen(f_name.c_str(), "w");
  if (fw == nullptr) {
    printf("file was not opened\n");
    free(edges);
    return;
  }
  for (uint64_t i = 0; i < num_edges; i++) {
    fprintf(fw, "%u   %u\n", new_node_ids[std::get<0>(edges[i])],
            new_node_ids[std::get<1>(edges[i])]);
  }
  free(edges);
  // return 0;
  fclose(fw);
  printf("finished writing %s\n", (filename + "el.shuf").c_str());
}

template <int index_size, typename value_type>
void PMA_map_insert_test_templated(uint32_t el_count, bool check = false) {
  uint32_t max_key = UINT32_MAX;
  if constexpr (index_size < 4) {
    max_key = 1UL << (index_size * 8U);
  }
  value_type max_val = std::numeric_limits<value_type>::max();
  PMA<index_size, value_type> pma;
  el_count = std::min(max_key - 1, el_count);
  std::unordered_map<uint32_t, value_type> random_pairs;

  uint64_t start = get_usecs();
  for (uint32_t i = 0; i < el_count; i++) {
    uint32_t key = random_int() % max_key;
    value_type value = random_int();
    if constexpr (std::is_same<bool, value_type>::value) {
      value = true;
    }
    if (value > max_val) {
      if constexpr (std::is_integral<value_type>::value) {
        value = value % max_val;
      } else {
        value = max_val;
      }
    }
    random_pairs[key] = value;
    // std::cout << "trying to insert (" << key << ", " << +value << ")"
    //           << std::endl;

    pma.insert({key, value});
    if (check) {
      if (!pma.has(key)) {
        printf("FAILED: don't have something we inserted while inserting "
               "elements, %u\n",
               key);
        pma.print_pma();
        exit(-1);
      }
      if (pma.value(key) != value) {
        printf("FAILED: value doesn't match while inserting elements, key was "
               "%u\n",
               key);
        std::cout << +pma.value(key) << ", " << +value << std::endl;
        pma.print_pma();
        exit(-1);
      }
    }
  }
  for (auto &pair : random_pairs) {
    if (check) {
      if (!pma.has(pair.first)) {
        printf("FAILED: don't have something we inserted after inserting "
               "elements, "
               "key was %u\n",
               pair.first);
        pma.print_pma();
        exit(-1);
      }
      if (pma.value(pair.first) != pair.second) {
        printf("FAILED: value doesn't match after inserting elements\n");
        std::cout << "key is " << +pair.first << std::endl;
        std::cout << +pma.value(pair.first) << ", " << +pair.second
                  << std::endl;
        ;
        pma.print_pma();
        exit(-1);
      }
    }
  }
  uint64_t insert_duration = get_usecs() - start;
  printf("max, %.9u, index_size , %.2u, value_size %.2lu, el_count , %.9u, "
         "average_insert_time, "
         "%.10f\n",
         max_key, index_size, sizeof(value_type), el_count,
         ((float)insert_duration) / (float)(1000000 * el_count));
}

void PMA_map_insert_test(uint32_t el_count, bool check = false) {
  PMA_map_insert_test_templated<1, uint8_t>(el_count, check);
  PMA_map_insert_test_templated<2, uint32_t>(el_count, check);
  PMA_map_insert_test_templated<3, float>(el_count, check);
  PMA_map_insert_test_templated<4, double>(el_count, check);
}

template <int index_size, typename value_type>
void PMA_map_remove_test_templated(uint32_t el_count, bool check = false) {
  uint32_t max_key = UINT32_MAX;
  if constexpr (index_size < 4) {
    max_key = 1UL << (index_size * 8U);
  }
  value_type max_val = std::numeric_limits<value_type>::max();
  PMA<index_size, value_type> pma;
  el_count = std::min(max_key - 1, el_count);
  std::unordered_map<uint32_t, value_type> random_pairs = {};

  uint64_t start = get_usecs();
  for (uint32_t i = 0; i < el_count; i++) {
    uint32_t key = random_int() % max_key;
    value_type value = random_int();
    if constexpr (std::is_same<bool, value_type>::value) {
      value = true;
    }
    if (value > max_val) {
      if constexpr (std::is_integral<value_type>::value) {
        value = value % max_val;
      } else {
        value = max_val;
      }
    }
    random_pairs.insert_or_assign(key, value);
    // std::cout << "trying to insert (" << key << ", " << +value << ")"
    //           << std::endl;

    pma.insert({key, value});
    // pma.print_pma();
    if (check) {
      if (!pma.has(key)) {
        printf("FAILED: don't have something we inserted while inserting "
               "elements, %u\n",
               key);
        pma.print_pma();
        exit(-1);
      }
      if (pma.value(key) != value) {
        printf("FAILED: value doesn't match while inserting elements\n");
        std::cout << pma.value(key) << ", " << value << std::endl;
        pma.print_pma();
        exit(-1);
      }
    }
  }
  for (auto &pair : random_pairs) {
    if (check) {
      if (!pma.has(pair.first)) {
        printf("FAILED: don't have something we inserted after inserting "
               "elements, "
               "key was %u\n",
               pair.first);
        exit(-1);
      }
      if (pma.value(pair.first) != pair.second) {
        printf("FAILED: value doesn't match after inserting elements\n");
        std::cout << "key is " << +pair.first << std::endl;
        std::cout << +pma.value(pair.first) << ", " << +pair.second
                  << std::endl;
        ;
        pma.print_pma();
        exit(-1);
      }
    }
  }
  uint64_t insert_duration = get_usecs() - start;
  start = get_usecs();
  for (auto &pair : random_pairs) {
    // std::cout << "trying to remove: " << pair.first << std::endl;
    pma.remove(pair.first);
    // pma.print_pma();
    if (check) {
      if (pma.has(pair.first)) {
        printf("have something we removed while removing elements, tried to "
               "remove %u\n",
               pair.first);
        exit(-1);
      }
    }
  }
  if (check && pma.get_n() != 0) {
    printf("still have elements when we shouldn't\n");
    exit(-1);
  }
  uint64_t remove_duration = get_usecs() - start;
  printf("insert duration = %lu, remove duration = %lu\n", insert_duration,
         remove_duration);
}

void PMA_map_remove_test(uint32_t el_count, bool check = false) {
  PMA_map_remove_test_templated<1, uint8_t>(el_count, check);
  PMA_map_remove_test_templated<2, uint32_t>(el_count, check);
  PMA_map_remove_test_templated<3, float>(el_count, check);
  PMA_map_remove_test_templated<4, double>(el_count, check);
}

template <typename value_type>
void tinyset_map_insert_test_templated(uint32_t el_count, uint32_t max_key,
                                       bool check = false) {
  value_type max_val = std::numeric_limits<value_type>::max();
  TinySetV<value_type> ts(max_key);
  std::unordered_map<uint32_t, value_type> random_pairs;

  uint64_t start = get_usecs();
  for (uint32_t i = 0; i < el_count; i++) {
    uint32_t key = random_int() % max_key;
    value_type value = random_int();
    if constexpr (std::is_same<bool, value_type>::value) {
      value = true;
    }
    if (value > max_val) {
      if constexpr (std::is_integral<value_type>::value) {
        value = value % max_val;
      } else {
        value = max_val;
      }
    }
    random_pairs[key] = value;
    // std::cout << "trying to insert (" << key << ", " << +value << ")"
    //           << std::endl;

    ts.insert(key, value);
    // ts.print_pmas();
    if (check) {
      if (!ts.has(key)) {
        printf("FAILED: don't have something we inserted while inserting "
               "elements, %u\n",
               key);
        ts.print();
        ts.print_pmas();
        exit(-1);
      }
      if (ts.value(key) != value) {
        printf("FAILED: value doesn't match while inserting elements\n");
        std::cout << +ts.value(key) << ", " << value << std::endl;
        ts.print();
        ts.print_pmas();
        exit(-1);
      }
    }
  }
  for (auto &pair : random_pairs) {
    if (check) {
      if (!ts.has(pair.first)) {
        printf("FAILED: don't have something we inserted after inserting "
               "elements, "
               "key was %u\n",
               pair.first);
        exit(-1);
      }
      if (ts.value(pair.first) != pair.second) {
        printf("FAILED: value doesn't match after inserting elements\n");
        std::cout << "key is " << +pair.first << std::endl;
        std::cout << +ts.value(pair.first) << ", " << +pair.second << std::endl;
        ;
        ts.print();
        exit(-1);
      }
    }
  }
  uint64_t insert_duration = get_usecs() - start;
  start = get_usecs();
  uint64_t sum_keys = ts.sum_keys();
  value_type sum_values = ts.sum_values();
  printf("sum_keys = %lu: ", sum_keys);
  std::cout << "sum_keys = " << +sum_values << ": ";
  uint64_t sum_duration = get_usecs() - start;
  if (check) {
    uint64_t correct_sum_keys = 0;
    value_type correct_sum_values = 0;
    for (auto &pair : random_pairs) {
      correct_sum_keys += pair.first;
      correct_sum_values += pair.second;
    }
    if (sum_keys != correct_sum_keys) {
      printf("\nFAILED: sum_keys doesn't match, got %lu, expected %lu\n",
             sum_keys, correct_sum_keys);
      exit(-1);
    }
    if (!approximatelyEqual(sum_values, correct_sum_values,
                            std::numeric_limits<float>::epsilon() * 10000)) {
      printf("\nFAILED: sum_values doesn't match\n");
      std::cout << "got " << +sum_values << " expected " << +correct_sum_values
                << std::endl;
      // ts.print();
      // for (auto &pair : random_pairs) {
      //   std::cout << "{" << pair.first << ", " << +pair.second << "}"
      //             << ", ";
      // }
      // std::cout << std::endl;
      exit(-1);
    }
  }

  printf("value_size %.2lu, el_count , %.9u, "
         "average_insert_time, "
         "%.10f, sum time = %lu\n",
         sizeof(value_type), el_count,
         ((float)insert_duration) / (float)(1000000 * el_count),
         sum_duration / 1000);
}

void tinyset_map_insert_test(
    uint32_t el_count, uint32_t max_key = std::numeric_limits<uint32_t>::max(),
    bool check = false) {
  tinyset_map_insert_test_templated<uint8_t>(el_count, max_key, check);
  tinyset_map_insert_test_templated<uint32_t>(el_count, max_key, check);
  tinyset_map_insert_test_templated<float>(el_count, max_key, check);
  tinyset_map_insert_test_templated<double>(el_count, max_key, check);
}

template <typename value_type>
void tinyset_map_remove_test_templated(uint32_t el_count, uint32_t max_key,
                                       bool check = false) {
  value_type max_val = std::numeric_limits<value_type>::max();
  TinySetV<value_type> ts(max_key);
  std::unordered_map<uint32_t, value_type> random_pairs;

  uint64_t start = get_usecs();
  for (uint32_t i = 0; i < el_count; i++) {
    uint32_t key = random_int() % max_key;
    value_type value = random_int();
    if constexpr (std::is_same<bool, value_type>::value) {
      value = true;
    }
    if (value > max_val) {
      if constexpr (std::is_integral<value_type>::value) {
        value = value % max_val;
      } else {
        value = max_val;
      }
    }
    random_pairs[key] = value;
    // std::cout << "trying to insert (" << key << ", " << +value << ")"
    //           << std::endl;

    ts.insert(key, value);
    if (check) {
      if (!ts.has(key)) {
        printf("FAILED: don't have something we inserted while inserting "
               "elements, %u\n",
               key);
        ts.print();
        exit(-1);
      }
      if (ts.value(key) != value) {
        printf("FAILED: value doesn't match while inserting elements\n");
        std::cout << +ts.value(key) << ", " << value << std::endl;
        ts.print();
        exit(-1);
      }
    }
  }
  for (auto &pair : random_pairs) {
    if (check) {
      if (!ts.has(pair.first)) {
        printf("FAILED: don't have something we inserted after inserting "
               "elements, "
               "key was %u\n",
               pair.first);
        exit(-1);
      }
      if (ts.value(pair.first) != pair.second) {
        printf("FAILED: value doesn't match after inserting elements\n");
        std::cout << "key is " << +pair.first << std::endl;
        std::cout << +ts.value(pair.first) << ", " << +pair.second << std::endl;
        ;
        ts.print();
        exit(-1);
      }
    }
  }
  uint64_t insert_duration = get_usecs() - start;
  start = get_usecs();
  for (auto &pair : random_pairs) {
    ts.remove(pair.first);
    if (check) {
      if (ts.has(pair.first)) {
        printf("have something we removed while removing elements, tried to "
               "remove %u\n",
               pair.first);
        exit(-1);
      }
    }
  }
  if (check && ts.get_n() != 0) {
    printf("still have elements when we shouldn't\n");
    exit(-1);
  }
  uint64_t remove_duration = get_usecs() - start;
  printf("insert duration = %lu, remove duration = %lu\n", insert_duration,
         remove_duration);
}

void tinyset_map_remove_test(
    uint32_t el_count, uint32_t max_key = std::numeric_limits<uint32_t>::max(),
    bool check = false) {
  tinyset_map_remove_test_templated<uint8_t>(el_count, max_key, check);
  tinyset_map_remove_test_templated<uint32_t>(el_count, max_key, check);
  tinyset_map_remove_test_templated<float>(el_count, max_key, check);
  tinyset_map_remove_test_templated<double>(el_count, max_key, check);
}

template <typename value_type>
void matrix_values_add_remove_test_templated(uint32_t el_count,
                                             uint32_t row_count,
                                             bool check = false) {
  value_type max_val = std::numeric_limits<value_type>::max();
  SparseMatrixV<true, value_type> mat(row_count, row_count);
  std::unordered_map<uint32_t, std::unordered_map<uint32_t, value_type>>
      correct_matrix;

  uint64_t start = get_usecs();
  for (uint32_t i = 0; i < el_count; i++) {
    uint32_t r = random_int() % row_count;
    uint32_t c = random_int() % row_count;
    value_type value = random_int();
    if constexpr (std::is_same<bool, value_type>::value) {
      value = true;
    }
    if (value > max_val) {
      if constexpr (std::is_integral<value_type>::value) {
        value = value % max_val;
      } else {
        value = max_val;
      }
    }
    correct_matrix[r][c] = value;

    mat.insert(r, c, value);
    if (check) {
      if (!mat.has(r, c)) {
        printf("FAILED: don't have something we inserted while inserting "
               "elements, %u, %u\n",
               r, c);
        mat.print_arrays();
        exit(-1);
      }
      if (mat.value(r, c) != value) {
        printf("FAILED: value doesn't match while inserting elements\n");
        std::cout << +mat.value(r, c) << ", " << value << std::endl;
        mat.print_arrays();
        exit(-1);
      }
    }
  }
  uint64_t correct_sum = 0;
  for (auto &row : correct_matrix) {
    for (auto &pair : row.second) {
      correct_sum += pair.first;
      if (check) {
        if (!mat.has(row.first, pair.first)) {
          printf("FAILED: don't have something we inserted after inserting "
                 "elements, "
                 "row was %u, col was %u\n",
                 row.first, pair.first);
          exit(-1);
        }
        if (mat.value(row.first, pair.first) != pair.second) {
          printf("FAILED: value doesn't match after inserting elements\n");
          mat.print_arrays();
          exit(-1);
        }
      }
    }
  }
  uint64_t insert_duration = get_usecs() - start;
  printf("sum is %lu\n", mat.touch_all_sum());
  if (mat.touch_all_sum() != correct_sum) {
    printf("FAILED: sum didn't match after inserting elements\n");
    mat.print_arrays();
    exit(-1);
  }

  start = get_usecs();
  for (auto &row : correct_matrix) {
    for (auto &pair : row.second) {
      mat.remove(row.first, pair.first);
      if (check) {
        if (mat.has(row.first, pair.first)) {
          printf("have something we removed while removing elements\n");
          exit(-1);
        }
      }
    }
  }
  if (check && mat.M() != 0) {
    printf("still have elements when we shouldn't\n");
    exit(-1);
  }
  uint64_t remove_duration = get_usecs() - start;
  printf("insert duration = %lu, remove duration = %lu\n", insert_duration,
         remove_duration);
}

void matrix_values_add_remove_test(
    uint32_t el_count,
    uint32_t row_count = std::numeric_limits<uint32_t>::max(),
    bool check = false) {
  matrix_values_add_remove_test_templated<uint8_t>(el_count, row_count, check);
  matrix_values_add_remove_test_templated<uint32_t>(el_count, row_count, check);
  matrix_values_add_remove_test_templated<float>(el_count, row_count, check);
  matrix_values_add_remove_test_templated<double>(el_count, row_count, check);
}
