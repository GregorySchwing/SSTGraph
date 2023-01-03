#include <algorithm> //std::sort
#include <iostream> //std::cout
#include <string> //std::string
#include <vector> //std::vector

#include "../Chen/CrownReduction.h"
// For DELETE_ALL_VERTICES_IN_VERTEX_SUBSET_F
#include "../Chen/Struction.h"

template <typename T> 
std::vector<T> intersection(std::vector<T> v1,
                                      std::vector<T> v2){
    std::vector<T> v3;

    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());

    std::set_intersection(v1.begin(),v1.end(),
                          v2.begin(),v2.end(),
                          back_inserter(v3));
    return v3;
}

template<typename T>
std::vector<T> findDiff(std::vector<T> x, std::vector<T> y) {   // no-ref, no-const
    std::vector<T> diff;
    std::sort(x.begin(), x.end());
    std::sort(y.begin(), y.end());
    std::set_difference(x.begin(), x.end(), y.begin(), y.end(), std::back_inserter(diff));
    return diff;
}

// C++ program to print all the cycles
// in an undirected graph
#include <bits/stdc++.h>
using namespace std;

template <typename SM> 
class CrownReduction {
  public:
    CrownReduction(
                    SM &_G,
                    VertexSubset & _remainingVertices,
                    int* _match,
                    int32_t* _Solution,
                    int32_t _k = -1) : 
                    G(_G),
                    V(_G.get_rows()),
                    E(_G.get_cols()),
                    remainingVertices(_remainingVertices),
                    match(_match),
                    Solution(_Solution),
                    k(_k)
    {
        // creates Parents array, initialized to all -1, except for start
        Cycles = (int32_t *)malloc(V * sizeof(int32_t));
        // creates Parents array, initialized to all -1, except for start
        Parents = (int32_t *)malloc(V * sizeof(int32_t));
        // creates Depth array, initialized to all -1, except for start  
        // necessary for the cycle detecter, to identify edges to other
        // vertices at the same depth as me.
        Depth = (int32_t *)malloc(V * sizeof(int32_t));
        // creates Pair array, initialized to all -1
        // when an edge is shared between two vertices in H or I,
        // they set each other, so they may backtrack until they converge.
        Pair = (int32_t *)malloc(V * sizeof(int32_t));

        parallel_for(int64_t i = 0; i < V; i++) { Cycles[i] = 0; }

        // creates inCover array, initialized to all -1, except for start
        b_size = 10000; 
        b_used = 0; 
        edgesToRemove = (std::tuple<el_t, el_t> *)malloc(b_size * sizeof(std::tuple<el_t, el_t>));

        // This is the new find crown method.
        //int32_t *parallel_cr_result = FindCrown();
    }
    ~CrownReduction(){
        free(Cycles);
        free(Parents);
        free(Depth);
        free(Pair);
        free(edgesToRemove);
    }
    void ResetCycles();
    bool FindCrown();

    private:
        SM &G;
        VertexSubset & remainingVertices;
        int32_t* Cycles;
        int32_t* Parents;
        int32_t* Depth;
        int32_t* Pair;
        int32_t* Solution;

        int V,E, v_0, k;
        int* match;

        int32_t b_size; 
        int32_t b_used; 
        std::tuple<el_t, el_t> *edgesToRemove;


};


template <typename SM> 
void CrownReduction<SM>::ResetCycles() {
    parallel_for(int64_t i = 0; i < V; i++) { Cycles[i] = 0; }
}

template <typename SM> 
bool CrownReduction<SM>::FindCrown() {
  bool vertexChanged = false;
  int64_t n = G.get_rows();
  int64_t remainingV = remainingVertices.get_n();

  int32_t i = 0;
  int32_t q;
  int32_t lock = -1;
  int32_t CycleEdge_u;
  int32_t CycleEdge_v;
  if (n == 0 || remainingV == 0)
    return vertexChanged;
  // choose first free vertex
  // 3. Pick a vertex v ∈V\(V(CY) ∪V(M))arbitrarily;
  // VertexSubset eligibleStartVertices = G.vertexMap(remainingVertices, GET_UNMATCHED_F(match, Cycles, G), true); // mark visited
  VertexSubset eligibleStartVertices = G.edgeMap(remainingVertices, GET_START_F(match, Cycles, G), true, 20);
  eligibleStartVertices.print();
  if (!eligibleStartVertices.get_n())
    return vertexChanged;

  int64_t start = eligibleStartVertices.pop();

  parallel_for(int64_t i = 0; i < n; i++) { Parents[i] = -1; }
  parallel_for(int64_t i = 0; i < n; i++) { Pair[i] = -1; }
  parallel_for(int64_t i = 0; i < n; i++) { Depth[i] = -1; }

  Parents[start] = start;
  Depth[start] = 0;
  VertexSubset frontier = VertexSubset(start, n); // creates initial frontier
  std::vector<VertexSubset> H_Set;
  std::vector<VertexSubset> I_Set;
  printf("Starting at %d\n", start);
  VertexSubset H;
  VertexSubset I;
  while (frontier.non_empty()) { // loop until frontier is empty
    H = G.edgeMap(frontier, H_F(Parents, Cycles, Depth), true, 20);
    H_Set.push_back(H);
    printf("H\n");
    H.print();
    lock = -1;
    CycleEdge_u = -1;
    CycleEdge_v = -1;
    VertexSubset H_Cy = G.edgeMap(H, CYCLE_DETECTION_F(Parents, &lock, Depth, &CycleEdge_u, &CycleEdge_v), true, 20);
    // Check for cycles in H
    if (H_Cy.get_n()){
        q = i;
        printf("H_Cy\n");
        H_Cy.print();
        printf("cycle edge %u %u\n", CycleEdge_u, CycleEdge_v);
        Cycles[CycleEdge_u] = 1;
        Cycles[CycleEdge_v] = 1;
        while (Parents[CycleEdge_u] != Parents[CycleEdge_v]){
            printf("%u -> %u\n", CycleEdge_u, Parents[CycleEdge_u]);
            printf("%u -> %u\n", CycleEdge_v, Parents[CycleEdge_v]);

            CycleEdge_u = Parents[CycleEdge_u];
            CycleEdge_v = Parents[CycleEdge_v];
            Cycles[CycleEdge_u] = 1;
            Cycles[CycleEdge_v] = 1;
        }         
        el_t xq = Parents[CycleEdge_u];  
        Cycles[xq] = 1; 
        printf("Xq %u\n", xq);
        printf("Cycles\n");
        for (int i = 0; i < V; ++i)
            printf("%d ", Cycles[i]);
        printf("\n");
        // xq is in an I level.
        // This handles the MM. Still need to remove CY from future calls
        // The point is just to flip the parity of the edges 
        // that are in the matching from xq,child[xq] to start
        el_t counter = 0;            
        while(xq != start){
            // Only on first time do I just turn myself off.
            if (!counter){
                match[xq] = -1;
                counter++;
                xq = Parents[xq];
            // Only on odd counters do I match with parent.
            } else {
                if (counter % 2 == 1){
                    el_t parent = Parents[xq];  
                    match[xq] = parent;
                    match[parent] = xq;
                    xq = parent;   
                    counter++; 
                } else {
                    counter++; 
                    el_t parent = Parents[xq]; 
                    xq = parent;   
                }
            }
        }   
        return FindCrown();
        // {M={M\{<xq, NM(xq)>}} ∪ {<NM(xq), xq−1>},
        // where x_q−1 ∈ I_q−1 ∩ N(NM(xq)); q = q−1;}
    } else {
        ++i;
        I = G.edgeMap(H, I_F(Parents, Cycles, Depth, match), true, 20);
        I_Set.push_back(I);
        printf("I\n");
        I.print();
    }

    // Check for cycles in I
    lock = -1;
    CycleEdge_u = -1;
    CycleEdge_v = -1;
    VertexSubset I_Cy = G.edgeMap(I, CYCLE_DETECTION_F(Parents, &lock, Depth, &CycleEdge_u, &CycleEdge_v), true, 20);
    if (I_Cy.get_n()){
        q = i;
        printf("I_Cy\n");
        I_Cy.print();
        printf("cycle edge %u %u\n", CycleEdge_u, CycleEdge_v);
        Cycles[CycleEdge_u] = 1;
        Cycles[CycleEdge_v] = 1;
        while (Parents[CycleEdge_u] != Parents[CycleEdge_v]){
            printf("%u -> %u\n", CycleEdge_u, Parents[CycleEdge_u]);
            printf("%u -> %u\n", CycleEdge_v, Parents[CycleEdge_v]);

            CycleEdge_u = Parents[CycleEdge_u];
            CycleEdge_v = Parents[CycleEdge_v];
            Cycles[CycleEdge_u] = 1;
            Cycles[CycleEdge_v] = 1;
        }         
        el_t xq = Parents[CycleEdge_u];  
        Cycles[xq] = 1; 
        printf("Xq %u\n", xq);
        printf("Cycles\n");
        for (int i = 0; i < V; ++i)
            printf("%d ", i);
        printf("\n");
        for (int i = 0; i < V; ++i)
            printf("%d ", Cycles[i]);
        printf("\n");
        // xq is in an I level.
        // This handles the MM. Still need to remove CY from future calls
        // The point is just to flip the parity of the edges 
        // that are in the matching from xq,child[xq] to start
        el_t counter = 0;            
        while(xq != start){
            // Only on first time do I just turn myself off.
            if (!counter){
                match[xq] = -1;
                counter++;
                xq = Parents[xq];
            // Only on odd counters do I match with parent.
            } else {
                if (counter % 2 == 1){
                    el_t parent = Parents[xq];  
                    match[xq] = parent;
                    match[parent] = xq;
                    xq = parent;   
                    counter++; 
                } else {
                    counter++; 
                    el_t parent = Parents[xq]; 
                    xq = parent;   
                }
            }
        }   
        // {M={M\{<xq, NM(xq)>}} ∪ {<NM(xq), xq−1>},
        // where x_q−1 ∈ I_q−2 ∩ N(NM(xq)); q = q−1;}
        return FindCrown();
    } else {
        I_Set.push_back(I);
    }
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
  // Add H to the solution
  for (auto H_i : H_Set){
    k -= H_i.get_n();
    G.vertexMap(H_i, SET_SOLUTION_H_F(Solution), false); // mark visited
  }
  // Remove H and I from Graph.
  VertexSubset vertices_to_delete;
  for (auto H_i : H_Set){
    do {
        b_used = 0;
        vertices_to_delete = G.edgeMap(H_i, DELETE_ALL_VERTICES_IN_VERTEX_SUBSET_F(edgesToRemove, &b_used, &b_size, G), true); // mark visited
        // Write phase
        G.remove_batch(edgesToRemove, min(b_used, b_size));
    } while(vertices_to_delete.non_empty());
  }
  for (auto I_i : I_Set){
    do {
        b_used = 0;
        vertices_to_delete = G.edgeMap(I_i, DELETE_ALL_VERTICES_IN_VERTEX_SUBSET_F(edgesToRemove, &b_used, &b_size, G), true); // mark visited
        // Write phase
        G.remove_batch(edgesToRemove, min(b_used, b_size));
    } while(vertices_to_delete.non_empty());
  }
  vertexChanged = true;
  // Return Find-CROWN(G\(I∪H), M\(I∪H), CY, k −|H|), where
  // where HSet and ISet form a crown.
  VertexSubset newRemainingVertices = G.vertexMap(remainingVertices, Update_Remaining_V_F(G), true); // mark visited
  remainingVertices = newRemainingVertices;
  return vertexChanged;
}
