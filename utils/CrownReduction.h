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
bool CrownReduction<SM>::FindCrown() {
  bool vertexChanged = false;
  int64_t n = G.get_rows();
  int64_t remainingV = remainingVertices.get_n();

  int32_t i = 0;
  int32_t q;
  if (n == 0 || remainingV == 0)
    return vertexChanged;
  // choose first free vertex
  // 3. Pick a vertex v ∈V\(V(CY) ∪V(M))arbitrarily;
  v_0 = -1;
  for (int i = 0; i < V; ++i){
    // Unmatched and not in a previously identified 
    // M alternating cycle
    if (match[i] == -1 && !Cycles[i]){
        v_0 = i;
        break;
    }
  }
  int64_t start = v_0;
  if (v_0 == -1)
    return vertexChanged;

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
    H = G.edgeMap(frontier, H_F(Parents, Depth), true, 20);
    H_Set.push_back(H);
    printf("H\n");
    H.print();
    VertexSubset H_Cy = G.edgeMap(H, CYCLE_DETECTION_F(Parents, Pair, Depth), true, 20);
    // Check for cycles in H
    if (H_Cy.get_n()){
        q = i;
        printf("H_Cy\n");
        H_Cy.print();
        // Might be unneccesary to copy.
        VertexSubset H_BT_Frontier = H_Cy;
        VertexSubset xq;
        // Has to be a do while since calling get_n
        // on an uninitted xq causes a segfault
        do { // loop until a cycle converges.
            //VertexSubset H_BT_Frontier_Int = G.edgeMap(H_BT_Frontier, CYCLE_BT_F(Parents, Pair), true, 20);
            //xq = G.vertexMap(H_BT_Frontier_Int, GET_XQ_F(Pair), true); // mark visited
            VertexSubset H_BT_Frontier_Int = G.edgeMap(H_BT_Frontier, CYCLE_BT_2_F(Parents, Pair), true, 20);
            // Important to call this on the previous frontier (I_BT_Frontier).
            // Not the new frontier (I_BT_Frontier_Int) which tenatively contains xq.
            // We actually find Xq's child, then get the parent outside the loop.
            xq = G.vertexMap(H_BT_Frontier, GET_XQ_2_F(Pair, Parents), true); // mark visited
            H_BT_Frontier = H_BT_Frontier_Int;
        } while (!xq.get_n());
        // This must be Xq since all the cycles have to be the same depth.
        printf("Xq's child\n");
        xq.print();
        el_t scalar_xq = xq.pop();
        el_t cycleTail = Pair[scalar_xq];
        el_t truexq = Parents[scalar_xq];

        printf("Xq %d cycleTail %d\n", truexq, cycleTail);
        VertexSubset cycleStart = VertexSubset(cycleTail, n); // creates initial frontier
        while (cycleStart.get_n()){ // loop until a cycle converges.
            VertexSubset H_BT_Frontier_Int = G.edgeMap(cycleStart, H_SET_CYCLE_F(Parents, Pair, Cycles, truexq), true, 20);
            printf("H_BT_Frontier_Int\n");
            H_BT_Frontier_Int.print();
            cycleStart = H_BT_Frontier_Int;
        }
        printf("Cycles\n");
        for (int i = 0; i < V; ++i)
            printf("%d ", Cycles[i]);
        printf("\n");
        // xq is in an I level.
        // This handles the MM. Still need to remove CY from future calls
        while(q != 0){
            match[scalar_xq] = -1;
            match[Parents[scalar_xq]] = Parents[Parents[scalar_xq]];
            match[Parents[Parents[scalar_xq]]] = Parents[scalar_xq];
            q = q-2;
        }   
        return FindCrown();
        // {M={M\{<xq, NM(xq)>}} ∪ {<NM(xq), xq−1>},
        // where x_q−1 ∈ I_q−1 ∩ N(NM(xq)); q = q−1;}
    } else {
        ++i;
        I = G.edgeMap(H, I_F(Parents, Depth, match), true, 20);
        I_Set.push_back(I);
        printf("I\n");
        I.print();
    }

    // Check for cycles in I
    VertexSubset I_Cy = G.edgeMap(I, CYCLE_DETECTION_F(Parents, Pair, Depth), true, 20);
    if (I_Cy.get_n()){
        q = i;
        printf("I_Cy\n");
        I_Cy.print();
        // Check for cycles in I
        VertexSubset I_BT_Frontier = I_Cy;
        VertexSubset xq;
        // Has to be a do while since calling get_n
        // on an uninitted xq causes a segfault
        do { // loop until a cycle converges.
            VertexSubset I_BT_Frontier_Int = G.edgeMap(I_BT_Frontier, CYCLE_BT_2_F(Parents, Pair), true, 20);
            //VertexSubset I_BT_Frontier_Int = G.edgeMap(I_BT_Frontier, CYCLE_BT_F(Parents, Pair), true, 20);
            //xq = G.vertexMap(I_BT_Frontier_Int, GET_XQ_F(Pair), true); // mark visited
            
            // Important to call this on the previous frontier (I_BT_Frontier).
            // Not the new frontier (I_BT_Frontier_Int) which tenatively contains xq.
            // We actually find Xq's child, then get the parent outside the loop.
            xq = G.vertexMap(I_BT_Frontier, GET_XQ_2_F(Pair, Parents), true); // mark visited

            I_BT_Frontier = I_BT_Frontier_Int;
        } while (!xq.get_n());
        // This must be Xq since all the cycles have to be the same depth.
        printf("Xq's child\n");
        xq.print();
        el_t scalar_xq = xq.pop();
        el_t cycleTail = Pair[scalar_xq];
        el_t truexq = Parents[scalar_xq];

        printf("Xq %d cycleTail %d\n", truexq, cycleTail);
        VertexSubset cycleStart = VertexSubset(cycleTail, n); // creates initial frontier
        while (cycleStart.get_n()){ // loop until a cycle converges.
            VertexSubset H_BT_Frontier_Int = G.edgeMap(cycleStart, H_SET_CYCLE_F(Parents, Pair, Cycles, truexq), true, 20);
            printf("H_BT_Frontier_Int\n");
            H_BT_Frontier_Int.print();
            cycleStart = H_BT_Frontier_Int;
        }
        printf("Cycles\n");
        for (int i = 0; i < V; ++i)
            printf("%d ", Cycles[i]);
        printf("\n");
        // xq is in an H level.
        // This handles the MM. Still need to remove CY from future calls
        while(q > 1){
            match[scalar_xq] = -1;
            match[Parents[scalar_xq]] = Parents[Parents[scalar_xq]];
            match[Parents[Parents[scalar_xq]]] = Parents[scalar_xq];
            q = q-2;
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
