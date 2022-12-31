#include <algorithm> //std::sort
#include <iostream> //std::cout
#include <string> //std::string
#include <vector> //std::vector

#include "../Chen/CrownReduction.h"

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
                    int32_t* _Solution) : 
                    G(_G),
                    V(_G.get_rows()),
                    E(_G.get_cols()),
                    remainingVertices(_remainingVertices),
                    match(_match),
                    Solution(_Solution)
    {
        color = (int *)malloc(V * sizeof(int));
        par = (int *)malloc(V * sizeof(int));

        // store the numbers of cycle
        int cyclenumber = 0;
        // call DFS to mark the cycles
	    dfs_cycle(1, 0, color, par, cyclenumber);
        // function to print the cycles
        printCycles(cyclenumber);
        int cy = MAlternatingCycles(cyclenumber);
        if (cy > -1){
            printf("Cycle %d is an M-alt cycle\n", cy);
        }

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

        // This is the new find crown method.
        int32_t *parallel_cr_result = FindCrown();
    }
    ~CrownReduction(){
        free(Cycles);
        free(Parents);
        free(Depth);
        free(Pair);

    }
    /*
    void FindCrown(){
        int v, i = 0, q, u_i, w_i;
        // choose first free vertex
        // 3. Pick a vertex v ∈V\(V(CY) ∪V(M))arbitrarily;
        for (int i = 0; i < V; ++i)
            if (match[i] == -1 && cy[i] == -1){
                v = i;
                break;
            }
        // 4. Let I0={v}, H0=N(I0), i =0;
        std::vector<el_t> h_0 = G.get_neighbors(v);
        std::vector<el_t> I_q;
        I_q.push_back(v);
        // 5. While(Hi!=∅) {
        while(h_0.size()){
            std::vector<el_t> n_u_q, n_w_q;
            // 5.1 If(there is an edge e =<ui, wi> ∈M in G[Hi])then{
            for (int i = 0; i < h_0.size(); ++i){
                if(match[i] >= 0){
                    q = i;
                    u_i = i;
                    n_u_q = G.get_neighbors(u_i);
                    for (int j = 0; i < n_u_q.size(); ++j){
                        if(match[i] == match[n_u_q[j]]){
                            w_i = n_u_q[j];
                            n_w_q = G.get_neighbors(u_i);
                            break;
                        }
                    }
                    break;
                }
            }
            // While(there are different neighbors u'q, w'q 
            // of uq, wq in Iq; respectively)
            std::vector<el_t> u_prime_q_intermediates, w_prime_q_intermediates;
            u_prime_q_intermediates = intersection(I_q, n_u_q);
            w_prime_q_intermediates = intersection(I_q, n_w_q);

            std::vector<el_t> u_prime_q_candidates, w_prime_q_candidates;
            std::vector<el_t> set_diff = findDiff(u_prime_q_intermediates, w_prime_q_intermediates)
            u_prime_q_candidates = intersection(set_diff, u_prime_q_intermediates);
            w_prime_q_candidates = intersection(set_diff, w_prime_q_intermediates);
            // These are all eligible u_prime_q and w_prime_q
            if (u_prime_q_candidates.size() && w_prime_q_candidates.size()){
                // {u_q−1=NM(u'q), w_q−1=NM(w'q), q =q−1;};
            }
            // Assume {xq} = N(uq) ∩ Iq = N(wq) ∩ Iq,
            // then cy =<xq, uq, NM(uq), ..., NM(ui−1),
            // ui, wi, NM(wi−1), ..., NM(wq), wq, xq> 
            // is an M-alternating odd cycle;
            std::vector<el_t> x_prime_q = intersection(set_diff, n_u_q);
        }
    }
    */

    private:
        int32_t * FindCrown();
        SM &G;
        VertexSubset & remainingVertices;
        int32_t* Cycles;
        int32_t* Parents;
        int32_t* Depth;
        int32_t* Pair;
        int32_t* Solution;

        int V,E, v_0;
        int* match;

};


template <typename SM> 
int32_t * CrownReduction<SM>::FindCrown() {
// choose first free vertex
// 3. Pick a vertex v ∈V\(V(CY) ∪V(M))arbitrarily;
  v_0 = -1;
  for (int i = 0; i < V; ++i)
    // Unmatched and not in a previously identified 
    // M alternating cycle
    if (match[i] == -1 && !Cycles[i]){
        v_0 = i;
        break;
    }
  int64_t start = v_0;
  int64_t n = G.get_rows();
  int32_t i = 0;
  int32_t q;
  if (n == 0 || v_0 == -1) {
    return nullptr;
  }

  parallel_for(int64_t i = 0; i < n; i++) { Parents[i] = -1; }
  parallel_for(int64_t i = 0; i < n; i++) { Pair[i] = -1; }
  parallel_for(int64_t i = 0; i < n; i++) { Depth[i] = -1; }

  if (n == 0) {
    return Parents;
  }
  Parents[start] = start;
  Depth[start] = 0;
  VertexSubset frontier = VertexSubset(start, n); // creates initial frontier
  VertexSubset H;
  VertexSubset I;
  while (frontier.non_empty()) { // loop until frontier is empty
    H = G.edgeMap(frontier, H_F(Parents, Depth), true, 20);
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
        I = G.edgeMap(H, I_F(Parents, match, Depth), true, 20);
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
  return Solution;
}
