
#pragma once
#include "../SparseMatrix.hpp"
#include "Match.h"
#include <limits>


#include <bits/stdc++.h>
using namespace std;

template <typename T, typename SM> 
struct GET_LEVEL_I_F {
  int64_t i;
  T *evenlevel;
  T *oddlevel;
  const SM &G;

  explicit GET_LEVEL_I_F(int64_t _i, T *_evenlevel, T *_oddlevel, const SM &_G) : i(_i), evenlevel(_evenlevel), oddlevel(_oddlevel), G(_G) {}
  inline bool operator()(uintE v) {
    return min(evenlevel[v], oddlevel[v]) == i;
  }
};



template <typename T, typename SM> struct EVEN_LEVEL_F {
  T * evenlevel;
  T * oddlevel;
  T * bridges;
  T * predecessors;
  T * anomalies;
  T * bridgesCount;
  T * predecessorsCount;
  T * anomaliesCount;
  bool *unvisited_v;
  bool *unused_e;
  bool *unvisited_e;
  int64_t i;
  SM &predecessors_SST;
  SM &anomalies_SST;
  SM &bridges_SST;
  const SM &G;
  EVEN_LEVEL_F(T * _evenlevel,
        T * _oddlevel,
        T * _bridges,
        T * _predecessors,
        T * _anomalies,
        T * _bridgesCount,
        T * _predecessorsCount,
        T * _anomaliesCount,
        bool *_unvisited_v,
        bool *_unused_e,
        bool *_unvisited_e,
        int64_t _i,
        SM &_predecessors_SST,
        SM &_anomalies_SST,
        SM &_bridges_SST,
        const SM &_G) :
  evenlevel(_evenlevel),
  oddlevel(_oddlevel),
  bridges(_bridges),
  predecessors(_predecessors),
  anomalies(_anomalies),
  unvisited_v(_unvisited_v),
  unused_e(_unused_e),
  unvisited_e(_unvisited_e),
  i(_i),
  predecessors_SST(_predecessors_SST),
  anomalies_SST(_anomalies_SST),
  bridges_SST(_bridges_SST),
  G(_G)  {}

  /*
    In dense mode, EdgeMap loops over all vertices and
    looks at incoming edges to see if the source is part of the
    vertex set. This does not require locking because each vertex
    only updates itself and is preferred when the vertex set is large.

    I am the destination.  I only update inCover[d].
    Degrees are const in this class.

  */

  inline bool update(uint32_t s, uint32_t d) { // Update
    // for each v with evenlevel == i
    if (evenlevel[d] == i){
      // unmatched, "unused" neighbors -- ??? Should "unused" be "unvisited"
      if (evenlevel[s] < 4){

      //if (match[s] < 4){
        if (evenlevel[s] < std::numeric_limits<int32_t>::max()){
          int64_t temp = (evenlevel[s] + evenlevel[d])/2;
          //bridges[temp] = ???
        } else {
          if (oddlevel[s] == std::numeric_limits<int32_t>::max()){
            oddlevel[s] = i + 1;
          }
          if (oddlevel[s] == i + 1){
            //predecessors[s] = ???
          }
          if (oddlevel[s] < i){
            //anomalies[s] = ???
          }
        } 
        return false;
      }
    }
  }
  /*
    In sparse mode, EdgeMap
    iterates over the outgoing edges of each vertex in the subset
    and updates the destination vertex for each edge. Because it_i is
    run in parallel, synchronization must be used when accessing
    the destination vertex data.

    I am the source.  I update inCover[d] using synchronization.
    Degrees are const in this class.

  */
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    //if (performStruction[s]){
    //  __sync_fetch_and_and(&numStructionNeighbors[d], performStruction[s]);
    //}
    return false;
  }
  // cond function checks if vertex in remaining vertices set has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
    //return G.getDegree(d) > 0;
  }
};


template <typename SM> 
class MaximumMatcherMicali {
  public:
    MaximumMatcherMicali(
                    SM &_G,
                    VertexSubset & _remainingVertices) : 
                    G(_G),
                    n(_G.get_rows()),
                    e(_G.get_cols()),
                    remainingVertices(_remainingVertices),
                    predecessors_SST(SparseMatrixV<true, bool>(_G.get_rows(), _G.get_rows())),
                    anomalies_SST(SparseMatrixV<true, bool>(_G.get_rows(), _G.get_rows())),
                    bridges_SST(SparseMatrixV<true, bool>(_G.get_rows(), _G.get_rows()))
{
  evenlevel = (int32_t *)malloc(n * sizeof(int32_t));
  oddlevel = (int32_t *)malloc(n * sizeof(int32_t));
  blossom = (int32_t *)malloc(n * sizeof(int32_t));
  predecessors = (int32_t *)malloc(n * sizeof(int32_t));
  anomalies = (int32_t *)malloc(n * sizeof(int32_t));
  predecessorsCount = (int32_t *)malloc(n * sizeof(int32_t));
  anomaliesCount = (int32_t *)malloc(n * sizeof(int32_t));

  unvisited_v = (bool *)malloc(n * sizeof(bool));

  unused_e = (bool *)malloc(e * sizeof(bool));
  unvisited_e = (bool *)malloc(e * sizeof(bool));

  bridges = (int32_t *)malloc(n * sizeof(int32_t));
  bridgesCount = (int32_t *)malloc(n * sizeof(int32_t));

  request = (int32_t *)malloc(n * sizeof(int32_t));
  match = (int32_t *)malloc(n * sizeof(int32_t));

  m.Match(G, remainingVertices, request, match);

  // 0
  for(int64_t i = 0; i < n; i++) { 
    evenlevel[i] = std::numeric_limits<int32_t>::max();
    oddlevel[i] = std::numeric_limits<int32_t>::max();
    blossom[i] = -1;
    predecessors[i] = -1;
    anomalies[i] = -1;
    unvisited_v[i] = true;
  }
  for(int64_t i = 0; i < e; i++) { 
    unused_e[i] = true;
    unvisited_e[i] = true;
  }
  i = -1;

}

    int Search() {
      // 1
      for(int64_t i = 0; i < n; i++)
        if (match[i] < 4) evenlevel[i] = 0;

      // 2 
      VertexSubset level_i_vertices = G.vertexMap(remainingVertices, GET_LEVEL_I_F(i, evenlevel, oddlevel, G), true); // mark visited
      if (!level_i_vertices.non_empty())
        return 0;

      // 3 
      if (i % 2 == 0){
        VertexSubset vertices_to_delete = G.edgeMap(remainingVertices, EVEN_LEVEL_F(evenlevel,
        oddlevel,
        bridges,
        predecessors,
        anomalies,
        bridgesCount,
        predecessorsCount,
        anomaliesCount,
        unvisited_v,
        unused_e,
        unvisited_e,
        i,
        predecessors_SST,
        anomalies_SST,
        bridges_SST,
        G), true, 20);
      }
    }

    private:
        SM &G;
        VertexSubset & remainingVertices;
        Matcher m;

        // |V| = n
        int64_t n;
        int64_t e;
        int64_t i;

        // Dynamic length
        int32_t *evenlevel;
        int32_t *oddlevel;
        int32_t *blossom;

        // Dynamic length
        SM predecessors_SST;
        SM anomalies_SST;
        SM bridges_SST;
        // Dynamic length

        int32_t *predecessors;
        int32_t *anomalies;
        int32_t *bridges;


        // Default count size == 5
        int32_t *predecessorsCount;
        int32_t *anomaliesCount;
        int32_t *bridgesCount;

        bool *unvisited_v;

        // Length E
        bool *unused_e;
        bool *unvisited_e;


        // Matcher variables
        int32_t *request;
        int32_t *match;

};

struct struct_edge{int v;struct_edge* n;};
typedef struct_edge* edge;

template <typename SM> 
class MaximumMatcherBlossom {
  public:
    MaximumMatcherBlossom(
                    SM &_G,
                    VertexSubset & _remainingVertices) : 
                    G(_G),
                    V(_G.get_rows()),
                    E(_G.get_cols()),
                    remainingVertices(_remainingVertices)
    {
      pool = (struct_edge *)malloc(V * V * sizeof(struct_edge));
      match = (int *)malloc(V * sizeof(int));
      q = (int *)malloc(V * sizeof(int));
      father = (int *)malloc(V * sizeof(int));
      base = (int *)malloc(V * sizeof(int));

      inq = (bool *)malloc(V * sizeof(bool));
      inb = (bool *)malloc(V * sizeof(bool));

      //ed = (bool**)malloc(V * sizeof(bool*));
      //for (int i = 0; i < V; i++)
      //  ed[i] = (bool*)malloc(V * sizeof(bool));

      
    }

    int edmonds();

    private:
        SM &G;
        VertexSubset & remainingVertices;

        struct_edge *pool;
        edge top=pool;
        int V,E,*match,qh,qt,*q,*father,*base;
        bool *inq,*inb;
        //struct_edge pool[M*M*2];
        //edge top=pool,adj[M];
        //int V,E,match[M],qh,qt,q[M],father[M],base[M];
        //bool inq[M],inb[M],ed[M][M];

        void add_edge(int u,int v);
        int LCA(int root,int u,int v);
        void mark_blossom(int lca,int u);
        void blossom_contraction(int s,int u,int v);
        int find_augmenting_path(int s);
        int augment_path(int s,int t);

};
template <typename SM> 
int MaximumMatcherBlossom<SM>::LCA(int root,int u,int v)
{
  static bool *inp;
  inp = (bool *)malloc(V * sizeof(bool)); 
  memset(inp,0,V*sizeof(bool));
  while(1)
    {
      inp[u=base[u]]=true;
      if (u==root) break;
      u=father[match[u]];
    }
  while(1)
    {
      if (inp[v=base[v]]) return v;
      else v=father[match[v]];
    }
    free(inp);
}
template <typename SM> 
void MaximumMatcherBlossom<SM>::mark_blossom(int lca,int u)
{
  while (base[u]!=lca)
    {
      int v=match[u];
      inb[base[u]]=inb[base[v]]=true;
      u=father[v];
      if (base[u]!=lca) father[u]=v;
    }
}
template <typename SM> 
void MaximumMatcherBlossom<SM>::blossom_contraction(int s,int u,int v)
{
  int lca=LCA(s,u,v);
  memset(inb,0,V*sizeof(bool));
  mark_blossom(lca,u);
  mark_blossom(lca,v);
  if (base[u]!=lca)
    father[u]=v;
  if (base[v]!=lca)
    father[v]=u;
  for (int u=0;u<V;u++)
    if (inb[base[u]])
      {
	base[u]=lca;
	if (!inq[u])
	  inq[q[++qt]=u]=true;
      }
}
template <typename SM> 
int MaximumMatcherBlossom<SM>::find_augmenting_path(int s)
{
  memset(inq,0,V*sizeof(bool));
  memset(father,-1,V*sizeof(int));
  for (int i=0;i<V;i++) base[i]=i;
  inq[q[qh=qt=0]=s]=true;
  while (qh<=qt)
    {
      int u=q[qh++];
      // adj list is 2d array of edges.
      // edges store v and a ptr to the next e in the array.
      // this avoids a counter for each row.
      // the final edge has a null ptr hence
      // e is false and the loop terminates.
      // therefore, this should iterate over all the 
      // edges belonging to u.
      // it's possible this could be parallelized in the future.

      // adj list
      // for (edge e=adj[u];e;e=e->n)
      std::vector<el_t> u_neighbors = G.get_neighbors(u);
      for(int v : u_neighbors) 
        {
          //int v=e->v;
          if (base[u]!=base[v]&&match[u]!=v)
            if ((v==s)||(match[v]!=-1 && father[match[v]]!=-1))
              blossom_contraction(s,u,v);
            else if (father[v]==-1)
            {
              father[v]=u;
              if (match[v]==-1)
                return v;
              else if (!inq[match[v]])
                inq[q[++qt]=match[v]]=true;
            }
        }
    }
  return -1;
}
template <typename SM> 
int MaximumMatcherBlossom<SM>::augment_path(int s,int t)
{
  int u=t,v,w;
  while (u!=-1)
    {
      v=father[u];
      w=match[v];
      match[v]=u;
      match[u]=v;
      u=w;
    }
  return t!=-1;
}
template <typename SM> 
int MaximumMatcherBlossom<SM>::edmonds()
{
  int matchc=0;
  memset(match,-1,V*sizeof(int));
  for (int u=0;u<V;u++)
    if (match[u]==-1)
      matchc+=augment_path(u,find_augmenting_path(u));
  for (int i=0;i<V;i++)
    if (i<match[i])
      cout<<i<<" "<<match[i]<<endl;
  return matchc;
}
