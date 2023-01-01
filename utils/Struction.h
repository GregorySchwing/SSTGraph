#include "../Chen/Struction.h"
#define PRINT_DEF 0

template <typename SM> 
class Struction {
  public:
    Struction(
                    SM &_G,
                    VertexSubset & _remainingVertices,
                    int32_t _b_size,
                    std::tuple<el_t, el_t> *_edgesToRemove,
                    std::tuple<el_t, el_t> *_edgesToInsert) : 
                    G(_G),
                    V(_G.get_rows()),
                    E(_G.get_cols()),
                    remainingVertices(_remainingVertices),
                    b_size(_b_size),
                    edgesToRemove(_edgesToRemove),
                    edgesToInsert(_edgesToInsert)
    {
      numberAntiEdges = (int32_t *)malloc(V * sizeof(int32_t));
      performStruction = (int32_t *)malloc(V * sizeof(int32_t));
      maxVertex = (int32_t *)malloc(V * sizeof(int32_t));
      numStructionNeighbors = (int32_t *)malloc(V * sizeof(int32_t));

      removeCounter = 0;
      insertCounter = 0;
    }
    ~Struction(){
      free(numberAntiEdges);
      free(performStruction);
      free(maxVertex);
      free(numStructionNeighbors);
    }

    bool FindStruction();

    private:
        SM &G;
        VertexSubset & remainingVertices;
        int V,E;

        int32_t *numberAntiEdges;
        int32_t *performStruction;
        int32_t *maxVertex;
        int32_t *numStructionNeighbors;

        int32_t removeCounter;
        int32_t insertCounter;

        int32_t b_size;
        std::tuple<el_t, el_t> *edgesToRemove;
        std::tuple<el_t, el_t> *edgesToInsert;
};

template <typename SM> 
bool Struction<SM>::FindStruction()
{
  removeCounter = 0;
  insertCounter = 0;
  bool structionPerformed = false;
  int64_t n = G.get_rows(); 

  parallel_for(int64_t i = 0; i < n; i++) { numberAntiEdges[i] = 0; }
  parallel_for(int64_t i = 0; i < n; i++) { performStruction[i] = 0; }
  parallel_for(int64_t i = 0; i < n; i++) { numStructionNeighbors[i] = 0; }

  parallel_for(int64_t i = 0; i < n; i++) { maxVertex[i] = 0; }

  VertexSubset struction = G.edgeMap(remainingVertices, SET_ANTI_EDGES_F(numberAntiEdges, G), true, 20);
  parallel_for(int64_t i = 0; i < n; i++) { numberAntiEdges[i] /= 2; }
  VertexSubset firstStructionSet = G.vertexMap(remainingVertices, SET_STRUCTION_F(numberAntiEdges, performStruction, G), true); // mark visited
  #ifdef NDEBUG 
  #if PRINT_DEF
  printf("Vertices\n");
  for (uint32_t j = 0; j < n; j++) {
    printf("%lu ", j);
  }
  printf("\nStruction flags\n");
  for (uint32_t j = 0; j < n; j++) {
    printf("%lu ", performStruction[j]);
  }
  printf("\n");
  #endif
  #endif
  VertexSubset structionMIS = G.edgeMap(remainingVertices, SOLVE_MIS_F(performStruction, G), true, 20);
  #ifdef NDEBUG 
  #if PRINT_DEF
  printf("MIS\n");
  printf("Vertices\n");
  for (uint32_t j = 0; j < n; j++) {
    printf("%lu ", j);
  }
  printf("\nStruction flags\n");
  for (uint32_t j = 0; j < n; j++) {
    printf("%lu ", performStruction[j]);
  }
  printf("\n");
  #endif
  #endif
  VertexSubset structionDeg = G.edgeMap(remainingVertices, SET_NUM_STRUCTION_NEIGHBORS_F(performStruction, numStructionNeighbors, G), true, 20);
  #ifdef NDEBUG 
  #if PRINT_DEF
  printf("Degree of struct\n");
  printf("Vertices\n");
  for (uint32_t j = 0; j < n; j++) {
    printf("%lu ", j);
  }
  printf("\nNum neighbors with set flags\n");
  for (uint32_t j = 0; j < n; j++) {
    printf("%lu ", numStructionNeighbors[j]);
  }
  printf("\n");
  #endif
  #endif
  VertexSubset maxDegree = G.edgeMap(remainingVertices, SET_LARGEST_VERTEX_STRUCT_F(maxVertex, numStructionNeighbors, performStruction, G), true, 20);
  #ifdef NDEBUG 
  #if PRINT_DEF
  printf("\nmaxDegree\n");
  for (uint32_t j = 0; j < n; j++) {
    printf("%lu ", maxVertex[j]);
  }
  printf("\n");
  #endif
  #endif
  VertexSubset fin = G.edgeMap(remainingVertices, RESOLVE_CONFLICTS_STRUCT_F(maxVertex, numStructionNeighbors, performStruction, G), true, 20);
  #ifdef NDEBUG 
  #if PRINT_DEF
  printf("\nResolve conflicts\n");
  for (uint32_t j = 0; j < n; j++) {
    printf("%lu ", performStruction[j]);
  }
  printf("\n");  
  #endif
  #endif
  VertexSubset structionSetAndNeighbors = G.vertexMap(remainingVertices, GET_STRUCTION_SET_AND_NEIGHBORS_F(performStruction, numStructionNeighbors, G), true); // mark visited
  //VertexSubset structionSetAndNeighborsDeleted = G.edgeMap(remainingVertices, DELETE_NEIGHBORHOOD_OF_STRUCTION_VERTEX_F(performStruction, maxVertex, numStructionNeighbors, edgesToRemove, &b_used, &b_size, G), true); // mark visited
  VertexSubset structionSetAndNeighborsDeleted = G.edgeMap(structionSetAndNeighbors, DELETE_ALL_VERTICES_IN_VERTEX_SUBSET_F(edgesToRemove, &removeCounter, &b_size, G), true); // mark visited

  VertexSubset structionSetOnly = G.vertexMap(remainingVertices, GET_STRUCTION_SET_F(performStruction, numStructionNeighbors, G), true); // mark visited

  // Assuming all the insertions fit in one batch, this should complete the struction op.
  while (structionSetOnly.non_empty()){
    structionPerformed = true;
    uint32_t v0 = structionSetOnly.pop();
    //printf("Perform struct operation on %lu\n", v0);
    //G.print_neighbors(v0);
    std::vector<el_t> v0_neighs = G.get_neighbors(v0);
    #ifdef NDEBUG 
  #if PRINT_DEF
    printf("Neighbors of %lu\n", v0);
    for (int i = 0; i < v0_neighs.size(); ++i)
      printf("%u \n", v0_neighs[i]);
    printf("\n");
    #endif
  #endif
    int usedVertexCounter = 0;
    std::map<std::tuple<el_t,el_t>,el_t> antiEdgeToNodeMap;
    for (int i = 0; i < v0_neighs.size(); ++i)
      for (int j = i+1; j < v0_neighs.size(); ++j){
        if (!G.has(v0_neighs[i], v0_neighs[j]))
          antiEdgeToNodeMap[std::tuple<el_t, el_t>{v0_neighs[i], v0_neighs[j]}] = v0_neighs[usedVertexCounter++];
      }
    for (auto& t : antiEdgeToNodeMap)
        std::cout << std::get<0>(t.first) << " " 
                  << std::get<1>(t.first) << " " 
                  << t.second << "\n";


    std::map<std::tuple<el_t,el_t>,el_t>::iterator it_i = antiEdgeToNodeMap.begin();
    std::map<std::tuple<el_t,el_t>,el_t>::iterator it_j;

    insertCounter = 0;

    while (it_i != antiEdgeToNodeMap.end())
    {
      // Condition 1 - 
      // remove the vertices {v0; v1; ... ; vp} from G and
      // introduce a new node vij for every anti-edge {vi; vj} in G
      // where 0 < i < j <= p;

      it_j = std::next(it_i, 1);
      while (it_j != antiEdgeToNodeMap.end())
      {
        std::cout << std::get<0>(it_i->first) << " " 
          << std::get<1>(it_i->first) << " " 
          << it_i->second << "\n"; 
        std::cout << std::get<0>(it_j->first) << " " 
          << std::get<1>(it_j->first) << " " 
          << it_j->second << "\n"; 

        // Condition 2 - 
        // add an edge (vir; vjs) if i = j and
        // (vr; vs) is an edge in G;
        if ((std::get<0>(it_i->first) == std::get<0>(it_j->first))
              && G.has(std::get<1>(it_i->first), std::get<1>(it_j->first))){
          edgesToInsert[insertCounter++] = std::tuple<el_t, el_t>{it_i->second, it_j->second};
          edgesToInsert[insertCounter++] = std::tuple<el_t, el_t>{it_j->second, it_i->second};
          printf("Adding edge (%lu, %lu)-(%lu, %lu) in new neighborhood\n",std::get<0>(it_i->first), std::get<1>(it_i->first),
          std::get<0>(it_j->first), std::get<1>(it_j->first));
        }

        // Condition 3 - 
        // if i != j, add an edge (vir; vjs)
        if (std::get<0>(it_i->first) != std::get<0>(it_j->first)){
          edgesToInsert[insertCounter++] = std::tuple<el_t, el_t>{it_i->second, it_j->second};
          edgesToInsert[insertCounter++] = std::tuple<el_t, el_t>{it_j->second, it_i->second};
          printf("Adding edge (%lu, %lu)-(%lu, %lu) in new neighborhood\n",std::get<0>(it_i->first), std::get<1>(it_i->first),
          std::get<0>(it_j->first), std::get<1>(it_j->first));
        }
        ++it_j;
      }
      // Condition 4 - 
      // for every u not in {v0; ... ; vp}, 
      // add the edge (vij ; u) if (vi; u)
      // or (vj ; u) is an edge in G.

      // exnoi
      std::cout << "Find external neighbors of i " << std::get<0>(it_i->first) << std::endl;

      std::vector<el_t> externalNeighborsOf_i = G.get_neighbors(std::get<0>(it_i->first));

      std::cout << "External neighbors of " << std::get<0>(it_i->first) << " before removal" << std::endl;
      for (auto element : externalNeighborsOf_i) {
        std::cout << element << " ";
      }
      std::cout << std::endl;

      externalNeighborsOf_i.erase( remove_if( begin(externalNeighborsOf_i),end(externalNeighborsOf_i),
          [&](auto x){return find(begin(v0_neighs),end(v0_neighs),x)!=end(v0_neighs);}), end(externalNeighborsOf_i) );

      std::cout << "External neighbors of " << std::get<0>(it_i->first) << " after N(v0) removal" << std::endl;
      for (auto element : externalNeighborsOf_i) {
        std::cout << element << " ";
      }
      std::cout << std::endl;

      externalNeighborsOf_i.erase(std::remove(externalNeighborsOf_i.begin(), externalNeighborsOf_i.end(), v0), externalNeighborsOf_i.end());

      std::cout << "External neighbors of " << std::get<0>(it_i->first) << " after v0 removal" << std::endl;
      for (auto element : externalNeighborsOf_i) {
        std::cout << element << " ";
      }
      std::cout << std::endl;

      // exnoj
      std::vector<el_t> externalNeighborsOf_j = G.get_neighbors(std::get<1>(it_i->first));
      std::cout << "Find external neighbors of j " << std::get<1>(it_i->first) << std::endl;
      std::cout << "External neighbors of " << std::get<1>(it_i->first) << " before removal" << std::endl;
      for (auto element : externalNeighborsOf_j) {
        std::cout << element << " ";
      }
      std::cout << std::endl;

      externalNeighborsOf_j.erase( remove_if( begin(externalNeighborsOf_j),end(externalNeighborsOf_j),
          [&](auto x){return find(begin(v0_neighs),end(v0_neighs),x)!=end(v0_neighs);}), end(externalNeighborsOf_j) );

      std::cout << "External neighbors of " << std::get<1>(it_i->first) << " after N(v0) removal" << std::endl;
      for (auto element : externalNeighborsOf_j) {
        std::cout << element << " ";
      }
      std::cout << std::endl;

      externalNeighborsOf_j.erase(std::remove(externalNeighborsOf_j.begin(), externalNeighborsOf_j.end(), v0), externalNeighborsOf_j.end());

      std::cout << "External neighbors of " << std::get<1>(it_i->first) << " after v0 removal" << std::endl;
      for (auto element : externalNeighborsOf_j) {
        std::cout << element << " ";
      }
      std::cout << std::endl;

      // Using default function
      std::vector<el_t> unionOfExternalNeighborsOfIAndJ;
      std::set_union(externalNeighborsOf_i.begin(), externalNeighborsOf_i.end(), externalNeighborsOf_j.begin(), externalNeighborsOf_j.end(), std::back_inserter(unionOfExternalNeighborsOfIAndJ));
      std::cout << "External neighbors of I and J are :\n";
      for (const auto &i : unionOfExternalNeighborsOfIAndJ) {
          std::cout << i << ' ';
          edgesToInsert[insertCounter++] = std::tuple<el_t, el_t>{it_i->second, i};
          edgesToInsert[insertCounter++] = std::tuple<el_t, el_t>{i, it_i->second};
      }   
      std::cout << '\n';

      ++it_i;
    }
  }

  //printf("\nBefore batch changes\n");
  //G.print_arrays();
  G.remove_batch(edgesToRemove, min(removeCounter, b_size));
  G.insert_batch(edgesToInsert, min(insertCounter, b_size));
  //printf("\nAfter batch changes\n");
  //G.print_arrays();

  return structionPerformed;

}