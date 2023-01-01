#include "../Chen/Dominated.h"
#include "../Chen/Struction.h"


template <typename SM> 
class Dominated {
  public:
    Dominated(
                    SM &_G,
                    VertexSubset & _remainingVertices,
                    int32_t _b_size,
                    std::tuple<el_t, el_t> *_edgesToRemove,
                    int32_t* _Solution
                    ) : 
                    G(_G),
                    V(_G.get_rows()),
                    E(_G.get_cols()),
                    remainingVertices(_remainingVertices),
                    b_size(_b_size),
                    edgesToRemove(_edgesToRemove),
                    Solution(_Solution)
    {
      vertexDominates = (int32_t *)malloc(V * sizeof(int32_t));

      removeCounter = 0;
    }
    ~Dominated(){
      free(vertexDominates);
    }

    bool FindDominated();

    private:
        SM &G;
        VertexSubset & remainingVertices;
        int V,E;

        int32_t *vertexDominates;

        int32_t removeCounter;

        int32_t b_size;
        std::tuple<el_t, el_t> *edgesToRemove;
        int32_t* Solution;
};

template <typename SM> 
bool Dominated<SM>::FindDominated()
{
  int64_t n = G.get_rows(); 
  parallel_for(int64_t i = 0; i < n; i++) { vertexDominates[i] = 0; }
  VertexSubset dominates = G.edgeMap(remainingVertices, SET_DOMINATED_F(vertexDominates, G), true, 20);
  bool vertexChanged = false;
  VertexSubset vertices_to_delete;
  while (dominates.non_empty()) { // loop until no dominates remain
    printf("DOMINATING VERTICES\n");
    dominates.print();
    vertexChanged = true;
    removeCounter = 0;
    //__sync_fetch_and_and(&b_used, 0);
    // returns vertices to delete
    vertices_to_delete = G.edgeMap(dominates, DELETE_ALL_VERTICES_IN_VERTEX_SUBSET_F(edgesToRemove, &removeCounter, &b_size, G), true); // mark visited

    // Write phase
    G.remove_batch(edgesToRemove, min(removeCounter, b_size));
    dominates = G.edgeMap(remainingVertices, SET_DOMINATED_F(vertexDominates, G), true, 20);
  }
  return vertexChanged;
}
