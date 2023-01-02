// For DELETE_ALL_VERTICES_IN_VERTEX_SUBSET_F
#include "../Chen/Struction.h"

template <typename SM> 
class MaxDegree {
  public:
    MaxDegree(
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

      removeCounter = 0;
    }
    ~MaxDegree(){
    }

    bool FindMaxDegree();

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
bool MaxDegree<SM>::FindMaxDegree()
{
    int32_t maxV;
    int32_t maxD = 0;
    bool hasEdges = true;
    if (!remainingVertices.get_n()){
        hasEdges = false;
        return hasEdges;
    }
    for (unsigned int i = 0; i < V; i++)
    {
        printf("%d \n", i);
        printf("%d \n", remainingVertices.get_n());
        printf("%d \n", G.getDegree(i));
        if (G.getDegree(i) > maxD)
        {
            maxV = i;
            maxD = G.getDegree(i);
        }
    }
    if (maxD == 0)
        hasEdges = false;
    else
    {
        VertexSubset vertices_to_delete;
	    VertexSubset maxVSet = VertexSubset(maxV, V); // creates initial frontier
        do {
            removeCounter = 0;
            vertices_to_delete = G.edgeMap(maxVSet, DELETE_ALL_VERTICES_IN_VERTEX_SUBSET_F(edgesToRemove, &removeCounter, &b_size, G), true); // mark visited
            // Write phase
            G.remove_batch(edgesToRemove, min(removeCounter, b_size));
        } while(vertices_to_delete.non_empty());
        Solution[maxV] = 1;
	    //++minimum;
	}
    return hasEdges;
}
