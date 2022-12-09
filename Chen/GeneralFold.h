
#define LEFTROTATE(a, b) (((a) << (b)) | ((a) >> (32 - (b))))

template <typename T, typename SM> 
struct SELECT_COLOR_F {
    T *match;
    const SM &G;
	uint random;
    bool keepMatchingTBB;
    const uint selectBarrier = 10;

    //Nothing-up-my-sleeve working constants from SHA-256.
    const uint32_t MD5K[64] = {0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
                    0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
                    0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
                    0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
                    0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
                    0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
                    0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
                    0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};

    //Rotations from MD5.
    const uint32_t MD5R[64] = {7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22,
                    5,  9, 14, 20, 5,  9, 14, 20, 5,  9, 14, 20, 5,  9, 14, 20,
                    4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23,
                    6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21};
                
    explicit SELECT_COLOR_F(T *_match, const SM &_G) : match(_match), G(_G) {}
    inline bool operator()(uintE i) {
        //This code should be the same as in matchgpu.cu!
        if (match[i] >= 2) return false;
        
        //There are still vertices to be matched.
        keepMatchingTBB = true;
        
        //Start hashing.
        uint h0 = 0x67452301, h1 = 0xefcdab89, h2 = 0x98badcfe, h3 = 0x10325476;
        uint a = h0, b = h1, c = h2, d = h3, e, f, g = i;

        for (int j = 0; j < 16; ++j)
        {
            f = (b & c) | ((~b) & d);

            e = d;
            d = c;
            c = b;
            b += LEFTROTATE(a + f + MD5K[j] + g, MD5R[j]);
            a = e;

            h0 += a;
            h1 += b;
            h2 += c;
            h3 += d;

            g *= random;
        }
        
        match[i] = ((h0 + h1 + h2 + h3) < selectBarrier ? 0 : 1);
        return true;
    }
};


template <typename T, typename SM> struct REQUEST_F {
  T *match;
  T *requests;
  int64_t nrVertices;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  REQUEST_F(T* _match, T* _requests, SM &_G) : 
  match(_match),
  requests(_requests),
  G(_G)  {
    nrVertices = G.get_rows(); 
  }
  // Assumes undirected graph
  inline bool update(uint32_t s, uint32_t d) { // Update
    //Look at all blue vertices and let them make requests.
    if (match[d] == 0)
    {
        //const int2 indices = neighbourRanges[i];
        int candidate = nrVertices;
        int dead = 1;

        //for (int j = indices.x; j < indices.y; ++j)
        //{
            //Only propose to red neighbours.
            //const int ni = neighbours[j];
            const int nm = match[s];
            //Do we have an unmatched neighbour?
            if (nm < 4)
            {
                //Is this neighbour red?
                if (nm == 1)
                {
                    //Propose to this neighbour.
                    requests[d] = s;
                    dead = 2;
                    //break;
                }
                else
                {
                    dead = 0;
                }
            }
        //}

        requests[d] = candidate + dead;
    }
    else
    {
        //Clear request value.
        requests[d] = nrVertices;
    }

    return true;
  }
  // ???
  // Assumes undirected graph
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    // Assuming I am not a neighbor to myself,
    // thus G.common_neighbors(s,d) > 0 indicates a triangle.
    uint32_t edgeIndex = __sync_fetch_and_add(&requests[d], (G.getDegree(d) - G.common_neighbors(s,d) - 1));
    return true;
  }
  
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};



template <typename T, typename SM> struct RESPOND_F {
  T *match;
  T *requests;
  int64_t nrVertices;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  RESPOND_F(T* _match, T* _requests, SM &_G) : 
  match(_match),
  requests(_requests),
  G(_G)  {
    nrVertices = G.get_rows(); 
  }
  // Assumes undirected graph
  inline bool update(uint32_t s, uint32_t d) { // Update
    //Look at all red vertices.
    if (match[d] == 1)
    {
        //const int2 indices = neighbourRanges[i];

        //Select first available proposer.
        //for (int j = indices.x; j < indices.y; ++j)
        //{
            //const int ni = neighbours[j];

            //Only respond to blue neighbours.
            if (match[s] == 0)
            {
                //Avoid data thrashing be only looking at the request value of blue neighbours.
                if (requests[s] == d)
                {
                    requests[d] = s;
                    //break;
                }
            }
        //}
    }

    return true;
  }
  // ???
  // Assumes undirected graph
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    // Assuming I am not a neighbor to myself,
    // thus G.common_neighbors(s,d) > 0 indicates a triangle.
    uint32_t edgeIndex = __sync_fetch_and_add(&requests[d], (G.getDegree(d) - G.common_neighbors(s,d) - 1));
    return true;
  }
  
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};


template <typename T, typename SM> 
struct MATCH_F {
    T *match;
    T *requests;
    int64_t nrVertices;
    const SM &G;
    explicit MATCH_F(T* _match, T* _requests, SM &_G) : 
    match(_match),
    requests(_requests),
    G(_G)  {
        nrVertices = G.get_rows(); 
    }
    inline bool operator()(uintE i) {
        uintE r = requests[i];

        //Only unmatched vertices make requests.
        if (r == nrVertices + 1)
        {
            //This is vertex without any available neighbours, discard it.
            match[i] = 2;
        }
        else if (r < nrVertices)
        {
            //This vertex has made a valid request.
            if (requests[r] == i)
            {
                //Match the vertices if the request was mutual.
                match[i] = 4 + min(i, r);
            }
        }
        return true;
    }
};