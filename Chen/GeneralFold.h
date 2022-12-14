#define NR_MATCH_ROUNDS 20
#define NR_MAX_MATCH_ROUNDS 256
#define LEFTROTATE(a, b) (((a) << (b)) | ((a) >> (32 - (b))))

    
template <typename T, typename SM> struct REQUEST_2_F {
  T *match;
  T *request;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  REQUEST_2_F(T *_numToEdgesRemove, T *_maxNumToEdgesRemove, SM &_G) : 
  match(_numToEdgesRemove),
  request(_maxNumToEdgesRemove),
  G(_G)  {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    //printf("s %u d %u\n",s, d);
    bool deletedEdge = false;
    //printf("s %u d %u\n", s, d);
    //printf("s %u match %u d %u match %u\n", s, match[s], d, match[d]);
    //Look at all blue vertices and let them make requests.
    // match[d] == 0 : blue
    // match[s] == 1 : red
    if (match[d] == 0 && match[s] == 1){
      // a match
      request[d] = s;
    }
    return deletedEdge;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    bool deletedEdge = false;
    //printf("s %u d %u\n",s, d);
    // If source is red and dest is blue, 
    if (match[d] == 0 && match[s] == 1){
      // Safe for multiple vertices trying to match with d.  
      // only the first will see match[d] == 0, and swap it's value with s. 
      //deletedEdge = CAS(&request[d], (T)0, (T)s);
      // not sure this is neccessary since a race condition of who wins the 
      // matching isn't neccessary a problem.
      request[d] = s;
    }
    return deletedEdge;
  }
  
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};


template <typename T, typename SM> struct RESPOND_2_F {
  T *match;
  T *request;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  RESPOND_2_F(T *_numToEdgesRemove, T *_maxNumToEdgesRemove, SM &_G) : 
  match(_numToEdgesRemove),
  request(_maxNumToEdgesRemove),
  G(_G)  {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    bool deletedEdge = false;
    //Look at all red vertices.
    //Only respond to blue neighbours.
    if (match[d] == 1 && match[s] == 0)
    {
      //Avoid data thrashing be only looking at the request value of blue neighbours.
      if (request[s] == d)
      {
          request[d] = s;
          //break;
      }
    }
    return deletedEdge;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    bool deletedEdge = false;    
    //Look at all red vertices.
    //Only respond to blue neighbours.
    if (match[d] == 1 && match[s] == 0)
    {
      //Avoid data thrashing be only looking at the request value of blue neighbours.
      if (request[s] == d)
      {
          request[d] = s;
          //break;
      }
    }
    return deletedEdge;
  }
  
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};

template <typename T, typename SM> 
struct SELECT_COLOR_F {
    T *match;
    const SM &G;
    const uint random;
    bool keepMatchingTBB;
    //const uint selectBarrier = 7;
    uint selectBarrier = 0x8000000;

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
                
    explicit SELECT_COLOR_F(T *_match, const SM &_G, uint _random = 5, uint barrier = 8) : match(_match), G(_G), random(_random) {
			const int b = max((uint)0, barrier);
      if (b >= 16) selectBarrier = 0xffffffff;
      else selectBarrier = (unsigned int)((long)b*0x10000000L);
		}
    
    inline bool operator()(uintE i) {
        //printf("called select vertex %u\n", i);

        //This code should be the same as in matchgpu.cu!
        if (match[i] >= 2 || !G.getDegree(i)) return false;
        
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
        //printf("vertex %u val %lu selBar %u\n", i, (h0 + h1 + h2 + h3), selectBarrier);
        match[i] = ((h0 + h1 + h2 + h3) < selectBarrier ? 0 : 1);
        return true;
    }
};


template <typename T, typename SM> 
struct SELECT_COLOR_AUX_F {
    T *match;
    T *auxMatch;

    const SM &G;

    explicit SELECT_COLOR_AUX_F(T *_match, T *_auxMatch, const SM &_G) : match(_match), auxMatch(_auxMatch), G(_G){}
    
    inline bool operator()(uintE i) {
        //printf("called select vertex %u\n", i);

        //This code should be the same as in matchgpu.cu!
        if (!G.getDegree(i)) return false;
        
        //printf("vertex %u val %lu selBar %u\n", i, (h0 + h1 + h2 + h3), selectBarrier);
        // This could be done inplace, unless the original M1 matching is needed in later steps.
        auxMatch[i] = match[i] >= 4 ? 0 : 1;
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
        bool unmatched = true;
        uintE r = requests[i];

        //Only unmatched vertices make requests.
        if (r < nrVertices)
        {
            //This vertex has made a valid request.
            if (requests[r] == i)
            {
                //Match the vertices if the request was mutual.
                match[i] = 4 + min(i, r);
                unmatched = false;
            }
        }
        return unmatched;
    }
};


template <typename T, typename SM> struct SET_NEIGHBORS_OF_UNMATCHED_O_F {
  T *match;
  T *H_n;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  SET_NEIGHBORS_OF_UNMATCHED_O_F(T *_numToEdgesRemove, T *_maxNumToEdgesRemove, SM &_G) : 
  match(_numToEdgesRemove),
  H_n(_maxNumToEdgesRemove),
  G(_G)  {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    bool deletedEdge = false;
    if (match[s] == 1)
    {
      H_n[d] = 1;
    }
    return deletedEdge;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    bool deletedEdge = false;    
    if (match[s] > 4)
    {
      H_n[d] = 1;
    }
    return deletedEdge;
  }
  
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};


template <typename T, typename SM> struct SET_NEIGHBORS_WITHIN_AUX_MATCHING_F {
  T *match;
  T *H_n;
  T *NM2HN;
  // vertex* V;
  // PR_F(double* _p_curr, double* _p_next, vertex* _V) :
  SM &G;
  SET_NEIGHBORS_WITHIN_AUX_MATCHING_F(T *_numToEdgesRemove, T *_maxNumToEdgesRemove, T *_maxNumToEdgesRemove2, SM &_G) : 
  match(_numToEdgesRemove),
  H_n(_maxNumToEdgesRemove),
  NM2HN(_maxNumToEdgesRemove2),
  G(_G)  {}
  inline bool update(uint32_t s, uint32_t d) { // Update
    bool deletedEdge = false;
    //if N(Hn) and Hn need to be disjoint
    //if (H_n[s] && match[d] >= 4)
    if (H_n[s] && match[d] >= 4 && !H_n[d])
    {
      NM2HN[d] = 1;
    }
    return deletedEdge;
  }
  inline bool updateAtomic(uint32_t s, uint32_t d) { // atomic version of Update
    bool deletedEdge = false;    
    //if N(Hn) and Hn need to be disjoint
    //if (H_n[s] && match[d] >= 4)
    if (H_n[s] && match[d] >= 4 && !H_n[d])
    {
      NM2HN[d] = 1;
    }
    return deletedEdge;
  }
  
  // cond function checks if vertex has non-zero degree
  inline bool cond(uint32_t d) {
    return true;
  }
};