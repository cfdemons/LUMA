// Header file to reference mapping functions
#ifndef HEADER_OPS_MAP
#define HEADER_OPS_MAP

// Forward declarations
int idxmap (int i, int j, int k, int vel, int M, int K, int nVels);	// Function: idxmap + 2 overloads
int idxmap (int i, int j, int vel, int M, int nVels);
int idxmap (int i, int j, int M);
std::vector<int> indmapref(int coarse_i, int x_start, int coarse_j, int y_start, int coarse_k, int z_start);
																	// Function: indmapref

#endif