/* File contains all the routines necessary to perform 
mapping of positions or indices.
*/

#include "stdafx.h"
#include "globalvars.h"

using namespace std;

// ***************************************************************************************************

// Routines to flatten 2/3/4D array indices to a single index

// 4D array index
int idxmap (int i, int j, int k, int v, int j_max, int k_max, int v_max) {

	// Loop over vel then k then j then i
	int idx = v + (k*v_max) + (j*v_max*k_max) + (i*v_max*k_max*j_max);

	return idx;
}

// 3D array index
int idxmap (int i, int j, int k, int j_max, int k_max) {

	// Loop over k then j then i
	int idx = k + (j*k_max) + (i*k_max*j_max);

	return idx;
}

// 2D array index
int idxmap (int i, int j, int j_max) {

	// Loop over j then i
	int idx = j + (i*j_max);

	return idx;
}

// ***************************************************************************************************

// Routine to map the index of a coarse grid site to a corresponding fine site on the level below
vector<int> indmapref(int coarse_i, int x_start, int coarse_j, int y_start, int coarse_k, int z_start) {

	// Initialise result
	vector<int> fine_ind;
	
	// Map indices
	fine_ind.insert(fine_ind.begin(), 2*(coarse_i - x_start + 1) - 2 );
	fine_ind.insert(fine_ind.begin() + 1, 2*(coarse_j - y_start + 1) - 2 );
	fine_ind.insert(fine_ind.begin() + 2, 2*(coarse_k - z_start + 1) - 2 );

	return fine_ind;
}