#include "stdafx.h"
#include "LBM_definitions.h"
#include "LBM_globalvars.h"
// File containing global variable initialisations

using namespace std; // Standard namespace in use

// Global initialisation as follows

// Lattice velocities
#if (dims == 3)

	// D3Q19 (defined as in Mawson 2013 thesis but with last column as the rest particle)
	const int c[3][nVels] =
		{
			{1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 0},
			{0, 0, 1, -1, 0, 0, 1, -1, -1, 1, 1, -1, 1, -1, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, -1, -1, 1, 1, -1, 1, -1, 0}
		};

#else

	// D2Q9 (anticlockwise numbering with 9th paricle as the rest particle)
	const int c[3][nVels] =
		{
			{1, 1, 0, -1, -1, -1, 0, 1, 0},
			{0, 1, 1, 1, 0, -1, -1,-1, 0},
			{0, 0, 0, 0, 0, 0, 0, 0, 0}
		};

#endif


// Weights
#if (dims == 3)

	const double w[nVels] = 
		{1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0, 1.0/18.0,
		1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0,
		1.0/3.0};

#else

	const double w[nVels] = 
		{1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 4.0/9.0};
	
#endif

// Lattice sound speed
const double cs = 1.0 / sqrt(3.0);

// Array of GridData structures for each grid level
GridData Grids[Nref+1];