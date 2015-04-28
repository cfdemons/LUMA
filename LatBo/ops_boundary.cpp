/* This file holds all the code for applying standard LBM boundary conditions such as bounce-back, inlet and outlet. */

#include "stdafx.h"
#include "LBM_definitions.h"
#include "LBM_globalvars.h"

using namespace std;

// ***************************************************************************************************

// Boundary condition application routine
// Supply an integer r indicating from which level the algorithm is to be executed plus a flag to specify
// which type of condition should be applied.
void LBM_boundary (int r, int bc_type_flag) {

	// Flag == 0 means apply all conditions in file
	if (bc_type_flag == 0) {

		// Apply all

	}


}

// ***************************************************************************************************