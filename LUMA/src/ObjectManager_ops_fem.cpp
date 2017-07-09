/*
 * --------------------------------------------------------------
 *
 * ------ Lattice Boltzmann @ The University of Manchester ------
 *
 * -------------------------- L-U-M-A ---------------------------
 *
 *  Copyright (C) The University of Manchester 2017
 *  E-mail contact: info@luma.manchester.ac.uk
 *
 * This software is for academic use only and not available for
 * further distribution commericially or otherwise without written consent.
 *
 */

#include "../inc/ObjectManager.h"

// *****************************************************************************
/// \brief	Construct the R vector for all bodies on current level.
///
///
/// \param level current grid level
void ObjectManager::fem_constructRVector(int level) {

	// If using MPI then need to gather forces from off rank markers
#ifdef L_BUILD_FOR_MPI

	// Buffer containing all off rank forces
	std::vector<std::vector<std::vector<double>>> forceBuffer;
	std::vector<std::vector<int>> markerIdx;

	// Get the off-rank marker forces
	fem_getOffRankForces(level, markerIdx, forceBuffer);
#endif

	// Vector of marker forces
	std::vector<std::vector<double>> IBMforce_xyz;

	// Required parameters
	double forceScale, length, a, b;
	std::vector<double> Rlocal, RGlobal;
	std::vector<double> F(L_DIMS, 0.0);
	std::vector<std::vector<double>> T(L_DIMS, std::vector<double>(L_DIMS, 0.0));
	int IBnode;

	// Loop through all bodies that this rank owns
	for (int ib = 0; ib < IdxFEM.size(); ib++) {

		// Check if body is on this level
		if (iBody[IdxFEM[ib]]._Owner->level == level) {

			// Set R vector to zero
			fill(iBody[IdxFEM[ib]].fBody->R.begin(), iBody[IdxFEM[ib]].fBody->R.end(), 0.0);

			// Resize vector
			Rlocal.resize(iBody[IdxFEM[ib]].fBody->DOFsPerElement, 0.0);

			// Get force scaling parameter
			forceScale = iBody[IdxFEM[ib]]._Owner->dm / SQ(iBody[IdxFEM[ib]]._Owner->dt);

			// Sort IBM forces into correct format (if using MPI this has already been done
#ifndef L_BUILD_FOR_MPI

			// Resize force vector
			IBMforce_xyz.resize(iBody[IdxFEM[ib]].markers.size(), std::vector<double>(L_DIMS, 0.0));

			// Insert forces into vector (and scale to IBM force)
			for (int m = 0; m < iBody[IdxFEM[ib]].markers.size(); m++) {
				for (int d = 0; d < L_DIMS; d++)
					IBMforce_xyz[m][d] = iBody[IdxFEM[ib]].markers[m].force_xyz[d] * iBody[IdxFEM[ib]].markers[m].epsilon * 1.0 * forceScale;
			}
#else
			// Get number of markers
			int nMarkers = iBody[IdxFEM[ib]].markers.size() + markerIdx[IdxFEM[ib]].size();

			// Resize force vector
			IBMforce_xyz.resize(nMarkers, std::vector<double>(L_DIMS, 0.0));

			// First on-rank markers in
			for (int m = 0; m < iBody[IdxFEM[ib]].markers.size(); m++) {
				for (int d = 0; d < L_DIMS; d++)
					IBMforce_xyz[iBody[IdxFEM[ib]].markers[m].id][d] = iBody[IdxFEM[ib]].markers[m].force_xyz[d] * iBody[IdxFEM[ib]].markers[m].epsilon * 1.0 * forceScale;
			}

			// Now fill in the rest of them
			int markerID;
			for (int m = 0; m < markerIdx[IdxFEM[ib]].size(); m++) {

				// Get marker ID
				markerID = markerIdx[IdxFEM[ib]][m];

				// Now fill in marker forces
				for (int d = 0; d < L_DIMS; d++)
					IBMforce_xyz[markerID][d] = forceBuffer[IdxFEM[ib]][m][d] * 1.0 * forceScale;
			}
#endif

			// Loop through all FEM elements
			for (int el = 0; el < iBody[IdxFEM[ib]].fBody->elements.size(); el++) {

				// Get element parameters
				length = iBody[IdxFEM[ib]].fBody->elements[el].length;

				// Now loop through all child IB nodes this element has
				for (int node = 0; node < iBody[IdxFEM[ib]].fBody->elements[el].IBChildNodes.size(); node++) {

					// Get IB node and integration ranges
					IBnode = iBody[IdxFEM[ib]].fBody->elements[el].IBChildNodes[node].nodeID;
					a = iBody[IdxFEM[ib]].fBody->elements[el].IBChildNodes[node].zeta1;
					b = iBody[IdxFEM[ib]].fBody->elements[el].IBChildNodes[node].zeta2;

					// Get subset of transpose matrix
					T = {{iBody[IdxFEM[ib]].fBody->elements[el].T[0][0], iBody[IdxFEM[ib]].fBody->elements[el].T[0][1]},
						 {iBody[IdxFEM[ib]].fBody->elements[el].T[1][0], iBody[IdxFEM[ib]].fBody->elements[el].T[1][1]}};

					// Convert force to local coordinates
					F = GridUtils::matrix_multiply(T, IBMforce_xyz[IBnode]);

					// Get the nodal values by integrating over range of IB point
					Rlocal[0] = F[0] * 0.5 * length * (0.5 * b - 0.5 * a + 0.25 * SQ(a) - 0.25 * SQ(b));
					Rlocal[1] = F[1] * 0.5 * length * (0.5 * b - 0.5 * a - SQ(a) * SQ(a) / 16.0 + SQ(b) * SQ(b) / 16.0 + 3.0 * SQ(a) / 8.0 - 3.0 * SQ(b) / 8.0);
					Rlocal[2] = F[1] * 0.5 * length * (length * (-SQ(a) * SQ(a) + SQ(b) * SQ(b)) / 32.0 - length * (-TH(a) + TH(b)) / 24.0 - length * (-SQ(a) + SQ(b)) / 16.0 + length * (b - a) / 8.0);
					Rlocal[3] = F[0] * 0.5 * length * (-0.25 * SQ(a) + 0.25 * SQ(b) + 0.5 * b - 0.5 * a);
					Rlocal[4] = F[1] * 0.5 * length * (0.5 * b - 0.5 * a + SQ(a) * SQ(a) / 16.0 - SQ(b) * SQ(b) / 16.0 - 3.0 * SQ(a) / 8.0 + 3.0 * SQ(b) / 8.0);
					Rlocal[5] = F[1] * 0.5 * length * (length * (-SQ(a) * SQ(a) + SQ(b) * SQ(b)) / 32.0 + length * (-TH(a) + TH(b)) / 24.0 - length * (-SQ(a) + SQ(b)) / 16.0 - length * (b - a) / 8.0);

					// Get element internal forces
					RGlobal = GridUtils::matrix_multiply(GridUtils::matrix_transpose(iBody[IdxFEM[ib]].fBody->elements[el].T), Rlocal);

					// Add to global vector
					GridUtils::assembleGlobalVec(el, iBody[IdxFEM[ib]].fBody->DOFsPerNode, RGlobal, iBody[IdxFEM[ib]].fBody->R);
				}
			}
		}
	}
}

/// \brief	Update the IBM markers using new FEM node vales
///
///
/// \param level current grid level
void ObjectManager::fem_updateIBMarkers(int level) {


//	// Get the positions
//	std::vector<double> uVec(DOFsPerElement, 0.0);
//	std::vector<double> uRes(L_DIMS, 0.0);
//	double Ux1, Uy1, Ux2, Uy2;
//	double Ux, Uy;
//	double angle;
//	int elID;
//
//	// Loop through all IBM nodes
//	for (int node = 0; node < IBPointer->IBnode.size(); node++) {
//
//		// Get owning element and angle
//		elID = IBPointer->IBnode[node].FEMparent.elementID;
//		angle = FEMelement[elID].angle;
//
//		// The displacements in global coordinates
//		Ux1 = U[nodeDOFs*elID];
//		Uy1 = U[nodeDOFs*elID+1];
//		Ux2 = U[nodeDOFs*(elID+1)];
//		Uy2 = U[nodeDOFs*(elID+1)+1];
//
//		// Convert to local coordinates
//		uVec[0] = Ux1 * cos(angle) + Uy1 * sin(angle);
//		uVec[1] = -Ux1 * sin(angle) + Uy1 * cos(angle);
//		uVec[2] = U[nodeDOFs*elID+2];
//		uVec[3] = Ux2 * cos(angle) + Uy2 * sin(angle);
//		uVec[4] = -Ux2 * sin(angle) + Uy2 * cos(angle);
//		uVec[5] = U[nodeDOFs*(elID+1)+2];
//
//		// Get velocities in local coordinates
//		shapeFunctions(uVec, IBPointer->IBnode[node].FEMparent.zeta, FEMelement[elID].length, uRes);
//
//		// Convert back to global
//		Ux = uRes[0] * cos(angle) - uRes[1] * sin(angle);
//		Uy = uRes[0] * sin(angle) + uRes[1] * cos(angle);
}

// ****************************************************************************
