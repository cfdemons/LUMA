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
	double forceScale, angle, length, a, b;
	double N, V, N1, V1, M1, N2, V2, M2;
	int IBnode;

	// Loop through all bodies that this rank owns
	for (int ib = 0; ib < IdxFEM.size(); ib++) {

		// Check if body is on this level
		if (iBody[IdxFEM[ib]]._Owner->level == level) {

			// Set R vector to zero
			fill(iBody[IdxFEM[ib]].fBody->R.begin(), iBody[IdxFEM[ib]].fBody->R.end(), 0.0);

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
				angle = iBody[IdxFEM[ib]].fBody->elements[el].angles;
				length = iBody[IdxFEM[ib]].fBody->elements[el].length;

				// Now loop through all child IB nodes this element has	// TODO Vectorise this
				for (int node = 0; node < iBody[IdxFEM[ib]].fBody->elements[el].IBChildNodes.size(); node++) {

					// Get IB node and integration ranges
					IBnode = iBody[IdxFEM[ib]].fBody->elements[el].IBChildNodes[node].nodeID;
					a = iBody[IdxFEM[ib]].fBody->elements[el].IBChildNodes[node].zeta1;
					b = iBody[IdxFEM[ib]].fBody->elements[el].IBChildNodes[node].zeta2;

					// Convert force to local coordinates
					N = IBMforce_xyz[IBnode][eXDirection] * cos(angle) + IBMforce_xyz[IBnode][eYDirection] * sin(angle);
					V = -IBMforce_xyz[IBnode][eXDirection] * sin(angle) + IBMforce_xyz[IBnode][eYDirection] * cos(angle);

					// Get the nodal values by integrating over range of IB point
					N1 = 0.5 * N * length * (0.5 * b - 0.5 * a + 0.25 * SQ(a) - 0.25 * SQ(b));
					N2 = 0.5 * N * length * (-0.25 * SQ(a) + 0.25 * SQ(b) + 0.5 * b - 0.5 * a);
					V1 = 0.5 * V * length * (0.5 * b - 0.5 * a - SQ(a) * SQ(a) / 16.0 + SQ(b) * SQ(b) / 16.0 + 3.0 * SQ(a) / 8.0 - 3.0 * SQ(b) / 8.0);
					V2 = 0.5 * V * length * (0.5 * b - 0.5 * a + SQ(a) * SQ(a) / 16.0 - SQ(b) * SQ(b) / 16.0 - 3.0 * SQ(a) / 8.0 + 3.0 * SQ(b) / 8.0);
					M1 = 0.5 * V * length * (length * (-SQ(a) * SQ(a) + SQ(b) * SQ(b)) / 32.0 - length * (-TH(a) + TH(b)) / 24.0 - length * (-SQ(a) + SQ(b)) / 16.0 + length * (b - a) / 8.0);
					M2 = 0.5 * V * length * (length * (-SQ(a) * SQ(a) + SQ(b) * SQ(b)) / 32.0 + length * (-TH(a) + TH(b)) / 24.0 - length * (-SQ(a) + SQ(b)) / 16.0 - length * (b - a) / 8.0);

					// Add to load vector (and convert back to global coordinates)
					iBody[IdxFEM[ib]].fBody->R[el*iBody[IdxFEM[ib]].fBody->DOFsPerNode] += N1 * cos(angle) - V1 * sin(angle);
					iBody[IdxFEM[ib]].fBody->R[el*iBody[IdxFEM[ib]].fBody->DOFsPerNode+1] += N1 * sin(angle) + V1 * cos(angle);
					iBody[IdxFEM[ib]].fBody->R[el*iBody[IdxFEM[ib]].fBody->DOFsPerNode+2] += M1;
					iBody[IdxFEM[ib]].fBody->R[(el+1)*iBody[IdxFEM[ib]].fBody->DOFsPerNode] += N2 * cos(angle) - V2 * sin(angle);
					iBody[IdxFEM[ib]].fBody->R[(el+1)*iBody[IdxFEM[ib]].fBody->DOFsPerNode+1] += N2 * sin(angle) + V2 * cos(angle);
					iBody[IdxFEM[ib]].fBody->R[(el+1)*iBody[IdxFEM[ib]].fBody->DOFsPerNode+2] += M2;
				}
			}
		}
	}
}

// ****************************************************************************
