/*
* --------------------------------------------------------------
*
* ------ Lattice Boltzmann @ The University of Manchester ------
*
* -------------------------- L-U-M-A ---------------------------
*
* Copyright 2019 The University of Manchester
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
* http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.*
*/

#include "../inc/stdafx.h"
#include "../inc/ObjectManager.h"
#include "../inc/IBInfo.h"


// *****************************************************************************
///	\brief	Build MPI-IBM comm classes
///
///	\param	level		current grid level
void ObjectManager::ibm_updateMPIComms(int level) {

	// Get MPI manager instance
	MpiManager *mpim = MpiManager::getInstance();

	// Build MPI comm classes for epsilon calculation
	mpim->mpi_buildMarkerComms(level);

	// Build MPI comm class for support communication
	mpim->mpi_buildSupportComms(level);
}


// *****************************************************************************
///	\brief	Pass velocity values from support site which exist off-rank
///
///	\param	level		current grid level
void ObjectManager::ibm_interpolateOffRankVels(int level) {

	// Get the mpi manager instance
	MpiManager *mpim = MpiManager::getInstance();

	// Perform interpolation communication
	std::vector<std::vector<double>> interpVels;
	mpim->mpi_interpolateComm(level, interpVels);

	// Create idx vector
	std::vector<int> idx(mpim->num_ranks, 0);

	// Now interpolate these remaining values onto the marker
	int ib, m, s, fromRank;
	for (int i = 0; i < mpim->supportCommMarkerSide[level].size(); i++) {

		// Get body idx
		ib = bodyIDToIdx[mpim->supportCommMarkerSide[level][i].bodyID];

		// Only do if body is on this grid level
		if (iBody[ib]._Owner->level == level) {

			// Get IDs of support site
			m = mpim->supportCommMarkerSide[level][i].markerIdx;
			s = mpim->supportCommMarkerSide[level][i].supportID;
			fromRank = mpim->supportCommMarkerSide[level][i].rankComm;

			// Interpolate density
			iBody[ib].markers[m].interpRho += interpVels[fromRank][idx[fromRank]] * iBody[ib].markers[m].deltaval[s] * iBody[ib].markers[m].local_area;

			// Shift index
			idx[fromRank]++;

			// Interpolate these values
			for (int dir = 0; dir < L_DIMS; dir++)
				iBody[ib].markers[m].interpMom[dir] += interpVels[fromRank][dir+idx[fromRank]] * iBody[ib].markers[m].deltaval[s] * iBody[ib].markers[m].local_area;

			// Shift index
			idx[fromRank] += L_DIMS;
		}
	}
}


// *****************************************************************************
///	\brief	Spread forces to support site which exist off-rank
///
///	\param	level		current grid level
void ObjectManager::ibm_spreadOffRankForces(int level) {

	// Get the mpi manager instance
	MpiManager *mpim = MpiManager::getInstance();

	// Perform interpolation communication
	std::vector<std::vector<double>> spreadForces;
	mpim->mpi_spreadComm(level, spreadForces);

	// Create idx vector
	std::vector<int> idx(mpim->num_ranks, 0);
	std::vector<int> suppIdx(3, 0);

	// Now interpolate these remaining values onto the marker
	int ib, fromRank;
	for (int i = 0; i < mpim->supportCommSupportSide[level].size(); i++) {

		// Get body idx
		ib = bodyIDToIdx[mpim->supportCommSupportSide[level][i].bodyID];

		// Only do if body is on this grid level
		if (iBody[ib]._Owner->level == level) {

			// Get grid sizes
			size_t M_lim = iBody[ib]._Owner->M_lim;
			size_t K_lim = iBody[ib]._Owner->K_lim;

			// Get IDs of support site
			fromRank = mpim->supportCommSupportSide[level][i].rankComm;
			suppIdx = mpim->supportCommSupportSide[level][i].supportIdx;

			// Interpolate these values
			for (int dir = 0; dir < L_DIMS; dir++)
				iBody[ib]._Owner->force_xyz(suppIdx[eXDirection], suppIdx[eYDirection], suppIdx[eZDirection], dir, M_lim, K_lim, L_DIMS) -=
						spreadForces[fromRank][idx[fromRank] + dir];

			// Shift index
			idx[fromRank] += L_DIMS;
		}
	}
}


// *****************************************************************************
///	\brief	Update new markers across all ranks
///
///	\param	level		current grid level
void ObjectManager::ibm_updateMarkers(int level) {

	// Get the mpi manager instance
	MpiManager *mpim = MpiManager::getInstance();

	// Loop through all flexible bodies that this rank owns
	for (auto ib : idxFEM) {

		// Only do if on this grid level
		if (iBody[ib]._Owner->level == level) {

			// Loop through all markers and assign rank
			for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
				iBody[ib].markers[m].owningRank = GridUtils::getRankfromPosition(iBody[ib].markers[m].position);
			}
		}
	}

	// Vectors of received markers
	std::vector<std::vector<int>> markerIDs(iBody.size(), std::vector<int>(0));
	std::vector<std::vector<std::vector<double>>> positions(iBody.size(), std::vector<std::vector<double>>(0, std::vector<double>(0)));
	std::vector<std::vector<std::vector<double>>> vels(iBody.size(), std::vector<std::vector<double>>(0, std::vector<double>(0)));

	// Do MPI comm for spreading markers
	mpim->mpi_spreadNewMarkers(level, markerIDs, positions, vels);

	// Loop through all iBodies
	for (size_t ib = 0; ib < iBody.size(); ib++) {

		// If body is on this level and flexible
		if (iBody[ib]._Owner->level == level && iBody[ib].isFlexible) {

			// Also if not owned by this rank
			if (iBody[ib].owningRank != mpim->my_rank) {

				// Markers to be deleted and inserted
				std::vector<int> deleteMarkers;
				std::vector<int> insertMarkers;
				std::vector<int> insertMarkersAt;

				// Loop through old markers
				for (size_t oldMarker = 0; oldMarker < iBody[ib].markers.size(); oldMarker++) {

					// Delete flag
					bool deleteFlag = true;

					// Loop through new markers
					for (size_t newMarker = 0; newMarker < markerIDs[ib].size(); newMarker++) {

						// Break early once new marker IDs are higher than old one
						if (markerIDs[ib][newMarker] > iBody[ib].markers[oldMarker].id) {
							break;
						}
						else if (iBody[ib].markers[oldMarker].id == markerIDs[ib][newMarker]) {

							// If IDs are same then just update position
							for (int d = 0; d < L_DIMS; d++) {
								iBody[ib].markers[oldMarker].position[d] = positions[ib][newMarker][d];
								iBody[ib].markers[oldMarker].markerVel[d] = vels[ib][newMarker][d];
							}
							deleteFlag = false;
							break;
						}
					}

					// If we get to hear on the last new marker then this old marker needs to be deleted
					if (deleteFlag == true)
						deleteMarkers.push_back(static_cast<int>(oldMarker));
				}

				// Now delete those markers
				for (size_t i = 0; i < deleteMarkers.size(); i++)
					iBody[ib].markers.erase(iBody[ib].markers.begin() + deleteMarkers[i] - i);

				// Now find the ones that need inserted
				for (size_t newMarker = 0; newMarker < markerIDs[ib].size(); newMarker++) {

					// Where to insert new marker
					int insertAt = 0;
					bool insertFlag = true;

					// Loop through old markers
					for (size_t oldMarker = 0; oldMarker < iBody[ib].markers.size(); oldMarker++) {

						// Increment the position to insert it
						if (iBody[ib].markers[oldMarker].id < markerIDs[ib][newMarker]) {
							insertAt++;
						}
						else if (iBody[ib].markers[oldMarker].id == markerIDs[ib][newMarker]) {
							insertFlag = false;
							break;
						} else {
							break;
						}
					}

					// If we get to hear on the last new marker then this old marker needs to be deleted
					if (insertFlag == true) {
						insertMarkers.push_back(static_cast<int>(newMarker));
						insertMarkersAt.push_back(insertAt);
					}
				}

				// Declare values
				double x = 0.0, y = 0.0, z = 0.0;

				// Now loop through and insert
				for (int i = 0; i < insertMarkers.size(); i++) {

					// Seperate into values
					x = positions[ib][insertMarkers[i]][eXDirection];
					y = positions[ib][insertMarkers[i]][eYDirection];
#if (L_DIMS == 3)
					z = positions[ib][insertMarkers[i]][eZDirection];
#endif

					// Create the marker to be inserted
					IBMarker insertMarker(x, y, z, markerIDs[ib][insertMarkers[i]], iBody[ib]._Owner);

					// Update velocity
					for (int d = 0; d < L_DIMS; d++)
						insertMarker.markerVel[d] = vels[ib][insertMarkers[i]][d];

					// Insert the marker
					iBody[ib].markers.insert(iBody[ib].markers.begin() + insertMarkersAt[i] + i, insertMarker);
				}
			}

			// Update valid markers
			iBody[ib].getValidMarkers();
		}
	}
}
