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
#include "../inc/GridObj.h"
#include "../inc/ObjectManager.h"



// *****************************************************************************
///	\brief	Write out marker positions for debugging
///
///	\param	ib		index of IBBody
void ObjectManager::ibm_debug_markerPosition(int ib) {

	// Only do if there is data to write
	if (iBody[ib].validMarkers.size() > 0) {

		// Get rank safely
		int rank = GridUtils::safeGetRank();

		// Open file and write header
		std::ofstream bodyout;
		bodyout.open(GridUtils::path_str + "/IBbody_" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out");
		bodyout << "Marker\tx\ty\tz" << std::endl;

		// Loop through markers
		for (auto i : iBody[ib].validMarkers) {
			bodyout << iBody[ib].markers[i].id << "\t" << iBody[ib].markers[i].position[eXDirection] << "\t" << iBody[ib].markers[i].position[eYDirection] << "\t" << iBody[ib].markers[i].position[eZDirection] << std::endl;
		}
		bodyout.close();
	}
}

// *****************************************************************************
///	\brief	Write out support point info
///
///	\param	ib		index of IBBody
void ObjectManager::ibm_debug_supportInfo(int ib) {

	// Only do if there is data to write
	if (iBody[ib].validMarkers.size() > 0) {

		// Get rank safely
		int rank = GridUtils::safeGetRank();

		// Open file and write header
		std::ofstream supportout;
		supportout.open(GridUtils::path_str + "/IBsupport_" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out", std::ios::app);

		// Write out the first (nearest) support marker
		supportout << "Marker\tRank\ti\tj\tk\tX\tY\tZ\tdeltaVal" << std::endl;

		// Loop through markers and its support points
		for (auto m : iBody[ib].validMarkers) {
			for (int s = 0; s < iBody[ib].markers[m].deltaval.size(); s++) {

				// Write out info
				supportout
					<< iBody[ib].markers[m].id << "\t"
					<< iBody[ib].markers[m].support_rank[s] << "\t"
					<< iBody[ib].markers[m].supp_i[s] << "\t"
					<< iBody[ib].markers[m].supp_j[s] << "\t"
					<< iBody[ib].markers[m].supp_k[s] << "\t"
					<< iBody[ib].markers[m].supp_x[s] << "\t"
					<< iBody[ib].markers[m].supp_y[s] << "\t"
					<< iBody[ib].markers[m].supp_z[s] << "\t"
					<< iBody[ib].markers[m].deltaval[s] << std::endl;
			}
		}

		// Close file
		supportout.close();
	}
}

// *****************************************************************************
///	\brief	Write out epsilon info
///
///	\param	ib		index of IBBody
void ObjectManager::ibm_debug_epsilon(int ib) {

	// Only do if there is data to write
	if (iBody[ib].validMarkers.size() > 0) {

		// Get rank safely
		int rank = GridUtils::safeGetRank();

		// Open file and write header
		std::ofstream epout;
		epout.open(GridUtils::path_str + "/Epsilon_" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out",std::ios::app);
		epout << "NEW TIME STEP" << std::endl;
		epout << "Marker\tEpsilon\tds" << std::endl;

		// Loop through markers
		for (auto m : iBody[ib].validMarkers) {
			epout << iBody[ib].markers[m].id << "\t" << iBody[ib].markers[m].epsilon << "\t" << iBody[ib].markers[m].ds << std::endl;
		}
		epout << std::endl;
		epout.close();
	}
}




// *****************************************************************************
///	\brief	Write out interpolate velocities
///
///	\param	ib		index of IBBody
void ObjectManager::ibm_debug_interpVel(int ib) {

	// Only do if there is data to write
	if (iBody[ib].validMarkers.size() > 0) {

		// Get rank safely
		int rank = GridUtils::safeGetRank();

		// Open file and write header
		std::ofstream predout;
		predout.open(GridUtils::path_str + "/interpVel_" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out",std::ios::app);
		predout << "NEW TIME STEP" << std::endl;
		predout << "Marker\tVelX\tVelY\tVelZ" << std::endl;

		// Loop through markers
		for (auto m : iBody[ib].validMarkers) {
			predout << iBody[ib].markers[m].id << "\t" 
				<< iBody[ib].markers[m].interpMom[eXDirection] / iBody[ib].markers[m].interpRho << "\t" 
				<< iBody[ib].markers[m].interpMom[eYDirection] / iBody[ib].markers[m].interpRho << "\t"
#if (L_DIMS == 3)
				<< iBody[ib].markers[m].interpMom[eZDirection] / iBody[ib].markers[m].interpRho
#else
				<< "0.0"
#endif
				<< std::endl;
		}
		predout << std::endl;
		predout.close();
	}
}

// *****************************************************************************
///	\brief	Write out marker forces
///
///	\param	ib		index of IBBody
void ObjectManager::ibm_debug_markerForce(int ib) {

	// Only do if there is data to write
	if (iBody[ib].validMarkers.size() > 0) {

		// Get rank safely
		int rank = GridUtils::safeGetRank();

		// Open file and write header
		std::ofstream forceout;
		forceout.open(GridUtils::path_str + "/force_xyz_" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out",std::ios::app);
		forceout << "NEW TIME STEP" << std::endl;
		forceout << "Marker\tFx\tFy\tFz" << std::endl;

		// Loop through markers
		for (auto m : iBody[ib].validMarkers) {
			forceout << iBody[ib].markers[m].id << "\t" << iBody[ib].markers[m].force_xyz[eXDirection] << "\t" << iBody[ib].markers[m].force_xyz[eYDirection] << "\t"
#if (L_DIMS == 3)
				<< iBody[ib].markers[m].force_xyz[eZDirection]
#else
				<< "0"
#endif
				<< std::endl;
		}
		forceout << std::endl;
		forceout.close();
	}
}


// *****************************************************************************
///	\brief	Write out support point velocities
///
///	\param	ib		index of IBBody
void ObjectManager::ibm_debug_supportVel(int ib) {

	// Only do if there is data to write
	if (iBody[ib].validMarkers.size() > 0) {

		// Get rank safely
		int rank = GridUtils::safeGetRank();

		// Open file and write header
		std::ofstream testout;
		testout.open(GridUtils::path_str + "/velSupp" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out", std::ios::app);
		testout << "NEW TIME STEP" << std::endl;

		// Get size of grid
		size_t M_lim = iBody[ib]._Owner->M_lim;
#if (L_DIMS == 3)
		size_t K_lim = iBody[ib]._Owner->K_lim;
#endif

		// Loop through markers and support sites
		for (auto m : iBody[ib].validMarkers) {
			for (int s = 0; s < iBody[ib].markers[m].deltaval.size(); s++) {

				// Write out the first (nearest) support marker
				if (m == 0 && s == 0)
					testout << "Marker\tSupport\tVelX\tVelY\tVelZ" << std::endl;

				// Only write out the support sites that this rank owns
				if (rank == iBody[ib].markers[m].support_rank[s]) {
					testout << m << "\t"
							<< s << "\t"
#if (L_DIMS == 2)
							<< iBody[ib]._Owner->u(iBody[ib].markers[m].supp_i[s], iBody[ib].markers[m].supp_j[s], eXDirection, M_lim, L_DIMS) << "\t"
							<< iBody[ib]._Owner->u(iBody[ib].markers[m].supp_i[s], iBody[ib].markers[m].supp_j[s], eYDirection, M_lim, L_DIMS) << std::endl;
#elif (L_DIMS == 3)
							<< iBody[ib]._Owner->u(iBody[ib].markers[m].supp_i[s], iBody[ib].markers[m].supp_j[s], iBody[ib].markers[m].supp_k[s], eXDirection, M_lim, K_lim, L_DIMS) << "\t"
							<< iBody[ib]._Owner->u(iBody[ib].markers[m].supp_i[s], iBody[ib].markers[m].supp_j[s], iBody[ib].markers[m].supp_k[s], eYDirection, M_lim, K_lim, L_DIMS) << "\t"
							<< iBody[ib]._Owner->u(iBody[ib].markers[m].supp_i[s], iBody[ib].markers[m].supp_j[s], iBody[ib].markers[m].supp_k[s], eZDirection, M_lim, K_lim, L_DIMS) << std::endl;
#endif
				}
			}
		}
		testout << std::endl;
		testout.close();
	}
}


// *****************************************************************************
///	\brief	Write out support point forces
///
///	\param	ib		index of IBBody
void ObjectManager::ibm_debug_supportForce(int ib) {


	// Only do if there is data to write
	if (iBody[ib].validMarkers.size() > 0) {

		// Get rank safely
		int rank = GridUtils::safeGetRank();

		// Open file and write header
		std::ofstream testout;
		testout.open(GridUtils::path_str + "/force_xyz_supp" + std::to_string(iBody[ib].id) + "_rank" + std::to_string(rank) + ".out", std::ios::app);
		testout << "NEW TIME STEP" << std::endl;

		// Get size of grid
		size_t M_lim = iBody[ib]._Owner->M_lim;
		size_t K_lim = iBody[ib]._Owner->K_lim;

		// Loop through markers and support sites
		for (auto m : iBody[ib].validMarkers) {
			for (int s = 0; s < iBody[ib].markers[m].deltaval.size(); s++) {

				// Write out the first (nearest) support marker
				if (m == 0 && s == 0)
					testout << "Marker\tSupport\tFx\t\tFy\t\tFz" << std::endl;

				// Only write out the support sites that this rank owns
				if (rank == iBody[ib].markers[m].support_rank[s]) {
					testout << m << "\t"
							<< s << "\t"
							<< iBody[ib]._Owner->force_xyz(iBody[ib].markers[m].supp_i[s], iBody[ib].markers[m].supp_j[s], iBody[ib].markers[m].supp_k[s], eXDirection, M_lim, K_lim, L_DIMS) << "\t"
							<< iBody[ib]._Owner->force_xyz(iBody[ib].markers[m].supp_i[s], iBody[ib].markers[m].supp_j[s], iBody[ib].markers[m].supp_k[s], eYDirection, M_lim, K_lim, L_DIMS) << "\t"
							<< iBody[ib]._Owner->force_xyz(iBody[ib].markers[m].supp_i[s], iBody[ib].markers[m].supp_j[s], iBody[ib].markers[m].supp_k[s], eZDirection, M_lim, K_lim, L_DIMS) << std::endl;
				}
			}
		}
		testout << std::endl;
		testout.close();
	}
}
