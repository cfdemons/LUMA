/*
 * --------------------------------------------------------------
 *
 * ------ Lattice Boltzmann @ The University of Manchester ------
 *
 * -------------------------- L-U-M-A ---------------------------
 *
 *  Copyright (C) 2015, 2016
 *  E-mail contact: info@luma.manchester.ac.uk
 *
 * This software is for academic use only and not available for
 * distribution without written consent.
 *
 */

/* This file contains the body wrapper building methods for GridObj which in turn call
 * a suitable IBBody.makebody constructor.
 */

#include "../inc/stdafx.h"
#include "../inc/ObjectManager.h"

//// ************************************************************************** //
///// \brief Builds a prefab immersed boundary body.
///// \param	body_type	type of body to be built.
//void ObjectManager::ibm_buildBody(int body_type) {

//	// Declarations
//	std::vector<double> dimensions, centrepoint, start_position, end_position;
//
//	// Allocate the owning grid level and region for this immersed boundary
//	GridObj* g = NULL;
//	GridUtils::getGrid(this->_Grids, L_IB_ON_LEV, L_IB_ON_REG, g);
//
//	// Increase the iBody array by single IBBody object
//	iBody.emplace_back(g, iBody.size());
//
//	// Check that the owning grid was found and exists
//	if ( iBody.back()._Owner == NULL ) {
//		L_ERROR("Could not find the subgrid where the immersed boundary is meant to lie.", GridUtils::logfile);
//	}
//
//
//	if (body_type == 1) {
//		// =========================== Build a rectangle/cuboid ===========================  //
//
//		// Dimensions
//		dimensions.push_back(L_IBB_W);
//		dimensions.push_back(L_IBB_L);
//		dimensions.push_back(L_IBB_D);
//		centrepoint.push_back(L_IBB_X);
//		centrepoint.push_back(L_IBB_Y);
//		centrepoint.push_back(L_IBB_Z);
//
//		// Angle vector
//		std::vector<double> angles;
//		angles.push_back(L_IBB_ANGLE_VERT);
//#if (L_DIMS == 3)
//		angles.push_back(L_IBB_ANGLE_HORZ);
//#else
//		angles.push_back(0.0);
//#endif
//
//		// Build
//		iBody.back().makeBody(dimensions, angles, centrepoint, false, L_IBB_MOVABLE, static_cast<int>(iBody.size() - 1));
//
//
//	} else if (body_type == 2) {
//		// =========================== Build a circle/sphere ===========================  //
//
//		// Dimensions
//		centrepoint.push_back(L_IBB_X);
//		centrepoint.push_back(L_IBB_Y);
//		centrepoint.push_back(L_IBB_Z);
//
//		// Build
//		iBody.back().makeBody(L_IBB_R, centrepoint, false, L_IBB_MOVABLE, static_cast<int>(iBody.size() - 1));
//
//
//	} else if (body_type == 3) {
//		// =========================== Build both with custom dimensions ===========================  //
//
//		// Dimensions
//		centrepoint.push_back( 2*(double)L_BX/5 );
//		centrepoint.push_back( 2*(double)L_BY/5 );
//		centrepoint.push_back(L_IBB_Z);
//		dimensions.push_back(L_IBB_W);
//		dimensions.push_back(L_IBB_L);
//		dimensions.push_back(L_IBB_D);
//
//		// Build circle first
//		iBody.back().makeBody(L_IBB_R, centrepoint, false, L_IBB_MOVABLE, static_cast<int>(iBody.size() - 1));
//
//		// Grow vector
//		iBody.emplace_back(g, iBody.size());
//
//		// Angle vector
//		std::vector<double> angles;
//		angles.push_back(L_IBB_ANGLE_VERT);
//#if (L_DIMS == 3)
//		angles.push_back(L_IBB_ANGLE_HORZ);
//#else
//		angles.push_back(0.0);
//#endif
//		// Position
//		centrepoint[0] = 3*(double)L_BX/5;
//		centrepoint[1] = 3*(double)L_BY/5;
//		centrepoint[2] = L_IBB_Z;
//
//		// Build rectangle
//		iBody.back().makeBody(dimensions, angles, centrepoint, false, L_IBB_MOVABLE, static_cast<int>(iBody.size() - 1));
//
//
//	} else if (body_type == 4) {
//		// =========================== Build flexible filament ===========================  //
//
//		// End points
//		start_position.push_back(L_IBB_FILAMENT_START_X);
//		start_position.push_back(L_IBB_FILAMENT_START_Y);
//		start_position.push_back(L_IBB_FILAMENT_START_Z);
//
//		// Angle vector
//		std::vector<double> angles;
//		angles.push_back(L_IBB_ANGLE_VERT);
//#if (L_DIMS == 3)
//		angles.push_back(L_IBB_ANGLE_HORZ);
//#else
//		angles.push_back(0.0);
//#endif
//
//		// Pack BC info
//		std::vector<int> BCs;
//		BCs.push_back(L_FILAMENT_START_BC);
//		BCs.push_back(L_FILAMENT_END_BC);
//
//		// Build body
//		iBody.back().makeBody(L_NUM_MARKERS, start_position, L_IBB_FILAMENT_LENGTH, angles, BCs, true, true, static_cast<int>(iBody.size() - 1));
//
//
//	} else if (body_type == 5) {
//		// =========================== Build array of flexible filaments ===========================  //
//
//		std::vector<double> width_length, centrepoint, angles;
//
//
//		// If 2D, just build a filament
//#if (L_DIMS == 2)
//		*GridUtils::logfile << "Building plate as a filament as only 2D" << std::endl;
//		iBody.resize(iBody.size()-1); // Reduce size of iBody vector
//		ibm_buildBody(4);	// Call add again but as a filament
//#else
//		// In 3D, build an array of filaments
//		int group = 999;
//
//		// Centre point
//		centrepoint.push_back(L_IBB_X);
//		centrepoint.push_back(L_IBB_Y);
//		centrepoint.push_back(L_IBB_Z);
//
//		// Dimensions
//		width_length.push_back(L_IBB_W);
//		width_length.push_back(L_IBB_D);
//
//		// Angles
//		angles.push_back(L_IBB_ANGLE_VERT);
//		angles.push_back(L_IBB_ANGLE_HORZ);
//
//		// Pack BC info
//		std::vector<int> BCs;
//		BCs.push_back(L_FILAMENT_START_BC);
//		BCs.push_back(L_FILAMENT_END_BC);
//
//		// Declarations
//		std::vector<double> start_point; start_point.resize(3);
//
//		// Get length projected onto horizontal plane
//		double width_h = width_length[0]*cos(angles[0] * L_PI / 180);
//
//		// Define spacing for filaments
//		double  lengthwise_spacing;
//		// Must use odd number of filaments to ensure ones on edges and in centre
//		if (L_NUM_MARKERS % 2 == 0) {
//			lengthwise_spacing = (width_length[1]) / (L_NUM_MARKERS-1);
//		} else {
//			lengthwise_spacing = (width_length[1]) / L_NUM_MARKERS;
//		}
//
//		// Loop
//		int counter = (int)(width_length[1] / lengthwise_spacing);
//		for (int i = 0; i < counter; i++ ) {
//
//			// Get starting point
//			start_point[0] = (centrepoint[0] - (width_h*cos(angles[1]* L_PI / 180))/2 ) + i*( lengthwise_spacing*sin(angles[1]* L_PI / 180) );
//			start_point[1] = centrepoint[1] - (width_length[0] * sin(angles[0]* L_PI / 180))/2 ;
//			start_point[2] = (centrepoint[2] - (width_length[1]*cos(angles[1]* L_PI / 180))/2 ) + i*( lengthwise_spacing*cos(angles[1]* L_PI / 180) );
//
//			if (i != 0) {
//				// Increase the iBody array by single IBBody object when not on first loop
//				iBody.emplace_back(g, iBody.size());
//			}
//
//			// Build filament with fixed group ID making centre filament flexible
//			if (i == (counter-1)/2) {
//				iBody.back().makeBody(L_NUM_MARKERS, start_point, width_length[0], angles, BCs, true, true, group);
//			} else {
//				iBody.back().makeBody(L_NUM_MARKERS, start_point, width_length[0], angles, BCs, false, true, group);
//			}
//
//
//		}
//
//#endif
//
//
//	} else if (body_type == 6) {
//		// =========================== Build 2D rigid plate ===========================  //
//
//		// End points
//		start_position.push_back( L_IBB_FILAMENT_START_X - (L_IBB_FILAMENT_LENGTH/2)*cos(-L_IBB_ANGLE_VERT * L_PI / 180) );
//		start_position.push_back( L_IBB_FILAMENT_START_Y - (L_IBB_FILAMENT_LENGTH/2)*sin(-L_IBB_ANGLE_VERT * L_PI / 180) );
//		start_position.push_back( 0.0 );
//
//		// Angle vector
//		std::vector<double> angles;
//		angles.push_back(-L_IBB_ANGLE_VERT);
//		angles.push_back(0.0);
//
//		// Type of boundary condition at filament start:	0 == free; 1 = simply supported; 2 == clamped
//		// Pack BC info
//		std::vector<int> BCs;
//		BCs.push_back(L_FILAMENT_START_BC);
//		BCs.push_back(L_FILAMENT_END_BC);
//
//		// Build body as rigid filament
//		iBody.back().makeBody(L_NUM_MARKERS, start_position, L_IBB_FILAMENT_LENGTH, angles, BCs, false, false, static_cast<int>(iBody.size() - 1));
//
//
//	} else if (body_type == 7) {
//		// =========================== Build 2D rigid plate + flexible flap ===========================  //
//
//		// End points
//		start_position.push_back( L_IBB_FILAMENT_START_X - (L_IBB_FILAMENT_LENGTH/2)*cos(-L_IBB_ANGLE_VERT * L_PI / 180) );
//		start_position.push_back( L_IBB_FILAMENT_START_Y - (L_IBB_FILAMENT_LENGTH/2)*sin(-L_IBB_ANGLE_VERT * L_PI / 180) );
//		start_position.push_back( 0.0 );
//
//		// Angle vector
//		std::vector<double> angles;
//		angles.push_back(-L_IBB_ANGLE_VERT);
//		angles.push_back(0.0);
//
//		// Type of boundary condition at filament start:	0 == free; 1 = simply supported; 2 == clamped
//		// Pack BC info
//		std::vector<int> BCs;
//		BCs.push_back(L_FILAMENT_START_BC);
//		BCs.push_back(L_FILAMENT_END_BC);
//
//		// Build body with bit chopped off end to allow flap
//		double lspace = L_IBB_FILAMENT_LENGTH / (L_NUM_MARKERS-2);
//		iBody.back().makeBody(L_NUM_MARKERS - 1, start_position, L_IBB_FILAMENT_LENGTH - lspace, angles, BCs, false, false, static_cast<int>(iBody.size() - 1));
//
//		// Add flap of half the length (slightly lower resolution)
//		iBody.emplace_back(g, iBody.size());
//		angles[0] = 0.0;
//		start_position[0] = ( L_IBB_FILAMENT_START_X - (L_IBB_FILAMENT_LENGTH/2)*cos(-L_IBB_ANGLE_VERT * L_PI / 180) ) + (L_IBB_FILAMENT_LENGTH)*cos(-L_IBB_ANGLE_VERT * L_PI / 180);
//		start_position[1] = ( L_IBB_FILAMENT_START_Y - (L_IBB_FILAMENT_LENGTH/2)*sin(-L_IBB_ANGLE_VERT * L_PI / 180) ) + (L_IBB_FILAMENT_LENGTH)*sin(-L_IBB_ANGLE_VERT * L_PI / 180);
//		iBody.back().makeBody((L_NUM_MARKERS / 2) - 3, start_position, L_IBB_FILAMENT_LENGTH / 2, angles, BCs, true, true, static_cast<int>(iBody.size() - 1));
//
//
//
//	} else if (body_type == 8) {
//		// =========================== Build 3D rigid plate ===========================  //
//
//		// Dimensions
//		dimensions.push_back(L_IBB_W);
//		dimensions.push_back(L_IBB_D);
//		centrepoint.push_back(L_IBB_X);
//		centrepoint.push_back(L_IBB_Y);
//		centrepoint.push_back(L_IBB_Z);
//
//		// Inclination
//		std::vector<double> angles;
//		angles.push_back(-L_IBB_ANGLE_VERT);
//
//		// Build using plate flag
//		iBody.back().makeBody(dimensions, -L_IBB_ANGLE_VERT, centrepoint, false, false, static_cast<int>(iBody.size() - 1), true);
//
//
//	} else if (body_type == 9) {
//		// =========================== Build 3D rigid plate + flexible flap ===========================  //
//
//		// Dimensions
//		dimensions.push_back(L_IBB_W);
//		dimensions.push_back(L_IBB_D);
//		centrepoint.push_back(L_IBB_X);
//		centrepoint.push_back(L_IBB_Y);
//		centrepoint.push_back(L_IBB_Z);
//
//		// Build a rigid plate
//		double lspace = iBody.back().makeBody(dimensions, -L_IBB_ANGLE_VERT, centrepoint, false, false, static_cast<int>(iBody.size() - 1), true);
//
//		// Add filament array of half the length starting horizontal
//		// Trailing edge
//		double start_z = centrepoint[2] - L_IBB_D/2,
//		start_x = centrepoint[0] + (L_IBB_W/2)*cos(-L_IBB_ANGLE_VERT * L_PI / 180),
//		start_y = centrepoint[1] + (L_IBB_W/2)*sin(-L_IBB_ANGLE_VERT * L_PI / 180);
//
//		// Delete TE markers
//		std::vector<int> logvec;
//		for (size_t m = 0; m < iBody.back().markers.size(); m++) {
//			if (fabs(iBody.back().markers[m].position[1] - start_y) < L_SMALL_NUMBER) {
//				logvec.push_back(static_cast<int>(m));
//			}
//		}
//		iBody.back().markers.erase(
//			iBody.back().markers.begin() + logvec[0],
//			iBody.back().markers.begin() + logvec[logvec.size()-1]+1 );
//
//
//		// Angles
//		std::vector<double> angles;
//		angles.push_back(0.0);
//		angles.push_back(0.0);
//
//		// Pack BC info
//		std::vector<int> BCs;
//		BCs.push_back(L_FILAMENT_START_BC);
//		BCs.push_back(L_FILAMENT_END_BC);
//
//		// Define number of filaments (slightly lower resolution)
//		int lengthwise_fils;
//		// Must use odd number of filaments to ensure ones on edges and in centre
//		lengthwise_fils = (int) (L_IBB_D / lspace);
//		if (lengthwise_fils % 2 == 0) {
//			lengthwise_fils = lengthwise_fils + 1;
//		}
//
//		// Starting points
//		std::vector<double> start_point; start_point.resize(3);
//		start_point[0] = start_x;
//		start_point[1] = start_y;
//
//		// Loop
//		for (int i = 0; i < lengthwise_fils; i++ ) {
//
//			// Increase the iBody array by single IBBody object
//			iBody.emplace_back(g, iBody.size());
//
//			// Start point in z
//			start_point[2] = start_z + i*lspace;
//
//			// Build filament with fixed group ID making centre filament flexible
//			if (i == (lengthwise_fils-1)/2) {
//				iBody.back().makeBody((L_NUM_MARKERS/2)-6, start_point, L_IBB_W/2, angles, BCs, true, true, 888);
//			} else {
//				iBody.back().makeBody((L_NUM_MARKERS/2)-6, start_point, L_IBB_W/2, angles, BCs, false, true, 888);
//			}
//
//
//		}
//
//
//	}

//}
// ************************************************************************** //
