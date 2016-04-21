/* This file contains the body building methods for GridObj which in turn call a suitable IBBody.makebody constructor.
*/

#include "../inc/stdafx.h"
#include "../inc/ObjectManager.h"
#include "../inc/definitions.h"

// ***************************************************************************************************
// Method to call body building routine for IBM
void ObjectManager::ibm_build_body(int body_type) {

	// Declarations
	std::vector<double> dimensions, centrepoint, start_position, end_position;

	// Allocate the owning grid level and region for this immersed boundary
	GridObj* g = NULL;
	GridUtils::getGrid(this->_Grids, IB_Lev, IB_Reg, g);

	// Increase the iBody array by single IBBody object
	iBody.emplace_back(g);

	// Check that the owning grid was found and exists
	if ( iBody.back()._Owner == NULL ) {
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "Could not find the subgrid where the immersed boundary is meant to lie." << std::endl;
		exit(LATBO_FAILED);
	}


	if (body_type == 1) {
		// =========================== Build a rectangle/cuboid ===========================  //

		// Dimensions
		dimensions.push_back(ibb_w);
		dimensions.push_back(ibb_l);
		dimensions.push_back(ibb_d);
		centrepoint.push_back(ibb_x);
		centrepoint.push_back(ibb_y);
		centrepoint.push_back(ibb_z);

		// Angle vector
		std::vector<double> angles;
		angles.push_back(ibb_angle_vert);
#if (dims == 3)
		angles.push_back(ibb_angle_horz);
#else
		angles.push_back(0.0);
#endif

		// Build
		iBody.back().makeBody(dimensions,angles,centrepoint,false,ibb_deform,iBody.size()-1);


	} else if (body_type == 2) {
		// =========================== Build a circle/sphere ===========================  //

		// Dimensions
		centrepoint.push_back(ibb_x);
		centrepoint.push_back(ibb_y);
		centrepoint.push_back(ibb_z);

		// Build
		iBody.back().makeBody(ibb_r,centrepoint,false,ibb_deform,iBody.size()-1);


	} else if (body_type == 3) {
		// =========================== Build both with custom dimensions ===========================  //

		// Dimensions
		centrepoint.push_back( 2*(double)(b_x-a_x)/5 );
		centrepoint.push_back( 2*(double)(b_y - a_y)/5 );
		centrepoint.push_back(ibb_z);
		dimensions.push_back(ibb_w);
		dimensions.push_back(ibb_l);
		dimensions.push_back(ibb_d);

		// Build circle first
		iBody.back().makeBody(ibb_r,centrepoint,false,ibb_deform,iBody.size()-1);
		
		// Grow vector
		iBody.emplace_back();

		// Angle vector
		std::vector<double> angles;
		angles.push_back(ibb_angle_vert);
#if (dims == 3)
		angles.push_back(ibb_angle_horz);
#else
		angles.push_back(0.0);
#endif
		// Position
		centrepoint[0] = 3*(double)(b_x - a_x)/5;
		centrepoint[1] = 3*(double)(b_y - a_y)/5;
		centrepoint[2] = ibb_z;

		// Build rectangle
		iBody.back().makeBody(dimensions,angles,centrepoint,false,ibb_deform,iBody.size()-1);


	} else if (body_type == 4) {
		// =========================== Build flexible filament ===========================  //

		// End points
		start_position.push_back(ibb_start_x);
		start_position.push_back(ibb_start_y);
		start_position.push_back(ibb_start_z);

		// Angle vector
		std::vector<double> angles;
		angles.push_back(ibb_angle_vert);
#if (dims == 3)
		angles.push_back(ibb_angle_horz);
#else
		angles.push_back(0.0);
#endif

		// Pack BC info
		std::vector<int> BCs;
		BCs.push_back(start_BC);
		BCs.push_back(end_BC);

		// Build body
		iBody.back().makeBody(num_markers, start_position, ibb_length, angles, BCs, true, true, iBody.size()-1);

		
	} else if (body_type == 5) {
		// =========================== Build array of flexible filaments ===========================  //

		std::vector<double> width_length, centrepoint, angles;


		// If 2D, just build a filament
#if (dims == 2)
		*GridUtils::logfile << "Building plate as a filament as only 2D" << std::endl;
		iBody.resize(iBody.size()-1); // Reduce size of iBody vector
		ibm_build_body(4);	// Call add again but as a filament
#else
		// In 3D, build an array of filaments
		int group = 999;

		// Centre point
		centrepoint.push_back(ibb_x);
		centrepoint.push_back(ibb_y);
		centrepoint.push_back(ibb_z);		

		// Dimensions
		width_length.push_back(ibb_w);
		width_length.push_back(ibb_d);

		// Angles
		angles.push_back(ibb_angle_vert);
		angles.push_back(ibb_angle_horz);

		// Pack BC info
		std::vector<int> BCs;
		BCs.push_back(start_BC);
		BCs.push_back(end_BC);

		// Declarations
		std::vector<double> start_point; start_point.resize(3);

		// Get length projected onto horizontal plane
		double width_h = width_length[0]*cos(angles[0] * PI / 180);
	
		// Define spacing for filaments
		double  lengthwise_spacing;
		// Must use odd number of filaments to ensure ones on edges and in centre
		if (num_markers % 2 == 0) {
			lengthwise_spacing = (width_length[1]) / (num_markers-1);
		} else {
			lengthwise_spacing = (width_length[1]) / num_markers;
		}

		// Loop
		int counter = (int)(width_length[1] / lengthwise_spacing);
		for (int i = 0; i < counter; i++ ) {

			// Get starting point
			start_point[0] = (centrepoint[0] - (width_h*cos(angles[1]* PI / 180))/2 ) + i*( lengthwise_spacing*sin(angles[1]* PI / 180) );
			start_point[1] = centrepoint[1] - (width_length[0] * sin(angles[0]* PI / 180))/2 ;
			start_point[2] = (centrepoint[2] - (width_length[1]*cos(angles[1]* PI / 180))/2 ) + i*( lengthwise_spacing*cos(angles[1]* PI / 180) );

			if (i != 0) {
				// Increase the iBody array by single IBBody object when not on first loop
				iBody.emplace_back();
			}

			// Build filament with fixed group ID making centre filament flexible
			if (i == (counter-1)/2) {
				iBody.back().makeBody(num_markers, start_point, width_length[0], angles, BCs, true, true, group);
			} else {
				iBody.back().makeBody(num_markers, start_point, width_length[0], angles, BCs, false, true, group);
			}


		}

#endif


	} else if (body_type == 6) {
		// =========================== Build 2D rigid plate ===========================  //

		// End points
		start_position.push_back( ibb_start_x - (ibb_length/2)*cos(-ibb_angle_vert * PI / 180) );
		start_position.push_back( ibb_start_y - (ibb_length/2)*sin(-ibb_angle_vert * PI / 180) );
		start_position.push_back( 0.0 );

		// Angle vector
		std::vector<double> angles;
		angles.push_back(-ibb_angle_vert);
		angles.push_back(0.0);

		// Type of boundary condition at filament start:	0 == free; 1 = simply supported; 2 == clamped
		// Pack BC info
		std::vector<int> BCs;
		BCs.push_back(start_BC);
		BCs.push_back(end_BC);

		// Build body as rigid filament
		iBody.back().makeBody(num_markers, start_position, ibb_length, angles, BCs, false, false, iBody.size()-1);


	} else if (body_type == 7) {
		// =========================== Build 2D rigid plate + flexible flap ===========================  //
	
		// End points
		start_position.push_back( ibb_start_x - (ibb_length/2)*cos(-ibb_angle_vert * PI / 180) );
		start_position.push_back( ibb_start_y - (ibb_length/2)*sin(-ibb_angle_vert * PI / 180) );
		start_position.push_back( 0.0 );

		// Angle vector
		std::vector<double> angles;
		angles.push_back(-ibb_angle_vert);
		angles.push_back(0.0);

		// Type of boundary condition at filament start:	0 == free; 1 = simply supported; 2 == clamped
		// Pack BC info
		std::vector<int> BCs;
		BCs.push_back(start_BC);
		BCs.push_back(end_BC);

		// Build body with bit chopped off end to allow flap
		double lspace = ibb_length / (num_markers-2);
		iBody.back().makeBody(num_markers-1, start_position, ibb_length - lspace, angles, BCs, false, false, iBody.size()-1);

		// Add flap of half the length (slightly lower resolution)
		iBody.emplace_back();
		angles[0] = 0.0;
		start_position[0] = ( ibb_start_x - (ibb_length/2)*cos(-ibb_angle_vert * PI / 180) ) + (ibb_length)*cos(-ibb_angle_vert * PI / 180);
		start_position[1] = ( ibb_start_y - (ibb_length/2)*sin(-ibb_angle_vert * PI / 180) ) + (ibb_length)*sin(-ibb_angle_vert * PI / 180);
		iBody.back().makeBody((num_markers/2)-3, start_position, ibb_length/2, angles, BCs, true, true, iBody.size()-1);

		

	} else if (body_type == 8) {
		// =========================== Build 3D rigid plate ===========================  //

		// Dimensions
		dimensions.push_back(ibb_w);
		dimensions.push_back(ibb_d);
		centrepoint.push_back(ibb_x);
		centrepoint.push_back(ibb_y);
		centrepoint.push_back(ibb_z);

		// Inclination
		std::vector<double> angles;
		angles.push_back(-ibb_angle_vert);

		// Build using plate flag
		iBody.back().makeBody(dimensions,-ibb_angle_vert,centrepoint,false,false,iBody.size()-1,true);


	} else if (body_type == 9) {
		// =========================== Build 3D rigid plate + flexible flap ===========================  //

		// Dimensions
		dimensions.push_back(ibb_w);
		dimensions.push_back(ibb_d);
		centrepoint.push_back(ibb_x);
		centrepoint.push_back(ibb_y);
		centrepoint.push_back(ibb_z);

		// Build a rigid plate
		double lspace = iBody.back().makeBody(dimensions,-ibb_angle_vert,centrepoint,false,false,iBody.size()-1,true);
		
		// Add filament array of half the length starting horizontal
		// Trailing edge
		double start_z = centrepoint[2] - ibb_d/2, 
		start_x = centrepoint[0] + (ibb_w/2)*cos(-ibb_angle_vert * PI / 180), 
		start_y = centrepoint[1] + (ibb_w/2)*sin(-ibb_angle_vert * PI / 180);

		// Delete TE markers
		std::vector<int> logvec;
		for (size_t m = 0; m < iBody.back().markers.size(); m++) {
			if (fabs(iBody.back().markers[m].position[1] - start_y) < 1e-6) {				
				logvec.push_back(m);
			}
		}
		iBody.back().markers.erase( 
			iBody.back().markers.begin() + logvec[0], 
			iBody.back().markers.begin() + logvec[logvec.size()-1]+1 );


		// Angles
		std::vector<double> angles;
		angles.push_back(0.0);
		angles.push_back(0.0);

		// Pack BC info
		std::vector<int> BCs;
		BCs.push_back(start_BC);
		BCs.push_back(end_BC);
	
		// Define number of filaments (slightly lower resolution)
		int lengthwise_fils;
		// Must use odd number of filaments to ensure ones on edges and in centre
		lengthwise_fils = (int) (ibb_d / lspace);
		if (lengthwise_fils % 2 == 0) {
			lengthwise_fils = lengthwise_fils + 1;
		}

		// Starting points
		std::vector<double> start_point; start_point.resize(3);
		start_point[0] = start_x;
		start_point[1] = start_y;

		// Loop
		for (int i = 0; i < lengthwise_fils; i++ ) {

			// Increase the iBody array by single IBBody object
			iBody.emplace_back();

			// Start point in z
			start_point[2] = start_z + i*lspace;

			// Build filament with fixed group ID making centre filament flexible
			if (i == (lengthwise_fils-1)/2) {
				iBody.back().makeBody((num_markers/2)-6, start_point, ibb_w/2, angles, BCs, true, true, 888);
			} else {
				iBody.back().makeBody((num_markers/2)-6, start_point, ibb_w/2, angles, BCs, false, true, 888);
			}


		}


	}

}

// ***************************************************************************************************
