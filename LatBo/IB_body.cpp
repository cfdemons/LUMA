#include "stdafx.h"
#include "IB_body.h"
#include "definitions.h"
#include "ops_generic.h"
#include <math.h>

// ***************************************************************************************************
// Default constructor and destructor
IB_body::IB_body(void) { }
IB_body::~IB_body(void) { }

// ***************************************************************************************************
// Method to add marker
void IB_body::addMarker(double x, double y, double z, bool flex_rigid) {
	
	// Extend array of particles by 1 and construct a new IB_marker object
	markers.emplace_back( x, y, z, flex_rigid );

}

// ***************************************************************************************************
// Method to build a sphere / circle
void IB_body::makeBody(double radius, std::vector<double> centre, bool flex_rigid) {

	// Work out circumference
	double circum = 2 * PI * radius;

#if (dims == 3)

	// Sphere //

	// Following code for point generation on unit sphere actually seeds
	// using Fibonacci sphere technique. Code is not my own but works.
	double inc = PI * (3 - sqrt(5));
    double off = 2.0 / (float)Lspace ;
    for (unsigned int k = 0; k < Lspace; k++) { 
        double y = k * off - 1 + (off / 2); 
        double r = sqrt(1 - y*y); 
        double phi = k * inc;
		
		// Add Lagrange marker to body (scale by radius)
        addMarker(centre[0] + (cos(phi)*r * radius), y*radius + centre[1], centre[2] + (sin(phi)*r*radius), flex_rigid);
	}

	// Spacing (assuming all Lagrange markers are uniformly spaced)
	std::vector<double> diff;
	for (unsigned int d = 0; d < dims; d++) {
		diff.push_back ( markers[1].position[d] - markers[0].position[d] );
	}
	spacing = vecnorm( diff );

	
	
#else
	// Circle -- find theta
	std::vector<double> theta = linspace(0, 2*PI - (2*PI / Lspace), Lspace);
	for (size_t i = 0; i < theta.size(); i++) {

		// Add Lagrange marker to body
		addMarker(	centre[0] + radius * cos(theta[i]),
					centre[1] + radius * sin(theta[i]),
					0.0, flex_rigid	);
	}

	// Spacing
	std::vector<double> diff;
	for (unsigned int d = 0; d < dims; d++) {
		diff.push_back ( markers[1].position[d] - markers[0].position[d] );
	}
	spacing = vecnorm( diff );
#endif

}

// ***************************************************************************************************
// Method to build a cuboid/rectangle
void IB_body::makeBody(std::vector<double> width_length_depth, std::vector<double> centre, bool flex_rigid) {
	
#if (dims == 3)

	// Cube //

	// Check side lengths to make sure we can ensure points on the corners
	if (fmod(ibb_w,ibb_l) != 0 && fmod(ibb_l,ibb_w) != 0 && fmod(ibb_l,ibb_d) != 0 && fmod(ibb_d,ibb_l) != 0) {
			std::cout << "IB body cannot be built with uniform points. Change its dimensions. Exiting." << std::endl;
			system("pause");
			exit(EXIT_FAILURE);
		}
	
	// Get ratio of sides and degree of point refinement
	int side_ratio_1, side_ratio_2;
	double smallest_side;
	if (ibb_w > ibb_l) {
		if (ibb_d > ibb_l) { // Length is smallest
			smallest_side = ibb_l;
			side_ratio_1 = (int)(ibb_w / ibb_l);
			side_ratio_2 = (int)(ibb_d / ibb_l);
		} else { // Depth is smallest
			smallest_side = ibb_d;
			side_ratio_1 = (int)(ibb_w / ibb_d);
			side_ratio_2 = (int)(ibb_l / ibb_d);
		}
	} else if (ibb_d > ibb_w) { // Width is smallest
		smallest_side = ibb_w;
		side_ratio_1 = (int)(ibb_l / ibb_w);
		side_ratio_2 = (int)(ibb_d / ibb_w);
	} else { // Depth is smallest
		smallest_side = ibb_d;
		side_ratio_1 = (int)(ibb_w / ibb_d);
		side_ratio_2 = (int)(ibb_l / ibb_d);
	}

	// The following is derived by inspection of the problem and ensures a point on each corner of the body
	int ref = (int)floor( log( (Lspace + 4) / (4 * (1+side_ratio_1+side_ratio_2)) ) / log(2) );

	

	// Check to see if enough points
	if (ref == 0) {
		// Advisory of number of points
		unsigned int advisory_num_points = (unsigned int)(8 + (4 * (pow(2,1) -1) ) + 
		(4 * ( (side_ratio_1 * pow(2,1)) -1) ) + (4 * ( (side_ratio_2 * pow(2,1)) -1) ) );
		std::cout << "IB body does not have enough points. Need " << advisory_num_points << " to build body. Exiting." << std::endl;
		system("pause");
		exit(EXIT_FAILURE);
	}

	// Number of points required to get uniform distribution and points on corners
	unsigned int num_points = (unsigned int)(8 + (4 * (pow(2,ref) -1) ) + 
		(4 * ( (side_ratio_1 * pow(2,ref)) -1) ) + (4 * ( (side_ratio_2 * pow(2,ref)) -1) ) );

	// Spacing from smallest side
	spacing = smallest_side / pow(2,ref);
	for (double x = centre[0]-width_length_depth[0]/2; x <= centre[0]+width_length_depth[0]/2; x+=spacing) {
		for (double y = centre[1]-width_length_depth[1]/2; y <= centre[1]+width_length_depth[1]/2; y+=spacing) {
			for (double z = centre[2]-width_length_depth[2]/2; z <= centre[2]+width_length_depth[2]/2; z+=spacing) {

				// Only add if on surface so check to see if interior point
				if (	(x > centre[0]-width_length_depth[0]/2 && x < centre[0]+width_length_depth[0]/2) &&
						(y > centre[1]-width_length_depth[1]/2 && y < centre[1]+width_length_depth[1]/2) &&
						(z > centre[2]-width_length_depth[2]/2 && z < centre[2]+width_length_depth[2]/2)	) {
				} else {
					addMarker(x,y,z,flex_rigid);
				}
			}
		}
	}
	
#else
	// Square //

	// Check side lengths to make sure we can ensure points on the corners
	if ((fmod(ibb_w,ibb_l) != 0) && (fmod(ibb_l,ibb_w) != 0)) {
			std::cout << "IB body cannot be built with uniform points. Change its dimensions. Exiting." << std::endl;
			system("pause");
			exit(EXIT_FAILURE);
		}
	
	// Get ratio of sides and degree of point refinement
	int side_ratio;
	if (ibb_w > ibb_l) { // Length is shortest
		side_ratio = (int)(ibb_w / ibb_l);
	} else { // Width is shortest
		side_ratio = (int)(ibb_l / ibb_w);
	}

	// The following is derived by inspection of the problem and ensures a point on each corner of the body
	int ref = (int)floor( log( Lspace / (2 * (1+side_ratio)) ) / log(2) );

	// Check to see if enough points
	if (ref == 0) {
		// Advisory of number of points
		unsigned int advisory_num_points = (unsigned int)(unsigned int)(4 + (2 * (pow(2,1) -1) ) + (2 * ( (side_ratio * pow(2,1)) -1) ) );
		std::cout << "IB body does not have enough points. Need " << advisory_num_points << " to build body. Exiting." << std::endl;
		system("pause");
		exit(EXIT_FAILURE);
	}

	// Number of points required to get uniform distribution and points on corners
	unsigned int num_points = (unsigned int)(4 + (2 * (pow(2,ref) -1) ) + (2 * ( (side_ratio * pow(2,ref)) -1) ));

	// Find spacing of nodes
	spacing = (2 * (ibb_w + ibb_l) ) / num_points;
	for (double x = centre[0]-width_length_depth[0]/2; x <= centre[0]+width_length_depth[0]/2; x+=spacing) {
		for (double y = centre[1]-width_length_depth[1]/2; y <= centre[1]+width_length_depth[1]/2; y+=spacing) {
			double z = (b_z - a_z)/2;
			// Only add if on surface so check to see if interior point
			if (	(x > centre[0]-width_length_depth[0]/2 && x < centre[0]+width_length_depth[0]/2) &&
					(y > centre[1]-width_length_depth[1]/2 && y < centre[1]+width_length_depth[1]/2)	) {

			} else {
				addMarker(x,y,z,flex_rigid);
			}
		}
	}

#endif

}

// ***************************************************************************************************