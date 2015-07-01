#include "stdafx.h"
#include "IB_body.h"
#include "definitions.h"
#include "ops_generic.h"
#include <math.h>
#include "GridObj.h"

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

	// Designate body as beiong flexible or rigid
	this->flex_rigid = flex_rigid;

	// Work out circumference
	double circum = 2 * PI * radius;

#if (dims == 3)

	// Sphere //

	// Following code for point generation on unit sphere actually seeds
	// using Fibonacci sphere technique. Code is not my own but works.
	double inc = PI * (3 - sqrt(5));
    double off = 2.0 / (float)num_markers ;
    for (unsigned int k = 0; k < num_markers; k++) { 
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
	std::vector<double> theta = linspace(0, 2*PI - (2*PI / num_markers), num_markers);
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
	
	// Designate body as being flexible or rigid
	this->flex_rigid = flex_rigid;

	// Shorter variable names for convenience
	double wid = width_length_depth[0];
	double len = width_length_depth[1];
	double dep = width_length_depth[2];

#if (dims == 3)

	// Cube //

	// Check side lengths to make sure we can ensure points on the corners
	if (fmod(ibb_w,ibb_l) != 0 && fmod(len,wid) != 0 && fmod(len,ibb_d) != 0 && fmod(dep,len) != 0) {
			std::cout << "IB body cannot be built with uniform points. Change its dimensions. Exiting." << std::endl;
			system("pause");
			exit(EXIT_FAILURE);
		}
	
	// Get ratio of sides and degree of point refinement
	int side_ratio_1, side_ratio_2;
	double smallest_side;
	if (wid > len) {
		if (dep > len) { // Length is smallest
			smallest_side = len;
			side_ratio_1 = (int)(wid / len);
			side_ratio_2 = (int)(dep / len);
		} else { // Depth is smallest
			smallest_side = dep;
			side_ratio_1 = (int)(wid / dep);
			side_ratio_2 = (int)(len / dep);
		}
	} else if (dep > wid) { // Width is smallest
		smallest_side = wid;
		side_ratio_1 = (int)(len / wid);
		side_ratio_2 = (int)(dep / wid);
	} else { // Depth is smallest
		smallest_side = dep;
		side_ratio_1 = (int)(wid / dep);
		side_ratio_2 = (int)(len / dep);
	}

	// The following is derived by inspection of the problem and ensures a point on each corner of the body
	// Since we have points on the faces as well as the edges we have a quadratic expression to be solved to
	// find the level of refinement required to get the nearest number of points to num_markers.
	int ref;

	int a = 2*side_ratio_1 + 2*side_ratio_1 + 2*side_ratio_1*side_ratio_2;
	int b = 4 + 4*side_ratio_1 + 4*side_ratio_2;
	int c = 8-num_markers;
	
	double P = (-b + sqrt(pow(b,2) - (4*a*c)) )/(2*a) + 1;
		
	ref = (int)ceil( log( P ) / log(2) );
		

	// Check to see if enough points
	if (ref == 0) {
		// Advisory of number of points
		unsigned int advisory_num_points = (unsigned int)(8 + 
			(4 * (pow(2,1) -1) ) + 
			(4 * ( (side_ratio_1 * pow(2,1)) -1) ) + 
			(4 * ( (side_ratio_2 * pow(2,1)) -1) ) +
			(2 * ( (pow(2,1) -1)*side_ratio_1 * (pow(2,1) -1) )) + 
			(2 * ( (pow(2,1) -1)*side_ratio_2 * (pow(2,1) -1) )) + 
			(2 * ( (pow(2,1) -1)*side_ratio_1 * (pow(2,1) -1)*side_ratio_2 ))
		);
		std::cout << "IB body does not have enough points. Need " << advisory_num_points << " to build body. Exiting." << std::endl;
		system("pause");
		exit(EXIT_FAILURE);
	}

	// Number of points required to get uniform distribution and points on corners
	unsigned int num_points = (unsigned int)(8 + 
			(4 * (pow(2,ref) -1) ) + 
			(4 * ( (side_ratio_1 * pow(2,ref)) -1) ) + 
			(4 * ( (side_ratio_2 * pow(2,ref)) -1) ) +
			(2 * ( (pow(2,ref) -1)*side_ratio_1 * (pow(2,ref) -1) )) + 
			(2 * ( (pow(2,ref) -1)*side_ratio_2 * (pow(2,ref) -1) )) + 
			(2 * ( (pow(2,ref) -1)*side_ratio_1 * (pow(2,ref) -1)*side_ratio_2 ))
		);

	// Spacing from smallest side
	spacing = smallest_side / pow(2,ref);
	
	// Start locations of point generator
	double start_x = centre[0]-wid/2, start_y = centre[1]-len/2, start_z = centre[2]-dep/2, x, y, z;
	
	// Number of points in each direction
	unsigned int np_x = (unsigned int)(len / spacing)+1, np_y = (unsigned int)(wid / spacing)+1, np_z = (unsigned int)(dep / spacing)+1;

	// Loop over matrix of points
	for (unsigned int i = 0; i < np_x; i++) {
		for (unsigned int j = 0; j < np_y; j++) {
			for (unsigned int k = 0; k < np_z; k++) {

				// x, y and z positions
				x = start_x + i*spacing;
				y = start_y + j*spacing;
				z = start_z + k*spacing;

				// Only add the marker if the point is on an edge
				if (	(i == 0 || i == np_x - 1) ||
						(j == 0 || j == np_y - 1) ||
						(k == 0 || k == np_z - 1)
					) {

					addMarker(x,y,z,flex_rigid);

				}

			}
		}
	}

	
#else
	// Square //

	// Check side lengths to make sure we can ensure points on the corners
	if ((fmod(wid,len) != 0) && (fmod(len,wid) != 0)) {
			std::cout << "IB body cannot be built with uniform points. Change its dimensions. Exiting." << std::endl;
			system("pause");
			exit(EXIT_FAILURE);
		}
	
	// Get ratio of sides and degree of point refinement
	int side_ratio;
	if (wid > len) { // Length is shortest
		side_ratio = (int)(wid / len);
	} else { // Width is shortest
		side_ratio = (int)(len / wid);
	}

	// The following is derived by inspection of the problem and ensures a point on each corner of the body
	// Double the number of points to be used as can be rounded down to only a relatively small number
	int ref = (int)ceil( log( num_markers / (2 * (1+side_ratio)) ) / log(2) );

	// Check to see if enough points
	if (ref == 0) {
		// Advisory of number of points
		unsigned int advisory_num_points = (unsigned int)(4 + (2 * (pow(2,1) -1) ) + (2 * ( (side_ratio * pow(2,1)) -1) ) );
		std::cout << "IB body does not have enough points. Need " << advisory_num_points << " to build body. Exiting." << std::endl;
		system("pause");
		exit(EXIT_FAILURE);
	}

	// Number of points required to get uniform distribution and points on corners
	unsigned int num_points = (unsigned int)(4 + (2 * (pow(2,ref) -1) ) + (2 * ( (side_ratio * pow(2,ref)) -1) ));

	// Find spacing of nodes
	spacing = (2 * (wid + len) ) / num_points;

	// Start locations of point generator
	double start_x = centre[0]-wid/2, start_y = centre[1]-len/2, x, y, z;
	
	// Number of points in each direction
	unsigned int np_x = (unsigned int)(len / spacing)+1, np_y = (unsigned int)(wid / spacing)+1;

	// Loop over matrix of points
	for (unsigned int i = 0; i < np_x; i++) {
		for (unsigned int j = 0; j < np_y; j++) {
			// 2D specification of z
			z = (b_z - a_z)/2;

			// x and y positions
			x = start_x + i*spacing;
			y = start_y + j*spacing;

			// Only add the marker if the point is on an edge
			if (	(i == 0 || i == np_x - 1) ||
					(j == 0 || j == np_y - 1)
				) {

				addMarker(x,y,z,flex_rigid);

			}				

		}
	}

#endif

	// Just in case anything goes wrong here...
	if (markers.size() != num_points) {
		std::cout << "Body is not closed. Exiting." << std::endl;
		system("pause");
		exit(EXIT_FAILURE);
	}

}

// ***************************************************************************************************
// Method to build a flexible filament
void IB_body::makeBody(std::vector<double> start_point, std::vector<double> end_point, std::vector<int> BCs, bool flex_rigid) {

	// Designate the body as being flexible or rigid
	this->flex_rigid = flex_rigid;

	// Assign delta_rho and flexural rigidity
	this->delta_rho = ibb_delta_rho;
	this->flexural_rigidity = ibb_EI;

	// Create zero tension vector for filament
	std::fill(tension.begin(), tension.end(), 0.0);

	// Find physical length of filament
	std::vector<double> L_vector;
	L_vector.push_back(end_point[0]-start_point[0]);
	L_vector.push_back(end_point[1]-start_point[1]);
	L_vector.push_back(end_point[2]-start_point[2]);
	double L = vecnorm(L_vector);

	// Compute angles of filament (in vertical and horizontal planes)
	double body_angle_h, body_angle_v;
	
	// Take into account if the filament is perpendicular to X axis
	if (L_vector[0] == 0) { body_angle_h = PI / 2; }
	else { body_angle_h = atan( L_vector[2] / L_vector[0] ); }

	// Take into account if the filament is perpendicular to Y axis
	if (vecnorm(L_vector[0],L_vector[2]) == 0) { body_angle_v = PI / 2; }
	else { body_angle_v = atan( L_vector[1] / vecnorm(L_vector[0],L_vector[2]) ); }

	// Compute spacing
	spacing = L / (num_markers - 1);				// Physical spacing between markers
	double spacing_h = spacing*cos(body_angle_v);	// Local spacing projected onto the horizontal plane

	// Add all markers
	for (size_t i = 0; i < num_markers; i++) {
		addMarker(	start_point[0] + i*spacing_h*cos(body_angle_h),
					start_point[1] + i*spacing*sin(body_angle_v),
					start_point[2] + i*spacing_h*sin(body_angle_h),
					true);
	}


	// Add dynamic positions in physical units
	for (size_t i = 0; i < num_markers; i++) {
		markers[i].position_old.push_back(markers[i].position[0]);
		markers[i].position_old.push_back(markers[i].position[1]);
		markers[i].position_old.push_back(markers[i].position[2]);
	}

	// Assign initial tension as zero
	tension.resize(markers.size());
	std::fill(tension.begin(), tension.end(), 0.0);

	// Correct for endpoint boundary conditions if necessary
	if (BCs[0] == 0) {
		markers[0].flex_rigid = false;
	}
	if (BCs[1] == 0) {
		markers[markers.size()-1].flex_rigid = false;
	}
	
}
// ***************************************************************************************************