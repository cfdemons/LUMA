/* This file defines the constructors and methods for the immersed boundary body object.
*/

#include "../inc/stdafx.h"
#include "../inc/IBBody.h"
#include "../inc/definitions.h"
#include <math.h>
#include "../inc/GridObj.h"

// ***************************************************************************************************
// Constructor and destructor
IBBody::IBBody()
{
	this->groupID = 0;	// Default ID
}
IBBody::~IBBody(void) { }

// ***************************************************************************************************
// Method to add marker
void IBBody::addMarker(double x, double y, double z, bool flex_rigid) {

	// Extend array of particles by 1 and construct a new IBMarker object
	markers.emplace_back( x, y, z, flex_rigid );

}

// ***************************************************************************************************
// Method to seed markers for a sphere / circle
void IBBody::makeBody(double radius, std::vector<double> centre,
					   bool flex_rigid, bool deform, unsigned int group) {

	// Designate body as beiong flexible or rigid and a closed surface
	this->flex_rigid = flex_rigid;
	this->closed_surface = true;

	// Designate deformable flag and groupID
	if (flex_rigid) {
		this->deformable = true;	// A flexible body must be deformable by definition
	} else {
		this->deformable = deform;
	}
	this->groupID = group;

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
	spacing = GridUtils::vecnorm( diff );



#else
	// Circle -- find theta
	std::vector<double> theta = GridUtils::linspace(0, 2*PI - (2*PI / num_markers), num_markers);
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
	spacing = GridUtils::vecnorm( diff );
#endif

}

// ***************************************************************************************************
// Method to seed markers for a cuboid/rectangle
void IBBody::makeBody(std::vector<double> width_length_depth, std::vector<double> angles, std::vector<double> centre,
					   bool flex_rigid, bool deform, unsigned int group) {

	// Designate body as being flexible or rigid and a closed surface
	this->flex_rigid = flex_rigid;
	this->closed_surface = true;

	// Designate deformable flag and groupID
	if (flex_rigid) {
		this->deformable = true;	// A flexible body must be deformable by definition
	} else {
		this->deformable = deform;
	}
	this->groupID = group;

	// Shorter variable names for convenience
	double wid = width_length_depth[0];
	double len = width_length_depth[1];

#if (dims == 3)

	// Cube //

	double dep = width_length_depth[2];

	// Check side lengths to make sure we can ensure points on the corners
	if (fmod(ibb_w,ibb_l) != 0 && fmod(len,wid) != 0 && fmod(len,ibb_d) != 0 && fmod(dep,len) != 0) {
			std::cout << "Error: See Log File" << std::endl;
			*GridUtils::logfile << "IB body cannot be built with uniform points. Change its dimensions. Exiting." << std::endl;
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
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "IB body does not have enough points. Need " << advisory_num_points << " to build body. Exiting." << std::endl;
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
	double start_x = centre[0]-wid/2, start_y = centre[1]-len/2, start_z = centre[2]-dep/2, x, y, z, xdash, ydash, zdash;

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

					// Transform x y and z based on rotation
					xdash = (x*cos(angles[0] * PI / 180) - y*sin(angles[0] * PI / 180))*cos(angles[1] * PI / 180)
						- z*sin(angles[1] * PI / 180);
					ydash = y*cos(angles[0] * PI / 180) + x*sin(angles[0] * PI / 180);
					zdash = (x*cos(angles[0] * PI / 180) - y*sin(angles[0] * PI / 180))*sin(angles[1] * PI / 180)
						+ z*cos(angles[1] * PI / 180);

					// Add marker
					addMarker(xdash,ydash,zdash,flex_rigid);

				}

			}
		}
	}


#else
	// Square //

	// Check side lengths to make sure we can ensure points on the corners
	if ((fmod(wid,len) != 0) && (fmod(len,wid) != 0)) {
			std::cout << "Error: See Log File" << std::endl;
			*GridUtils::logfile << "IB body cannot be built with uniform points. Change its dimensions. Exiting." << std::endl;
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
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "IB body does not have enough points. Need " << advisory_num_points << " to build body. Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}

	// Number of points required to get uniform distribution and points on corners
	unsigned int num_points = (unsigned int)(4 + (2 * (pow(2,ref) -1) ) + (2 * ( (side_ratio * pow(2,ref)) -1) ));

	// Find spacing of nodes
	spacing = (2 * (wid + len) ) / num_points;

	// Start locations of point generator
	double start_x = centre[0]-wid/2, start_y = centre[1]-len/2, x, y, z, xdash, ydash, zdash;

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

				// Transform x y and z based on rotation
				xdash = x*cos(angles[0] * PI / 180) - y*sin(angles[0] * PI / 180);
				ydash = x*sin(angles[0] * PI / 180) + y*cos(angles[0] * PI / 180);
				zdash = z;

				// Add marker
				addMarker(xdash,ydash,zdash,flex_rigid);

			}

		}
	}

#endif

	// Just in case anything goes wrong here...
	if (markers.size() != num_points) {
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "Body is not closed. Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}

}

// ***************************************************************************************************
// Method to seed markers for a flexible filament
void IBBody::makeBody(unsigned int nummarkers, std::vector<double> start_point, double fil_length, std::vector<double> angles,
					   std::vector<int> BCs, bool flex_rigid, bool deform, unsigned int group) {


	// **** Currently only allows start end to be simply supported or clamped and other end to be free ****
	if ( BCs[1] != 0  || BCs[0] == 0 ) {
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "Only allowed to have a fixed starting end and a free ending end of a filament at the minute. Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}

	// Designate BCs
	this->BCs = BCs;

	// Designate the body as being flexible or rigid and an open surface
	this->flex_rigid = flex_rigid;
	this->closed_surface = false;

	// Designate deformable flag and groupID
	if (flex_rigid) {
		this->deformable = true;	// A flexible body must be deformable by definition
	} else {
		this->deformable = deform;
	}
	this->groupID = group;

	// Assign delta_rho and flexural rigidity
	this->delta_rho = ibb_delta_rho;
	this->flexural_rigidity = ibb_EI;

	// Create zero tension vector for filament
	std::fill(tension.begin(), tension.end(), 0.0);

	// Get angles
	double body_angle_v = angles[0];
	double body_angle_h = angles[1];

	// Compute spacing
	spacing = fil_length / (nummarkers - 1);		// Physical spacing between markers
	double spacing_h = spacing*cos(body_angle_v * PI / 180);	// Local spacing projected onto the horizontal plane

	// Add all markers
	for (size_t i = 0; i < nummarkers; i++) {
		addMarker(	start_point[0] + i*spacing_h*cos(body_angle_h * PI / 180),
					start_point[1] + i*spacing*sin(body_angle_v * PI / 180),
					start_point[2] + i*spacing_h*sin(body_angle_h * PI / 180),
					true);
	}


	// Add dynamic positions in physical units
	for (size_t i = 0; i < nummarkers; i++) {
		markers[i].position_old.push_back(markers[i].position[0]);
		markers[i].position_old.push_back(markers[i].position[1]);
		markers[i].position_old.push_back(markers[i].position[2]);
	}

	// Assign initial tension as zero
	tension.resize(markers.size());
	std::fill(tension.begin(), tension.end(), 0.0);


}
// ***************************************************************************************************
// Method to seed markers for a 3D plate inclined from the xz plane
double IBBody::makeBody(std::vector<double> width_length, double angle, std::vector<double> centre,
	bool flex_rigid, bool deform, unsigned int group, bool plate) {

	// Exit if called in 2D
	if ( dims == 2 ) {
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "Plate builder must only be called in 3D. To build a 2D plate, use a rigid filament. Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}

	// Designate body as being flexible or rigid and an open surface
	this->flex_rigid = flex_rigid;
    this->closed_surface = false;

	// Designate deformable flag and groupID
	if (flex_rigid) {
		this->deformable = true;	// A flexible body must be deformable by definition
	} else {
		this->deformable = deform;
	}
	this->groupID = group;

	// Assign delta_rho and flexural rigidity
	this->delta_rho = ibb_delta_rho;
	this->flexural_rigidity = ibb_EI;

	// Shorter variable names for convenience
	double len_z, len_x;
	len_z = width_length[1];
	len_x = width_length[0];

	// Square //

	// Check side lengths to make sure we can ensure points on the corners
	if ((fmod(len_z,len_x) != 0) && (fmod(len_x,len_z) != 0)) {
			std::cout << "Error: See Log File" << std::endl;
			*GridUtils::logfile << "IB body cannot be built with uniform points. Change its dimensions. Exiting." << std::endl;
			exit(EXIT_FAILURE);
		}

	// Get ratio of sides and degree of point refinement
	int side_ratio;
	if (len_z > len_x) { // x length is shortest
		side_ratio = (int)(len_z / len_x);
	} else { // z length is shortest
		side_ratio = (int)(len_x / len_z);
	}

	// The following is derived by inspection of the problem and ensures a point on each corner of the body
	// Double the number of points to be used as can be rounded down to only a relatively small number
	int ref = (int)ceil( log( 2*num_markers / (2 * (1+side_ratio)) ) / log(2) );

	// Check to see if enough points
	if (ref == 0) {
		// Advisory of number of points
		unsigned int advisory_num_points = (unsigned int)(4 + (2 * (pow(2,1) -1) ) + (2 * ( (side_ratio * pow(2,1)) -1) ) );
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "IB body does not have enough points. Need " << advisory_num_points << " to build body. Exiting." << std::endl;
		exit(EXIT_FAILURE);
	}

	// Number of points required to get uniform distribution and points on corners
	unsigned int num_points = (unsigned int)(4 + (2 * (pow(2,ref) -1) ) + (2 * ( (side_ratio * pow(2,ref)) -1) ));

	// Find spacing of nodes
	spacing = (2 * (len_z + len_x) ) / num_points;

	// Start locations of point generator (bottom corner)
	double start_z = centre[2]-len_z/2,
		start_x = centre[0] - (len_x/2)*cos(angle * PI / 180),
		start_y = centre[1] - (len_x/2)*sin(angle * PI / 180),
		x, y, z;

	// Number of points in each direction
	unsigned int np_x = (unsigned int)(len_x / spacing)+1, np_z = (unsigned int)(len_z / spacing)+1;

	// Loop over matrix of points
	for (unsigned int i = 0; i < np_x; i++) {
		for (unsigned int j = 0; j < np_z; j++) {

			// x and y positions
			x = start_x + (i*spacing*cos(angle * PI / 180));
			y = start_y + (i*spacing*sin(angle * PI / 180));
			z = start_z + (j*spacing);

			// Add marker
			addMarker(x,y,z,flex_rigid);


		}
	}

	// Add dynamic positions in physical units
	for (size_t i = 0; i < np_x*np_z; i++) {
		markers[i].position_old.push_back(markers[i].position[0]);
		markers[i].position_old.push_back(markers[i].position[1]);
		markers[i].position_old.push_back(markers[i].position[2]);
	}

	// Assign initial tension as zero
	tension.resize(markers.size());
	std::fill(tension.begin(), tension.end(), 0.0);

	return spacing;


}

// ***************************************************************************************************
