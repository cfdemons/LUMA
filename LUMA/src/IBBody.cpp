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

/* This file defines the constructors and methods for the immersed boundary body object.
*/

#include "../inc/stdafx.h"
#include "../inc/IBBody.h"
#include "../inc/IBMarker.h"
#include "../inc/PCpts.h"
#include "../inc/ObjectManager.h"

// ***************************************************************************************************
/// \brief Constructor which sets group ID to zero by default.
IBBody::IBBody()
{
	this->_Owner = nullptr;
	this->id = 0;
}

/// \brief	Constructor which assigns the owner grid.
///
///			Also sets the group ID to zero.
///
/// \param g	pointer to owner grid
/// \param id	ID of body in array of bodies.
IBBody::IBBody(GridObj* g, size_t id)
{
	this->_Owner = g;
	this->id = id;
};

// ***************************************************************************************************
/// \brief	Constructor to build a body in place using a point cloud,
///
///			isFlexible and isMovable properties taken from definitions.
///
/// \param g		pointer to owner grid
/// \param id		ID of body in array of bodies.
/// \param _PCpts	pointer to point cloud data.
IBBody::IBBody(GridObj* g, size_t id, PCpts* _PCpts, eMoveableType moveProperty, bool clamped) : Body(g, id, _PCpts)
{

	// Set movable and flexible parameters
	if (moveProperty == eFlexible) {
		this->isFlexible = true;
		this->isMovable = true;
	}
	else if (moveProperty == eMovable) {
		this->isFlexible = false;
		this->isMovable = true;
	}
	else if (moveProperty == eRigid) {
		this->isFlexible = false;
		this->isMovable = false;
	}

	// Delete any markers which are on the receiver layer as not needed for IBM
	*GridUtils::logfile << "Deleting IB markers which exist on receiver layer..." << std::endl;
	deleteRecvLayerMarkers();
}


// ***************************************************************************************************
/// \brief	Constructor to build a body in place using a point cloud,
///
///			isFlexible and isMovable properties taken from definitions.
///
/// \param g		pointer to owner grid
/// \param id		ID of body in array of bodies.
/// \param _PCpts	pointer to point cloud data.
IBBody::IBBody(GridObj* g, size_t bodyID, int lev, int reg, std::vector<double> &start_position,
		double length, std::vector<double> &angles, eMoveableType moveProperty, bool clamped) : Body(g, bodyID, lev, reg, start_position, length, angles)
{

	// Set movable and flexible parameters
	if (moveProperty == eFlexible) {
		this->isFlexible = true;
		this->isMovable = true;
	}
	else if (moveProperty == eMovable) {
		this->isFlexible = false;
		this->isMovable = true;
	}
	else if (moveProperty == eRigid) {
		this->isFlexible = false;
		this->isMovable = false;
	}

	// Delete any markers which are on the receiver layer as not needed for IBM
	*GridUtils::logfile << "Deleting IB markers which exist on receiver layer..." << std::endl;
	deleteRecvLayerMarkers();

}

/// Default destructor
IBBody::~IBBody(void)
{
}


/*************** PREFAB MAKE BODY METHODS ********************/

// ***************************************************************************************************
/// \brief	Method to seed markers for a sphere / circle.
/// \param radius radius of circle/sphere.
/// \param centre position vector of circle/sphere centre.
/// \param isFlexible flag to indicate whether body is flexible and requires a structural calculation.
/// \param isMovable flag to indicate whether body is movable and requires relocation each time step.
/// \param group ID indicating which group the body is part of for collective operations.
void IBBody::makeBody(double radius, std::vector<double> centre,
	bool isFlexible, bool isMovable, int group) {

//	// Designate body as beiong flexible or rigid and a closed surface
//	this->isFlexible = isFlexible;
//	this->closed_surface = true;
//
//	// Designate isMovable flag and groupID
//	if (isFlexible) {
//		this->isMovable = true;	// A flexible body must be movable by definition
//	} else {
//		this->isMovable = isMovable;
//	}
//	this->groupID = group;
//
//#if (L_DIMS == 3)
//
//	// Sphere //
//
//	// Following code for point generation on unit sphere actually seeds
//	// using Fibonacci sphere technique. Code is not my own but works.
//	double inc = L_PI * (3 - sqrt(5));
//    double off = 2.0 / (float)L_NUM_MARKERS ;
//    for (int k = 0; k < L_NUM_MARKERS; k++) {
//        double y = k * off - 1 + (off / 2);
//        double r = sqrt(1 - y*y);
//        double phi = k * inc;
//
//		// Add Lagrange marker to body (scale by radius)
//        addMarker(centre[0] + (cos(phi)*r * radius), y*radius + centre[1], centre[2] + (sin(phi)*r*radius), isFlexible);
//	}
//
//	// Spacing (assuming all Lagrange markers are uniformly spaced)
//	std::vector<double> diff;
//	for (int d = 0; d < L_DIMS; d++) {
//		diff.push_back ( markers[1].position[d] - markers[0].position[d] );
//	}
//	spacing = GridUtils::vecnorm( diff );
//
//
//
//#else
//	// Circle -- find theta
//	std::vector<double> theta = GridUtils::linspace(0, 2*L_PI - (2*L_PI / L_NUM_MARKERS), L_NUM_MARKERS);
//	for (size_t i = 0; i < theta.size(); i++) {
//
//		// Add Lagrange marker to body
//		addMarker(	centre[0] + radius * cos(theta[i]),
//					centre[1] + radius * sin(theta[i]),
//					0.0, isFlexible	);
//	}
//
//	// Spacing
//	std::vector<double> diff;
//	for (int d = 0; d < L_DIMS; d++) {
//		diff.push_back ( markers[1].position[d] - markers[0].position[d] );
//	}
//	spacing = GridUtils::vecnorm( diff );
//#endif

}

// ***************************************************************************************************
/// \brief	Method to seed markers for a cuboid / rectangle.
/// \param width_length_depth principal dimensions of cuboid / rectangle.
/// \param angles principal orientation of cuboid / rectangle w.r.t. domain axes.
/// \param centre position vector of cuboid / rectangle centre.
/// \param isFlexible flag to indicate whether body is flexible and requires a structural calculation.
/// \param isMovable flag to indicate whether body is movable and requires relocation each time step.
/// \param group ID indicating which group the body is part of for collective operations.
void IBBody::makeBody(std::vector<double> width_length_depth, std::vector<double> angles, std::vector<double> centre,
					   bool isFlexible, bool isMovable, int group) {

//	// Designate body as being flexible or rigid and a closed surface
//	this->isFlexible = isFlexible;
//	this->closed_surface = true;
//
//	// Designate isMovable flag and groupID
//	if (isFlexible) {
//		this->isMovable = true;	// A flexible body must be movable by definition
//	} else {
//		this->isMovable = isFlexible;
//	}
//	this->groupID = group;
//
//	// Shorter variable names for convenience
//	double wid = width_length_depth[0];
//	double len = width_length_depth[1];
//
//#if (L_DIMS == 3)
//
//	// Cube //
//
//	double dep = width_length_depth[2];
//
//	// Check side lengths to make sure we can ensure points on the corners
//	if (fmod(L_IBB_W,L_IBB_L) != 0 && fmod(len,wid) != 0 && fmod(len,L_IBB_D) != 0 && fmod(dep,len) != 0) {
//		L_ERROR("IB body cannot be built with uniform points. Change its dimensions. Exiting.",
//			GridUtils::logfile);
//		}
//
//	// Get ratio of sides and degree of point refinement
//	int side_ratio_1, side_ratio_2;
//	double smallest_side;
//	if (wid > len) {
//		if (dep > len) { // Length is smallest
//			smallest_side = len;
//			side_ratio_1 = (int)(wid / len);
//			side_ratio_2 = (int)(dep / len);
//		} else { // Depth is smallest
//			smallest_side = dep;
//			side_ratio_1 = (int)(wid / dep);
//			side_ratio_2 = (int)(len / dep);
//		}
//	} else if (dep > wid) { // Width is smallest
//		smallest_side = wid;
//		side_ratio_1 = (int)(len / wid);
//		side_ratio_2 = (int)(dep / wid);
//	} else { // Depth is smallest
//		smallest_side = dep;
//		side_ratio_1 = (int)(wid / dep);
//		side_ratio_2 = (int)(len / dep);
//	}
//
//	// The following is derived by inspection of the problem and ensures a point on each corner of the body
//	// Since we have points on the faces as well as the edges we have a quadratic expression to be solved to
//	// find the level of refinement required to get the nearest number of points to L_NUM_MARKERS.
//	int ref;
//
//	int a = 2*side_ratio_1 + 2*side_ratio_1 + 2*side_ratio_1*side_ratio_2;
//	int b = 4 + 4*side_ratio_1 + 4*side_ratio_2;
//	int c = 8-L_NUM_MARKERS;
//
//	double P = (-b + sqrt(pow(b,2) - (4*a*c)) )/(2*a) + 1;
//
//	ref = (int)ceil( log( P ) / log(2) );
//
//
//	// Check to see if enough points
//	if (ref == 0) {
//		// Advisory of number of points
//		int advisory_num_points = (int)(8 +
//			(4 * (pow(2,1) -1) ) +
//			(4 * ( (side_ratio_1 * pow(2,1)) -1) ) +
//			(4 * ( (side_ratio_2 * pow(2,1)) -1) ) +
//			(2 * ( (pow(2,1) -1)*side_ratio_1 * (pow(2,1) -1) )) +
//			(2 * ( (pow(2,1) -1)*side_ratio_2 * (pow(2,1) -1) )) +
//			(2 * ( (pow(2,1) -1)*side_ratio_1 * (pow(2,1) -1)*side_ratio_2 ))
//		);
//		L_ERROR("IB body does not have enough points. Need " + std::to_string(advisory_num_points) + " to build body. Exiting.",
//			GridUtils::logfile);
//	}
//
//	// Number of points required to get uniform distribution and points on corners
//	int num_points = (int)(8 +
//			(4 * (pow(2,ref) -1) ) +
//			(4 * ( (side_ratio_1 * pow(2,ref)) -1) ) +
//			(4 * ( (side_ratio_2 * pow(2,ref)) -1) ) +
//			(2 * ( (pow(2,ref) -1)*side_ratio_1 * (pow(2,ref) -1) )) +
//			(2 * ( (pow(2,ref) -1)*side_ratio_2 * (pow(2,ref) -1) )) +
//			(2 * ( (pow(2,ref) -1)*side_ratio_1 * (pow(2,ref) -1)*side_ratio_2 ))
//		);
//
//	// Spacing from smallest side
//	spacing = smallest_side / pow(2,ref);
//
//	// Start locations of point generator
//	double start_x = centre[0]-wid/2, start_y = centre[1]-len/2, start_z = centre[2]-dep/2, x, y, z, xdash, ydash, zdash;
//
//	// Number of points in each direction
//	int np_x = (int)(len / spacing)+1, np_y = (int)(wid / spacing)+1, np_z = (int)(dep / spacing)+1;
//
//	// Loop over matrix of points
//	for (int i = 0; i < np_x; i++) {
//		for (int j = 0; j < np_y; j++) {
//			for (int k = 0; k < np_z; k++) {
//
//				// x, y and z positions
//				x = start_x + i*spacing;
//				y = start_y + j*spacing;
//				z = start_z + k*spacing;
//
//				// Only add the marker if the point is on an edge
//				if (	(i == 0 || i == np_x - 1) ||
//						(j == 0 || j == np_y - 1) ||
//						(k == 0 || k == np_z - 1)
//					) {
//
//					// Transform x y and z based on rotation
//					xdash = (x*cos(angles[0] * L_PI / 180) - y*sin(angles[0] * L_PI / 180))*cos(angles[1] * L_PI / 180)
//						- z*sin(angles[1] * L_PI / 180);
//					ydash = y*cos(angles[0] * L_PI / 180) + x*sin(angles[0] * L_PI / 180);
//					zdash = (x*cos(angles[0] * L_PI / 180) - y*sin(angles[0] * L_PI / 180))*sin(angles[1] * L_PI / 180)
//						+ z*cos(angles[1] * L_PI / 180);
//
//					// Add marker
//					addMarker(xdash,ydash,zdash,isFlexible);
//
//				}
//
//			}
//		}
//	}
//
//
//#else
//	// Square //
//
//	// Check side lengths to make sure we can ensure points on the corners
//	if ((fmod(wid,len) != 0) && (fmod(len,wid) != 0)) {
//		L_ERROR("IB body cannot be built with uniform points. Change its dimensions. Exiting.",
//			GridUtils::logfile);
//		}
//
//	// Get ratio of sides and degree of point refinement
//	int side_ratio;
//	if (wid > len) { // Length is shortest
//		side_ratio = (int)(wid / len);
//	} else { // Width is shortest
//		side_ratio = (int)(len / wid);
//	}
//
//	// The following is derived by inspection of the problem and ensures a point on each corner of the body
//	// Double the number of points to be used as can be rounded down to only a relatively small number
//	int ref = (int)ceil( log( L_NUM_MARKERS / (2 * (1+side_ratio)) ) / log(2) );
//
//	// Check to see if enough points
//	if (ref == 0) {
//		// Advisory of number of points
//		int advisory_num_points = (int)(4 + (2 * (pow(2,1) -1) ) + (2 * ( (side_ratio * pow(2,1)) -1) ) );
//		L_ERROR("IB body does not have enough points. Need " + std::to_string(advisory_num_points) + " to build body. Exiting.",
//			GridUtils::logfile);
//	}
//
//	// Number of points required to get uniform distribution and points on corners
//	int num_points = (int)(4 + (2 * (pow(2,ref) -1) ) + (2 * ( (side_ratio * pow(2,ref)) -1) ));
//
//	// Find spacing of nodes
//	spacing = (2 * (wid + len) ) / num_points;
//
//	// Start locations of point generator
//	double start_x = centre[0]-wid/2, start_y = centre[1]-len/2, x, y, z, xdash, ydash, zdash;
//
//	// Number of points in each direction
//	int np_x = (int)(len / spacing)+1, np_y = (int)(wid / spacing)+1;
//
//	// Loop over matrix of points
//	for (int i = 0; i < np_x; i++) {
//		for (int j = 0; j < np_y; j++) {
//			// 2D specification of z
//			z = L_BZ/2;
//
//			// x and y positions
//			x = start_x + i*spacing;
//			y = start_y + j*spacing;
//
//			// Only add the marker if the point is on an edge
//			if (	(i == 0 || i == np_x - 1) ||
//					(j == 0 || j == np_y - 1)
//				) {
//
//				// Transform x y and z based on rotation
//				xdash = x*cos(angles[0] * L_PI / 180) - y*sin(angles[0] * L_PI / 180);
//				ydash = x*sin(angles[0] * L_PI / 180) + y*cos(angles[0] * L_PI / 180);
//				zdash = z;
//
//				// Add marker
//				addMarker(xdash,ydash,zdash,isFlexible);
//
//			}
//
//		}
//	}
//
//#endif
//
//	// Just in case anything goes wrong here...
//	if (markers.size() != num_points) {
//		L_ERROR("Body is not closed. Exiting.", GridUtils::logfile);
//	}

}

// ***************************************************************************************************
/// \brief	Method to seed markers for a flexible filament.
/// \param nummarkers number of markers to use for filament.
/// \param start_point 3D position vector of the start of the filament.
/// \param fil_length length of filament in physical units.
/// \param angles two angles representing filament inclination w.r.t. domain axes 
///			(horizontal plane and vertical plane).
/// \param BCs vector containing start and end boundary condition types (see class definition for valid values).
/// \param isFlexible flag to indicate whether body is flexible and requires a structural calculation.
/// \param isMovable flag to indicate whether body is movable and requires relocation each time step.
/// \param group ID indicating which group the body is part of for collective operations.
void IBBody::makeBody(int nummarkers, std::vector<double> start_point, double fil_length, std::vector<double> angles,
	std::vector<int> BCs, bool isFlexible, bool isMovable, int group) {


//	// **** Currently only allows start end to be simply supported or clamped and other end to be free ****
//	if ( BCs[1] != 0  || BCs[0] == 0 ) {
//		L_ERROR("Only allowed to have a fixed starting end and a free ending end of a filament at the minute. Exiting.",
//			GridUtils::logfile);
//	}
//
//	// Designate BCs
//	this->BCs = BCs;
//
//	// Designate the body as being flexible or rigid and an open surface
//	this->isFlexible = isFlexible;
//	this->closed_surface = false;
//
//	// Designate isMovable flag and groupID
//	if (isFlexible) {
//		this->isMovable = true;	// A flexible body must be movable by definition
//	} else {
//		this->isMovable = isFlexible;
//	}
//	this->groupID = group;
//
//	// Assign delta_rho and flexural rigidity
//	this->delta_rho = L_IBB_DELTA_RHO;
//	this->flexural_rigidity = L_IBB_EI;
//
//	// Create zero tension vector for filament
//	std::fill(tension.begin(), tension.end(), 0.0);
//
//	// Get angles
//	double body_angle_v = angles[0];
//	double body_angle_h = angles[1];
//
//	// Compute spacing
//	spacing = fil_length / (nummarkers - 1);		// Physical spacing between markers
//	double spacing_h = spacing*cos(body_angle_v * L_PI / 180);	// Local spacing projected onto the horizontal plane
//
//	// Add all markers
//	for (int i = 0; i < nummarkers; i++) {
//		addMarker(	start_point[0] + i*spacing_h*cos(body_angle_h * L_PI / 180),
//					start_point[1] + i*spacing*sin(body_angle_v * L_PI / 180),
//					start_point[2] + i*spacing_h*sin(body_angle_h * L_PI / 180),
//					true);
//	}
//
//
//	// Add dynamic positions in physical units
//	for (int i = 0; i < nummarkers; i++) {
//		markers[i].position_old.push_back(markers[i].position[0]);
//		markers[i].position_old.push_back(markers[i].position[1]);
//		markers[i].position_old.push_back(markers[i].position[2]);
//	}
//
//	// Assign initial tension as zero
//	tension.resize(markers.size());
//	std::fill(tension.begin(), tension.end(), 0.0);


}

// ***************************************************************************************************
/// \brief	Method to seed markers for a 3D plate inclined from the XZ plane.
/// \param width_length 2D vector of principal dimensions of thin plate.
/// \param angle inclination angle from horizontal.
/// \param centre position vector of the plate centre.
/// \param isFlexible flag to indicate whether body is flexible and requires a structural calculation.
/// \param isMovable flag to indicate whether body is movable and requires relocation each time step.
/// \param group ID indicating which group the body is part of for collective operations.
/// \param plate arbitrary argument to allow overload otherwise would have the same signature as a filament builder.
double IBBody::makeBody(std::vector<double> width_length, double angle, std::vector<double> centre,
	bool isFlexible, bool isMovable, int group, bool plate) {

//	// Exit if called in 2D
//	if ( L_DIMS == 2 ) {
//		L_ERROR("Plate builder must only be called in 3D. To build a 2D plate, use a rigid filament. Exiting.",
//			GridUtils::logfile);
//	}
//
//	// Designate body as being flexible or rigid and an open surface
//	this->isFlexible = isFlexible;
//    this->closed_surface = false;
//
//	// Designate isMovable flag and groupID
//	if (isFlexible) {
//		this->isMovable = true;	// A flexible body must be movable by definition
//	} else {
//		this->isMovable = isFlexible;
//	}
//	this->groupID = group;
//
//	// Assign delta_rho and flexural rigidity
//	this->delta_rho = L_IBB_DELTA_RHO;
//	this->flexural_rigidity = L_IBB_EI;
//
//	// Shorter variable names for convenience
//	double len_z, len_x;
//	len_z = width_length[1];
//	len_x = width_length[0];
//
//	// Square //
//
//	// Check side lengths to make sure we can ensure points on the corners
//	if ((fmod(len_z,len_x) != 0) && (fmod(len_x,len_z) != 0)) {
//		L_ERROR("IB body cannot be built with uniform points. Change its dimensions. Exiting.",
//			GridUtils::logfile);
//		}
//
//	// Get ratio of sides and degree of point refinement
//	int side_ratio;
//	if (len_z > len_x) { // x length is shortest
//		side_ratio = (int)(len_z / len_x);
//	} else { // z length is shortest
//		side_ratio = (int)(len_x / len_z);
//	}
//
//	// The following is derived by inspection of the problem and ensures a point on each corner of the body
//	// Double the number of points to be used as can be rounded down to only a relatively small number
//	int ref = (int)ceil( log( 2*L_NUM_MARKERS / (2 * (1+side_ratio)) ) / log(2) );
//
//	// Check to see if enough points
//	if (ref == 0) {
//		// Advisory of number of points
//		int advisory_num_points = (int)(4 + (2 * (pow(2,1) -1) ) + (2 * ( (side_ratio * pow(2,1)) -1) ) );
//		L_ERROR("IB body does not have enough points. Need " + std::to_string(advisory_num_points) + " to build body. Exiting.",
//			GridUtils::logfile);
//	}
//
//	// Number of points required to get uniform distribution and points on corners
//	int num_points = (int)(4 + (2 * (pow(2,ref) -1) ) + (2 * ( (side_ratio * pow(2,ref)) -1) ));
//
//	// Find spacing of nodes
//	spacing = (2 * (len_z + len_x) ) / num_points;
//
//	// Start locations of point generator (bottom corner)
//	double start_z = centre[2]-len_z/2,
//		start_x = centre[0] - (len_x/2)*cos(angle * L_PI / 180),
//		start_y = centre[1] - (len_x/2)*sin(angle * L_PI / 180),
//		x, y, z;
//
//	// Number of points in each direction
//	int np_x = (int)(len_x / spacing)+1, np_z = (int)(len_z / spacing)+1;
//
//	// Loop over matrix of points
//	for (int i = 0; i < np_x; i++) {
//		for (int j = 0; j < np_z; j++) {
//
//			// x and y positions
//			x = start_x + (i*spacing*cos(angle * L_PI / 180));
//			y = start_y + (i*spacing*sin(angle * L_PI / 180));
//			z = start_z + (j*spacing);
//
//			// Add marker
//			addMarker(x,y,z,isFlexible);
//
//
//		}
//	}
//
//	// Add dynamic positions in physical units
//	for (int i = 0; i < np_x*np_z; i++) {
//		markers[i].position_old.push_back(markers[i].position[0]);
//		markers[i].position_old.push_back(markers[i].position[1]);
//		markers[i].position_old.push_back(markers[i].position[2]);
//	}
//
//	// Assign initial tension as zero
//	tension.resize(markers.size());
//	std::fill(tension.begin(), tension.end(), 0.0);
//
//	return spacing;


}
