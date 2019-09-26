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

#ifndef BODY_H
#define BODY_H

#include "stdafx.h"
class GridObj;
#include "PCpts.h"
#include "GridUtils.h"
#include "MarkerData.h"


/// \brief	Generic body class
///
///			Can consist of any type of Marker so templated.
template <typename MarkerType>
class Body
{

	// Make MPIManager a friend so it can access body data
	friend class MpiManager;

public:

	Body(void);					// Default Constructor
	virtual ~Body(void);		// Default destructor

	// Custom constructor which takes pointer to point cloud data and a pointer to the grid owner for the labelling
	Body(GridObj* g, int bodyID, PCpts* _PCpts);

	// Custom constructor for building prefab circle or sphere
	Body(GridObj* g, int bodyID, std::vector<double> &centre_point, double radius);

	// Custom constructor for building prefab square or cuboid
	Body(GridObj* g, int bodyID, std::vector<double> &centre_point, std::vector<double> &dimensions, std::vector<double> &angles);

	// Custom constructor for building prefab square or cuboid
	Body(GridObj* g, int bodyID, std::vector<double> &centre_point, double length, double width, std::vector<double> &angles);

	// Custom constructor for building prefab filament
	Body(GridObj* g, int bodyID, std::vector<double> &start_position, double length, std::vector<double> &angles);

	// ************************ Members ************************ //

protected:
	GridObj* _Owner;					///< Pointer to owning grid
	int id;								///< Unique ID of the body
	bool closed_surface;				///< Flag to specify whether or not it is a closed surface (i.e. last marker should link to first)
	int owningRank;						///< ID of the rank that owns this body (for epsilon and structural calculation)
	std::vector<MarkerType> markers;	///< Array of markers which make up the body
	int level;							///< Level on which body exists

	std::vector<int> validMarkers;		///< Vector of indices to valid markers within this body which actually exist on this rank


	// ************************ Methods ************************ //

	virtual void addMarker(double x, double y, double z, int markerID);		// Add a marker (can be overrriden)
	MarkerData* getMarkerData(double x, double y, double z);				// Retireve nearest marker data
	void passToVoxelFilter(double x, double y, double z, int markerID,
		int& curr_mark, std::vector<int>& counter);							// Voxelising marker adder
	void deleteRecvLayerMarkers();											// Delete any markers which are on receiver layer
	void deleteOffRankMarkers();											// Delete any markers which don't exist on this rank

private:
	bool isInVoxel(double x, double y, double z, int curr_mark);			// Check a point is inside an existing marker voxel
	bool isVoxelMarkerVoxel(double x, double y, double z);					// Check whether nearest voxel is a marker voxel
	int assignOwningRank(int id);											// Assign owning rank based on which ranks own which grids


protected:
	void buildFromCloud(PCpts *_PCpts);										// Method to build body from point cloud
	virtual void writeVtkPosition(int tval);								// VTK body writer
};


// ************************ Implementation of Body methods ************************ //

/// Default Constructor
template <typename MarkerType>
Body<MarkerType>::Body(void) 
	: _Owner(nullptr)
{
};

/// Default destructor
template <typename MarkerType>
Body<MarkerType>::~Body(void)
{
};

/*********************************************/
/// \brief	Custom constructor to populate body from array of points
///
///	\param g		hierarchy pointer to grid hierarchy
/// \param bodyID	ID of body in array of bodies
/// \param _PCpts	pointer to point cloud data
template <typename MarkerType>
Body<MarkerType>::Body(GridObj* g, int bodyID, PCpts* _PCpts)
	: _Owner(g), id(bodyID)
{
	// Set as unclosed surface by default
	this->closed_surface = false;

	// Set level
	this->level = _Owner->level;

	// Set the rank which owns this body
	this->owningRank = assignOwningRank(id);

	// Call method to build from point cloud
	this->buildFromCloud(_PCpts);
};

/*********************************************/
/// \brief	Custom constructor for building prefab filament
///
/// \param 	g				hierarchy pointer to grid hierarchy
/// \param 	bodyID			ID of body in array of bodies
/// \param 	start_position	start position of base of filament
/// \param 	length			length of filament
/// \param 	angles			angle of filament
template <typename MarkerType>
Body<MarkerType>::Body(GridObj* g, int bodyID, std::vector<double> &start_position, double length, std::vector<double> &angles)
{

	// Set the body base class parameters from constructor inputs
	this->_Owner = g;
	this->id = bodyID;
	this->closed_surface = false;

	// Set level
	this->level = _Owner->level;

	// Set the rank which owns this body
	this->owningRank = assignOwningRank(id);

	// Get horizontal and vertical angles
	double body_angle_v = angles[0];
#if (L_DIM == 3)
	double body_angle_h = angles[1];
#else
	double body_angle_h = 0.0;
#endif

	// Compute spacing
	int numMarkers = static_cast<int>(std::floor(length / g->dh)) + 1;
	double spacing = length / (numMarkers - 1);							// Physical spacing between markers
	double spacing_h = spacing * cos(body_angle_v * L_PI / 180);		// Local spacing projected onto the horizontal plane

	// Add all markers
	for (int i = 0; i < numMarkers; i++) {
		addMarker(	start_position[0] + i * spacing_h * cos(body_angle_h * L_PI / 180.0),
					start_position[1] + i * spacing * sin(body_angle_v * L_PI / 180.0),
					start_position[2] + i * spacing_h * sin(body_angle_h * L_PI / 180.0),
					i);
	}

	// Get rank
	int rank = GridUtils::safeGetRank();

	// Delete markers which exist off rank
	if (rank != owningRank) {
		*GridUtils::logfile << "Deleting markers which are not on this rank..." << std::endl;
		deleteOffRankMarkers();
	}
};


/*********************************************/
/// \brief	Custom constructor for building prefab circle/sphere
///
/// \param 	g				hierarchy pointer to grid hierarchy
/// \param 	bodyID			ID of body in array of bodies
/// \param 	centre			centre point of circle
/// \param 	radius			radius of circle
template <typename MarkerType>
Body<MarkerType>::Body(GridObj* g, int bodyID, std::vector<double> &centre, double radius)
{
	// Set the body base class parameters from constructor inputs
	this->_Owner = g;
	this->id = bodyID;
	this->closed_surface = true;

	// Set level
	this->level = _Owner->level;

	// Set the rank which owns this body
	this->owningRank = assignOwningRank(id);


	// Build sphere (3D)
#if (L_DIMS == 3)	// TODO Sort out 3D sphere builder

	// Sphere //

	// Following code for point generation on unit sphere actually seeds
	// using Fibonacci sphere technique. Code is not my own but works.
	int numMarkers = static_cast<int>(std::floor(4.0 * L_PI * pow(radius, 2.0) / pow(g->dh, 2.0)));
	double inc = L_PI * (3 - sqrt(5));
	double off = 2.0 / (float)numMarkers ;
	for (int k = 0; k < numMarkers; k++) {
		double y = k * off - 1 + (off / 2);
		double r = sqrt(1 - y*y);
		double phi = k * inc;

		// Add Lagrange marker to body (scale by radius)
		addMarker(centre[0] + (cos(phi)*r * radius), y*radius + centre[1], centre[2] + (sin(phi)*r*radius), k);
	}
#else
	
	// Build circle (2D)
	int numMarkers = static_cast<int>(std::floor(2.0 * L_PI * radius / g->dh));
	std::vector<double> theta = GridUtils::linspace(0.0, 2.0 * L_PI - (2.0 * L_PI / static_cast<double>(numMarkers)), numMarkers);
	
	for (size_t i = 1; i < theta.size(); i++) {

		// Add Lagrange marker to body
		addMarker(	centre[0] + radius * cos(theta[i]),
					centre[1] + radius * sin(theta[i]),
					centre[2],
					static_cast<int>(i)-1 );
	}
#endif

	
	// Get rank
	int rank = GridUtils::safeGetRank();

	// Delete markers which exist off rank
	if (rank != owningRank) {
		*GridUtils::logfile << "Deleting markers which are not on this rank..." << std::endl;
		deleteOffRankMarkers();
	}
};


/*********************************************/
/// \brief	Custom constructor for building square/cuboid
///
/// \param 	g					hierarchy pointer to grid hierarchy
/// \param 	bodyID				ID of body in array of bodies
/// \param 	centre				centre point of square
/// \param 	width_length_depth	dimensions of square
/// \param 	angles				angle of square
template <typename MarkerType>
Body<MarkerType>::Body(GridObj* g, int bodyID, std::vector<double> &centre,	
	std::vector<double> &width_length_depth, std::vector<double> &angles)
{

	// Set the body base class parameters from constructor inputs
	this->_Owner = g;
	this->id = bodyID;
	this->closed_surface = true;

	// Set level
	this->level = _Owner->level;

	// Set the rank which owns this body
	this->owningRank = assignOwningRank(id);


	// Shorter variable names for convenience
	double length = width_length_depth[0];
	double height = width_length_depth[1];
	double depth = width_length_depth[2];

	// Get number of markers in each direction
	int numMarkersLength = static_cast<int>(std::floor(length / g->dh)) + 1;
	int numMarkersHeight = static_cast<int>(std::floor(height / g->dh)) + 1;
	int numMarkersDepth = static_cast<int>(std::floor(depth / g->dh)) + 1;

	// Get total number of markers and check
#if (L_DIMS == 3)
	int numMarkers = numMarkersLength * numMarkersHeight * numMarkersDepth - 
		((numMarkersLength - 2) * (numMarkersHeight - 2) * (numMarkersDepth - 2));

	// Check side lengths to make sure we can ensure points on the corners
	if (numMarkers < 24) {
		L_ERROR("Resolution not high enough to build cuboid this small. Change its dimensions. Exiting.",
		GridUtils::logfile);
	}
#else
	int numMarkers = numMarkersLength * numMarkersHeight - ((numMarkersLength - 2) * (numMarkersHeight - 2));

	// Check side lengths to make sure we can ensure points on the corners
	if (numMarkers < 8) {
		L_ERROR("Resolution not high enough to build square this small. Change its dimensions. Exiting.",
		GridUtils::logfile);
	}
#endif

	// Start locations of point generator
	double x, y, z, xdash, ydash, zdash;

	// Get spacing between markers
	double spacingLength = length / (numMarkersLength - 1);
	double spacingHeight = height / (numMarkersHeight - 1);
	double spacingDepth = depth / (numMarkersDepth - 1);

	// Marker ID
	int markerID = 0;

#if (L_DIMS == 3)

	// Build cuboid
	for (int i = 0; i < numMarkersLength; i++) {
		for (int j = 0; j < numMarkersHeight; j++) {
			for (int k = 0; k < numMarkersDepth; k++) {


				// x and y position
				x = -length / 2.0 + i * spacingLength;
				y = -height / 2.0 + j * spacingHeight;
				z = -depth / 2.0 + k * spacingDepth;

				// Only add the marker if the point is on an edge
				if ((i == 0 || i == numMarkersLength - 1) || (j == 0 || j == numMarkersHeight - 1) || (k == 0 || k == numMarkersDepth - 1)) {

					// Transform x y and z based on rotation
					xdash = (x * cos(angles[0] * L_PI / 180) - y * sin(angles[0] * L_PI / 180)) * cos(angles[1] * L_PI / 180) - z * sin(angles[1] * L_PI / 180);
					ydash = (y * cos(angles[0] * L_PI / 180) + x * sin(angles[0] * L_PI / 180));
					zdash = (x * cos(angles[0] * L_PI / 180) - y * sin(angles[0] * L_PI / 180)) * sin(angles[1] * L_PI / 180) + z * cos(angles[1] * L_PI / 180);
					xdash = xdash + centre[0];
					ydash = ydash + centre[1];
					zdash = zdash + centre[2];

					// Add marker
					addMarker(xdash, ydash, zdash, markerID);
					markerID++;
				}
			}
		}
	}


#else

	// Build square
	for (int i = 0; i < numMarkersLength; i++) {
		for (int j = 0; j < numMarkersHeight; j++) {

			// x and y position
			x = -length / 2.0 + i * spacingLength;
			y = -height / 2.0 + j * spacingHeight;
			z = 0.0;

			// Only add the marker if the point is on an edge
			if ((i == 0 || i == numMarkersLength - 1) || (j == 0 || j == numMarkersHeight - 1)) {

				// Transform x y and z based on rotation
				xdash = x * cos(angles[0] * L_PI / 180) - y * sin(angles[0] * L_PI / 180);
				ydash = x * sin(angles[0] * L_PI / 180) + y * cos(angles[0] * L_PI / 180);
				zdash = z;
				xdash = xdash + centre[0];
				ydash = ydash + centre[1];
				zdash = zdash + centre[2];

				// Add marker
				addMarker(xdash, ydash, zdash, markerID);
				markerID++;
			}
		}
	}
#endif

	// Get rank
	int rank = GridUtils::safeGetRank();

	// Delete markers which exist off rank
	if (rank != owningRank) {
		*GridUtils::logfile << "Deleting markers which are not on this rank..." << std::endl;
		deleteOffRankMarkers();
	}
};



/*********************************************/
/// \brief	Custom constructor for building plate
///
/// \param 	g			hierarchy pointer to grid hierarchy
/// \param 	bodyID		ID of body in array of bodies
/// \param centre		centre point of plate
/// \param length		length of plate
/// \param width		width of plate
/// \param angles		angle of plate
template <typename MarkerType>
Body<MarkerType>::Body(GridObj* g, int bodyID, std::vector<double> &centre,
	double length, double width, std::vector<double> &angles) {

	// Set the body base class parameters from constructor inputs
	this->_Owner = g;
	this->id = bodyID;
	this->closed_surface = false;

	// Set level
	this->level = _Owner->level;

	// Set the rank which owns this body
	this->owningRank = assignOwningRank(id);

	// Get number of markers in each direction
	int numMarkersLength = static_cast<int>(std::floor(length / g->dh)) + 1;
	int numMarkersWidth = static_cast<int>(std::floor(width / g->dh)) + 1;

	// Get spacing
	double spacingLength = length / (numMarkersLength - 1);
	double spacingWidth = width / (numMarkersWidth - 1);

	// Get start point
	double xStart = -length / 2.0;
	double yStart = 0.0;
	double zStart = -width / 2.0;

	double thetaX = angles[eXDirection] * L_PI / 180.0;
	double thetaY = angles[eYDirection] * L_PI / 180.0;
	double thetaZ = angles[eZDirection] * L_PI / 180.0;

	// Set Tz rotation matrix
	std::vector<std::vector<double>> Tx = {{1.0, 0.0, 0.0},
										   {0.0, cos(thetaX), -sin(thetaX)},
										   {0.0, sin(thetaX), cos(thetaX)}};

	// Set Tz rotation matrix
	std::vector<std::vector<double>> Ty = {{cos(thetaY), 0.0, sin(thetaY)},
										   {0.0, 1.0, 0.0},
										   {-sin(thetaY), 0.0, cos(thetaY)}};

	// Set Tz rotation matrix
	std::vector<std::vector<double>> Tz = {{cos(thetaZ), -sin(thetaZ), 0.0},
										   {sin(thetaZ), cos(thetaZ), 0.0},
										   {0.0, 0.0, 1.0}};

	// Get rotation matrix
	std::vector<std::vector<double>> R =  GridUtils::matrix_multiply(GridUtils::matrix_multiply(Tz, Ty), Tx);


	// Marker ID
	int markerID = 0;

	// Now loop through all markers
	std::vector<double> position(L_DIMS, 0);
	for (int i = 0; i < numMarkersLength; i++) {
		for (int j = 0; j < numMarkersWidth; j++) {

			// Get position
			position[eXDirection] = xStart + i * spacingLength;
			position[eYDirection] = yStart;
			position[eZDirection] = zStart + j * spacingWidth;

			// Rotate about z and x axis
			position = GridUtils::matrix_multiply(Tz, position);
			position = GridUtils::matrix_multiply(Tx, position);

			// Add the centre position
			position[eXDirection] += centre[eXDirection];
			position[eYDirection] += centre[eYDirection];
			position[eZDirection] += centre[eZDirection];

			// Add marker
			addMarker(position[eXDirection], position[eYDirection], position[eZDirection], markerID);
			markerID++;
		}
	}


	// Get rank
	int rank = GridUtils::safeGetRank();

	// Delete markers which exist off rank
	if (rank != owningRank) {
		*GridUtils::logfile << "Deleting markers which are not on this rank..." << std::endl;
		deleteOffRankMarkers();
	}
}

/*********************************************/
/// \brief	Add marker to the body
///
/// \param	x			global X-position of marker
/// \param	y			global Y-position of marker
/// \param	z 			global Z-position of marker
/// \param 	markerID	rank independent ID of marker within body
template <typename MarkerType>
void Body<MarkerType>::addMarker(double x, double y, double z, int markerID)
{

	// Add a new marker object to the array
	markers.emplace_back(x, y, z, markerID, _Owner);

};

/*********************************************/
/// \brief	Delete markers which have been built but exist off rank.
template <typename MarkerType>
void Body<MarkerType>::deleteRecvLayerMarkers()
{

	// Loop through markers in body (if any) and delete ones which are on receiver layer
	if (this->markers.size() > 0) {
		int a = 0;
		do {
			// If on receiver layer then delete that marker
			if (GridUtils::isOnRecvLayer(
				this->markers[a].position[eXDirection],
				this->markers[a].position[eYDirection],
				this->markers[a].position[eZDirection]))
			{
				this->markers.erase(this->markers.begin() + a);
			}
			// If not, keep and move onto next one
			else {

				// Increment counter
				a++;
			}
		} while (a < static_cast<int>(this->markers.size()));
	}
}

/*********************************************/
/// \brief	Delete markers which have been built but exist off rank
template <typename MarkerType>
void Body<MarkerType>::deleteOffRankMarkers()
{

	// Loop through markers in body and delete ones which are not on this rank
	eLocationOnRank loc = eNone;
	int a = 0;
	do {
		// If not on rank then delete that marker
		if (!GridUtils::isOnThisRank(
			this->markers[a].position[eXDirection],
			this->markers[a].position[eYDirection],
			this->markers[a].position[eZDirection],
			&loc, this->_Owner))
		{
			this->markers.erase(this->markers.begin() + a);
		}
		// If it is, keep and move onto next one
		else {

			// Increment counter
			a++;
		}
	} while (a < static_cast<int>(this->markers.size()));
};

/*********************************************/
/// \brief	Retrieve marker data
///
///			Return marker whose primary support data is nearest the
///			supplied global position.
///
/// \param	x 			X-position nearest to marker to be retrieved
/// \param	y 			Y-position nearest to marker to be retrieved
/// \param	z 			Z-position nearest to marker to be retrieved
/// \returns			marker data structure returned - if no marker found structure is marked as invalid by assigning an ID of -1
template <typename MarkerType>
MarkerData* Body<MarkerType>::getMarkerData(double x, double y, double z) {

	// Get indices of voxel associated with the supplied position
	std::vector<int> vox;
	eLocationOnRank *loc = nullptr;
	if (GridUtils::isOnThisRank(x, y, z, loc, _Owner, &vox))
	{

		// Find marker whose primary support point matches these indices
		for (int i = 0; i < static_cast<int>(markers.size()); ++i) {
			if (markers[i].supp_i[0] == vox[0] &&
				markers[i].supp_j[0] == vox[1] &&
				markers[i].supp_k[0] == vox[2]) {

				// Indice represents the target ID so create new MarkerData store on the heap
				MarkerData* m_MarkerData = new MarkerData(
					markers[i].supp_i[0],
					markers[i].supp_j[0],
					markers[i].supp_k[0],
					markers[i].position[0],
					markers[i].position[1],
					markers[i].position[2],
					i
					);
				return m_MarkerData;	// Return the pointer to store information
			}

		}
	}

	// If not found then create empty MarkerData store using default constructor
	MarkerData* m_MarkerData = new MarkerData();
	return m_MarkerData;

};

/*********************************************/
/// \brief	Downsampling voxel-grid filter to take a point and add it to current body
///
///			This method attempts to add a marker to body at the global location 
///			but obeys the rules of a 1 marker per cell voxel-grid filter to 
///			ensure markers are distributed such that their spacing roughly matches 
///			the background lattice. It is usually called inside a loop and requires
///			a few extra pieces of information to be tracked throughout.
///
/// \param	x			desired global X-position of new marker
/// \param	y			desired global Y-position of new marker
/// \param	z			desired global Z-position of new marker
///	\param	markerID	requested rank independent ID of marker within body
/// \param	curr_mark	is a reference to the index of last marker added
///	\param	counter		is a reference to the total number of markers in the body
template <typename MarkerType>
void Body<MarkerType>::passToVoxelFilter(double x, double y, double z, int markerID, int& curr_mark, std::vector<int>& counter) {

	// If point in current voxel
	if (isInVoxel(x, y, z, curr_mark)) {

		// Increment point counter
		counter[curr_mark]++;

		// Update position of marker in current voxel
		markers[curr_mark].position[0] =
			((markers[curr_mark].position[0] * (counter[curr_mark] - 1)) + x) / counter[curr_mark];
		markers[curr_mark].position[1] =
			((markers[curr_mark].position[1] * (counter[curr_mark] - 1)) + y) / counter[curr_mark];
		markers[curr_mark].position[2] =
			((markers[curr_mark].position[2] * (counter[curr_mark] - 1)) + z) / counter[curr_mark];


		// If point is in an existing voxel
	}
	else if (isVoxelMarkerVoxel(x, y, z)) {

		// Recover voxel number
		MarkerData* m_data = getMarkerData(x, y, z);
		curr_mark = m_data->ID;

		// Increment point counter
		counter[curr_mark]++;

		// Update position of marker in current voxel
		markers[curr_mark].position[0] =
			((markers[curr_mark].position[0] * (counter[curr_mark] - 1)) + x) / counter[curr_mark];
		markers[curr_mark].position[1] =
			((markers[curr_mark].position[1] * (counter[curr_mark] - 1)) + y) / counter[curr_mark];
		markers[curr_mark].position[2] =
			((markers[curr_mark].position[2] * (counter[curr_mark] - 1)) + z) / counter[curr_mark];

		delete m_data;
			
	}
	// Must be in a new marker voxel
	else {

		// Reset counter and increment voxel index
		curr_mark = static_cast<int>(counter.size());
		counter.push_back(1);

		// Create new marker as this is a new marker voxel
		addMarker(x, y, z, markerID);
	}

};

/*********************************************/
/// \brief	Determines whether a point is inside another marker's support voxel
///
///			Typically called indirectly by the voxel-grid filter method and not directly.
///
/// \param	x				X-position of point
/// \param	y				Y-position of point
/// \param	z				Z-position of point
/// \param	curr_mark		ID of the marker
/// \returns				true or false
template <typename MarkerType>
bool Body<MarkerType>::isInVoxel(double x, double y, double z, int curr_mark) {

	try {

		// Try to retrieve the position of the voxel centre belonging to <curr_mark>
		// Assume that first support point is the closest to the marker position
		double vx = _Owner->XPos[markers[curr_mark].supp_i[0]];
		double vy = _Owner->YPos[markers[curr_mark].supp_j[0]];
		double vz = _Owner->ZPos[markers[curr_mark].supp_k[0]];

		// Test within
		if ((x >= vx - (_Owner->dh / 2.0) && x < vx + (_Owner->dh / 2.0)) &&
			(y >= vy - (_Owner->dh / 2.0) && y < vy + (_Owner->dh / 2.0)) &&
			(z >= vz - (_Owner->dh / 2.0) && z < vz + (_Owner->dh / 2.0))
			) return true;

		// Catch all
	}
	catch (...) {

		// If failed, marker probably doesn't exist so return false
		return false;

	}

	return false;

};

/*********************************************/
/// \brief	Determines whether a point is inside an existing marker's support voxel
///
///			Typically called indirectly by the voxel-grid filter method and not directly.
///
/// \param	x				X-position of point
/// \param	y				Y-position of point
/// \param	z				Z-position of point
/// \returns				true or false
template <typename MarkerType>
bool Body<MarkerType>::isVoxelMarkerVoxel(double x, double y, double z) {

	// Try get the MarkerData store
	MarkerData* m_data = getMarkerData(x, y, z);

	// True if the data store is not empty
	if (m_data->isValid()) {

		delete m_data;	// Deallocate the store before leaving scope
		return true;

	}

	delete m_data;
	return false;

};


/*********************************************/
/// \brief	Method to build a body from point cloud data
///
/// \param	_PCpts	point cloud data from which to build body
template <typename MarkerType>
void Body<MarkerType>::buildFromCloud(PCpts *_PCpts)
{

	// Declare local variables
	std::vector<int> locals;

	// Voxel grid filter //

	*GridUtils::logfile << "ObjectManager: Applying voxel grid filter..." << std::endl;

	// Place first marker
	if (!_PCpts->x.empty())
		addMarker(_PCpts->x[0], _PCpts->y[0], _PCpts->z[0], _PCpts->id[0]);

	// Increment counters
	int curr_marker = 0;
	std::vector<int> counter;
	counter.push_back(1);

	// Loop over array of points
	for (size_t a = 1; a < _PCpts->x.size(); a++)
	{
		// Pass to point builder
		passToVoxelFilter(_PCpts->x[a], _PCpts->y[a], _PCpts->z[a], _PCpts->id[a], curr_marker, counter);
	}

	*GridUtils::logfile << "ObjectManager: Object represented by " << std::to_string(markers.size()) <<
		" markers using 1 marker / voxel voxelisation." << std::endl;
};

// ************************************************************************** //
/// \brief	VTK body writer
///
/// \param	tval	time value at which the write out is being performed.
template <typename MarkerType>
void Body<MarkerType>::writeVtkPosition(int tval)
{
	// Get rank
	int rank = GridUtils::safeGetRank();

	// Create file name then output file stream
	std::stringstream fileName;
	fileName << GridUtils::path_str + "/vtk_out.Body" << id << "." << std::to_string(rank) << "." << (int)tval << ".vtk";

	std::ofstream fout;
	fout.open(fileName.str().c_str());

	// Add header information
	fout << "# vtk DataFile Version 3.0f\n";
	fout << "Output for body ID " << id << ", rank " << std::to_string(rank) << " at time t = " << (int)tval << "\n";
	fout << "ASCII\n";
	fout << "DATASET POLYDATA\n";

	// Write out the positions of each Lagrange marker
	fout << "POINTS " << markers.size() << " float\n";
	for (auto m : markers) {

		fout << m.position[0] << " "
			 << m.position[1] << " "
			 << m.position[2] << std::endl;
	}

	// Write out the connectivity of each Lagrange marker
	size_t nLines;
	if (closed_surface == true)
		nLines = markers.size();
	else
		nLines = markers.size() - 1;

	// Write out number of lines
	fout << "LINES " << nLines << " " << 3 * nLines << std::endl;

	// Write out connectivity
	for (size_t i = 0; i < nLines; i++) {
		fout << 2 << " " << i << " " << (i + 1) % markers.size() << std::endl;
	}

	// Close file
	fout.close();
};


/*********************************************/
/// \brief	Assigns owning rank of body based on which grids each rank has access to
///
/// \param	id	global body ID
/// \returns	rank number
template <typename MarkerType>
int Body<MarkerType>::assignOwningRank(int id) {

	// If serial just return 0
#ifndef L_BUILD_FOR_MPI
	return 0;
#else

	// Some work required if in parallel
	MpiManager *mpim = MpiManager::getInstance();

	// Vector of valid ranks which are allowed to own it
	std::vector<int> validRanks;

	// Set pointer to current grid owner
	GridObj *g = NULL;
	GridUtils::getGrid(0, 0, g);

	// Loop through rankGrids to see which ranks are allowed to own it
	for (int rank = 0; rank < mpim->num_ranks; rank++) {

		// Check if this rank has this level
		if (mpim->rankGrids[rank] >= level)
			validRanks.push_back(rank);
	}

	// Return
	return validRanks[id % validRanks.size()];
#endif
};

#endif
