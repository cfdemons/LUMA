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
#ifndef BODY_H
#define BODY_H

#include "stdafx.h"
class GridObj;
#include "GridUtils.h"

/// \brief	Container class to hold marker information.
class MarkerData {

public:

	/// \brief Constructor.
	/// \param i i-index of primary support site
	/// \param j j-index of primary support site
	/// \param k k-index of primary support site
	/// \param x x-position of marker
	/// \param y y-position of marker
	/// \param z z-position of marker
	/// \param ID marker number in a given body
	MarkerData(int i, int j, int k, double x, double y, double z, int ID) {

		// Custom constructor
		this->i = i;
		this->j = j;
		this->k = k;
		this->x = x;
		this->y = y;
		this->z = z;
		this->ID = ID;

	};

	/// \brief	Default Constructor.
	///
	///			Initialise with invalid marker indicator which 
	///			is to set the x position to NaN.
	MarkerData(void) {

		// Essentially null a double in the store making it invalid
		this->x = std::numeric_limits<double>::quiet_NaN();

	};

	/// Default destructor
	~MarkerData(void) {};

	// Voxel indices
	int i;	///< i-index of primary support site
	int j;	///< j-index of primary support site
	int k;	///< k-index of primary support site

	/// Marker ID (position in array of markers)
	int ID;

	// Marker position
	double x;	///< x-position of marker
	double y;	///< y-position of marker
	double z;	///< z-position of marker

};



/// \brief	Generic body class.
///
///			Can consist of any type of Marker so templated.
template <typename MarkerType>
class Body
{


public:

	Body(void);			// Default Constructor
	~Body(void);		// Default destructor
	Body(GridObj* g, size_t id);	// Custom constructor assigning grid hierarchy

	// ************************ Members ************************ //

protected:
	double spacing;						///< Spacing of the markers in physical units
	std::vector<MarkerType> markers;	///< Array of markers which make up the body
	bool closed_surface;				///< Flag to specify whether or not it is a closed surface (for output)
	GridObj* _Owner;					///< Pointer to owning grid
	size_t id;							///< ID of body in array of bodies


	// ************************ Methods ************************ //

	void addMarker(double x, double y, double z);					// Add a marker
	MarkerData* getMarkerData(double x, double y, double z);		// Retireve nearest marker data
	virtual void markerAdder(double x, double y, double z,
		int& curr_mark, std::vector<int>& counter, bool flag);		// Voxelising marker adder (can be overridden)
	bool isInVoxel(double x, double y, double z, int curr_mark);	// Check a point is inside an existing marker voxel
	bool isVoxelMarkerVoxel(double x, double y, double z);			// Check whether nearest voxel is a marker voxel


};


// ************************ Implementation of Body methods ************************ //

/// Default Constructor
template <typename MarkerType>
Body<MarkerType>::Body(void)
{
};
/// Default destructor
template <typename MarkerType>
Body<MarkerType>::~Body(void)
{
};

/// \brief Custom constructor setting owning grid.
///
/// \param g pointer to grid which owns this body.
/// \param id indicates position of body in array of bodies.
template <typename MarkerType>
Body<MarkerType>::Body(GridObj* g, size_t id)
{
	this->_Owner = g;
	this->id = id;
};

/// \brief	Add marker to the body.
/// \param	x global X-position of marker.
/// \param	y global Y-position of marker.
/// \param	z global Z-position of marker.
template <typename MarkerType>
void Body<MarkerType>::addMarker(double x, double y, double z)
{
	// Add a new marker object to the array
	markers.emplace_back(x, y, z);

	// Add nearest node as basic support
	std::vector<int> globals = GridUtils::getVoxInd(x, y, z, _Owner);
	this->markers.back().supp_i.push_back(globals[0]);
	this->markers.back().supp_j.push_back(globals[1]);
	this->markers.back().supp_k.push_back(globals[2]);
};

/*********************************************/
/// \brief	Retrieve marker data.
///
///			Return marker and voxel/primary support data associated with
///			supplied global position.
///
/// \param x global X-position nearest to marker to be retrieved.
/// \param y global Y-position nearest to marker to be retrieved.
/// \param z global Z-position nearest to marker to be retrieved.
/// \return MarkerData marker data structure returned. If no marker found, structure is marked as invalid.
template <typename MarkerType>
MarkerData* Body<MarkerType>::getMarkerData(double x, double y, double z) {

	// Get global indices of voxel associated with the supplied position
	std::vector<int> vox = GridUtils::getVoxInd(x, y, z, _Owner);

	// Find all marker IDs that match these indices
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

	// If not found then create empty MarkerData store using default constructor
	MarkerData* m_MarkerData = new MarkerData();
	return m_MarkerData;

};

/*********************************************/
/// \brief	Downsampling marker adding method
///
///			This method tries to add a marker to body at the global location 
///			given but obeys the rules of a voxel-grid filter to ensure markers are
///			distributed such that their spacing roughly matches the 
///			background lattice. Making it virtual allows different attirbutes
///			to be passed to the addMarker() methods to allow for different
///			marker types.
///
/// \param x desired global X-position of new marker.
/// \param y desired globalY-position of new marker.
/// \param z desired globalZ-position of new marker.
/// \param curr_mark is a reference to the ID of last marker.
///	\param counter is a reference to the total number of markers in the body.
/// \param flag is an optional flag argument this is used in overrides of this method.
template <typename MarkerType>
void Body<MarkerType>::markerAdder(double x, double y, double z, int& curr_mark, std::vector<int>& counter, bool flag) {

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
		addMarker(x, y, z);

	}


};

/*********************************************/
/// \brief	Determines whether a point is inside another marker's support voxel.
/// \param x X-position of point.
/// \param y Y-position of point.
/// \param z Z-position of point.
/// \param curr_mark ID of the marker.
/// \return true of false
template <typename MarkerType>
bool Body<MarkerType>::isInVoxel(double x, double y, double z, int curr_mark) {

	try {

		// Try to retrieve the position of the marker identified by <curr_mark>
		// Assume that first support point is the closest to the marker position
		int i = markers[curr_mark].supp_i[0];
		int j = markers[curr_mark].supp_j[0];
		int k = markers[curr_mark].supp_k[0];

		// Test within
		if ((x > i - (_Owner->dx / 2) && x <= i + (_Owner->dx / 2)) &&
			(y > j - (_Owner->dx / 2) && y <= j + (_Owner->dx / 2)) &&
			(z > k - (_Owner->dx / 2) && z <= k + (_Owner->dx / 2))
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
/// \brief	Determines whether a point is inside an existing marker's support voxel.
/// \param x global X-position of point.
/// \param y global Y-position of point.
/// \param z global Z-position of point.
/// \return true of false
// Returns boolean as to whether a given point is in an existing marker voxel
template <typename MarkerType>
bool Body<MarkerType>::isVoxelMarkerVoxel(double x, double y, double z) {

	// Try get the MarkerData store
	MarkerData* m_data = getMarkerData(x, y, z);

	// True if the data store is not empty
	if (!L_IS_NAN(m_data->x)) {

		delete m_data;	// Deallocate the store
		return true;

	}

	return false;

};

#endif