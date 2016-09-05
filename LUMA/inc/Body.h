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

#pragma once
#include <vector>

class GridObj;

// Class representing return structure for marker query
class MarkerData {

public:

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
	MarkerData(void) {

		// Essentially null a double in the store making it invalid
		this->x = std::numeric_limits<double>::quiet_NaN();

	};
	~MarkerData(void) {};

	// Voxel indices
	int i;
	int j;
	int k;

	// Marker ID (position in array of markers)
	int ID;

	// Marker position
	double x;
	double y;
	double z;

};



// Markers may be of different types depending on the body
template <typename MarkerType>

// Class representing a general body
class Body
{

	// Default Constructor / Destructor //
public:
	Body(void)
	{
	};
	~Body(void)
	{
	};

	Body(GridObj* g)
	{
		this->_Owner = g;
	};

	// Protected Members //
	
protected:
	double spacing;						// Spacing of the markers in physical units
	std::vector<MarkerType> markers;	// Array of markers which make up the body
	bool closed_surface;				// Flag to specify whether or not it is a closed surface (for output)
	GridObj* _Owner;					// Pointer to owning grid
	

	// ************************ Methods ************************ //

	// Add marker to the body //
	void addMarker(double x, double y, double z)
	{
		markers.emplace_back(x, y, z);	// Add a new marker object to the array
	}

	/*********************************************/
	// Return marker and voxel data associated with global position supplied //
	MarkerData* getMarkerData(double x, double y, double z) {

		// Get indices of voxel associated with the supplied position
		std::vector<int> vox = ObjectManager::getInstance()->getVoxInd(x, y, z);

		// Find all marker IDs that match these indices
		for (int i = 0; i < static_cast<int>(markers.size()); ++i) {
			if (markers[i].supp_i[0] == vox[0] &&
				markers[i].supp_j[0] == vox[1] &&
				markers[i].supp_k[0] == vox[2]) {

				// Indice represents the target ID so create new MarkerData store on the heap
				MarkerData* m_MarkerData = new MarkerData(markers[i].supp_i[0],
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

	}

	/*********************************************/
	// Routine to add markers to the current body //
	// (pass the marker index and the counter by reference so they are updated in the calling scope) //
	void markerAdder(double x, double y, double z, int& curr_mark, std::vector<int>& counter) {

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


	}

	/*********************************************/
	// Returns boolean as to whether a given point is in the voxel associated with <curr_mark>
	bool isInVoxel(double x, double y, double z, int curr_mark) {

		try {

			// Try to retrieve the position of the marker identified by <curr_mark>
			// Assume that first support point is the closest to the marker position
			int i = markers[curr_mark].supp_i[0];
			int j = markers[curr_mark].supp_j[0];
			int k = markers[curr_mark].supp_k[0];

			// Test within
			if ((x > i - 0.5 && x <= i + 0.5) &&
				(y > j - 0.5 && y <= j + 0.5) &&
				(z > k - 0.5 && z <= k + 0.5)
				) return true;

			// Catch all
		}
		catch (...) {

			// If failed, marker probably doesn't exist so return false
			return false;

		}

		return false;

	}

	/*********************************************/
	// Returns boolean as to whether a given point is in an existing marker voxel
	bool isVoxelMarkerVoxel(double x, double y, double z) {

		// Try get the MarkerData store
		MarkerData* m_data = getMarkerData(x, y, z);

		// True if the data store is not empty
		if (!is_nan(m_data->x)) {

			delete m_data;	// Deallocate the store
			return true;

		}

		return false;

	}

};
