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

#ifndef GRIDMAN_H
#define GRIDMAN_H

#include "stdafx.h"
class GridObj;
struct HDFstruct;

/// \brief	Grid Manager class.
///
///			Class to manage all information relating to GridObjs in the application.
///			This singleton class may be accessed to supply information about local grids.
class GridManager
{
	friend class MpiManager;
	friend class GridUtils;
	friend class GridObj;

public:
	/// Number of active cells in the calculation
	long activeCellCount;

	/// Number of active cell operations in a coarse time step
	long activeCellOps;

protected:
	/// Pointer to grid hierarchy
	GridObj *Grids;

	// Grid data
	/// \brief	Overall size of each grid (excluding halo of course).
	///
	///			Since L0 can only be region = 0 this array should be accessed as 
	///			[level + region_number * L_NUM_LEVELS] in a loop where level 
	///			cannot be 0. To retrieve L0 info, simply access [0]. The first
	///			index can be accessed using the eCartesianDirection enumeration.
	///
	int global_size[3][L_NUM_LEVELS * L_NUM_REGIONS + 1];

	/// \brief	Absolute position of grid edges (excluding halo of course).
	///
	///			Since L0 can only be region = 0 this array should be accessed as 
	///			[level + region_number * L_NUM_LEVELS] in a loop where level 
	///			cannot be 0. To retrieve L0 info, simply access [0]. The first index
	///			should be accessed using the eCartesianMinMax enumeration.
	///
	double global_edges[6][L_NUM_LEVELS * L_NUM_REGIONS + 1];

	/// \brief	Boolean flag array to indicate the presence of a TL on sub-grid edges.
	///
	///			It is not a given that a sub-grid has a TL on every edge of the grid.
	///			Specifically if we have a sub-grid which is perodic (or in future, which
	///			merges with another sub-grid?). The HDF5 writer needs to know whether 
	///			to exclude sites to account for TL or not so we store information here
	///			from the sub-grid initialisation. The first index should be accessed 
	///			using the enumerator eCartesianMinMax. If no sub-grids present then 
	///			adopts a default 6x1 size to avoid a compilation error.
	///
#if (L_NUM_LEVELS == 0)
	bool subgrid_tlayer_key[6][1];
#else
	bool subgrid_tlayer_key[6][L_NUM_LEVELS * L_NUM_REGIONS];
#endif

	/// Vector indicating periodicity of grid
	bool periodic_flags[3][L_NUM_LEVELS * L_NUM_REGIONS + 1];

	/// Dimensions of coarsest lattice represented on this rank (includes halo if using MPI).
	std::vector<int> local_size;

	/// Vector of structures containing writable region descriptors for block writing (HDF5)
	std::vector<HDFstruct> p_data;

	// METHODS //

public:
	// Singleton design
	static GridManager* getInstance();	// Get the pointer to the singleton instance (create it if necessary)
	static void destroyInstance();
	void setGridHierarchy(GridObj *const grids);


private:
	// Set the local size (MpiManager can set local size)
	void setLocalCoarseSize(const std::vector<int>& size_vector);

	// Allow grid initialisation to store its writable data information
	void createWritableDataStore(HDFstruct *& datastruct);
	bool createWritableDataStore(GridObj const * const targetGrid);

	// Get estimated active cell count within the global bounds supplied
	void updateGlobalCellCount();
	long getActiveCellCount(double *bounds, bool bCountAsOps);
	long getCellCount(int targetLevel, int targetRegion, double *bounds);


private:
	GridManager(void);		///< Private constructor
	~GridManager(void);		///< Private destructor
	static GridManager* me;	///< Pointer to self

};


#endif // GRIDMAN_H