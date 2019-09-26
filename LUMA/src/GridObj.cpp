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

/* This file contains the constructors for the Grid objects. */

#include "../inc/stdafx.h"
#include "../inc/GridObj.h"

// Static declarations
std::ofstream* GridUtils::logfile;

// ****************************************************************************
/// Default Constructor
GridObj::GridObj(void)
{
}

// ****************************************************************************
/// Default Destructor
GridObj::~GridObj(void)
{
	// Loop over subgrid array and destroy each one
	for (GridObj *g : subGrid) if (g) delete g;
}

// ****************************************************************************
/// \brief	Serial build constructor for top level grid.
///
///			Coarse limits are set to zero and then L0-specific initialiser called.
///
/// \param level always should be zero as top level grid.
GridObj::GridObj(int level)
	: t(0), level(level), region_number(0),
	timeav_mpi_overhead(0.0), timeav_timestep(0.0),
	refinement_ratio(1.0 / pow(2.0, static_cast<double>(level)))
{
	// Set limits of refinement to zero as top level
	for (int i = 0; i < 2; i++) {
		this->CoarseLimsX[i] = 0;
		this->CoarseLimsY[i] = 0;
		this->CoarseLimsZ[i] = 0;
	}

	L_INFO("Constructing Grid level " + std::to_string(level) + "...", GridUtils::logfile);

	// Call L0 non-MPI initialiser
	this->LBM_initGrid();
}

// ****************************************************************************
/// \brief	Constructor for a sub-grid.
///
///			This is not called directly but by the addSubGrid() method which 
///			first performs a check to see if a sub-grid is required.
///
/// \param RegionNumber		ID indicating the region of nested refinement to which 
///							this sub-grid belongs.
/// \param pGrid			pointer to parent grid.
GridObj::GridObj(int RegionNumber, GridObj& pGrid)
	: t(0), level(pGrid.level + 1), region_number(RegionNumber),
	parentGrid(&pGrid), refinement_ratio(1.0 / pow(2.0, static_cast<double>(pGrid.level + 1))),
	timeav_mpi_overhead(0.0), timeav_timestep(0.0)
{	
	// Notify user that grid constructor has been called
	L_INFO("Constructing Sub-Grid level " + std::to_string(level) +
		", region " + std::to_string(region_number) + "...", GridUtils::logfile);
}

// ****************************************************************************
/// \brief Wrapper method to add sub-grid to this grid.
/// \param RegionNumber ID indicating the region of nested refinement to which 
///						this sub-grid belongs.
void GridObj::LBM_addSubGrid(int RegionNumber)
{

	/* We check here to see if the extent of this rank covers any part of the 
	 * extent of the sub-grid (once user-specified extent in definitions file
	 * is rounded to the nearest discrete voxel). */

	// Return if no subgrid required for the current grid
	if (!GridUtils::intersectsRefinedRegion(*this, RegionNumber)) return;
	
	// Ok to proceed and add the subgrid passing parent grid as a reference for initialisation
	subGrid.push_back(new GridObj(RegionNumber, *this));	
	
	// Initialise the subgrid passing position of corner of the refined region on parent grid
	this->subGrid.back()->LBM_initSubGrid(*this);

#ifndef L_BUILD_FOR_MPI
	// Update the writable data in the grid manager
	// WHen using MPI this is done when building the communicators.
	GridManager::getInstance()->createWritableDataStore(subGrid.back());
#endif

	// Try to add another subgrid beneath the one just created if necessary
	if (this->subGrid.back()->level < L_NUM_LEVELS)
		this->subGrid.back()->LBM_addSubGrid(this->subGrid.back()->region_number);

}

// ****************************************************************************
// ****************************************************************************
// Other member methods are in their own files prefixed GridObj_
