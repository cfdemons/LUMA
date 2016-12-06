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

/* This file contains the constructors for the Grid objects. */

#include "../inc/stdafx.h"
#include "../inc/GridObj.h"
#include "../inc/MpiManager.h"
#include "../inc/GridUtils.h"

// Static declarations
std::ofstream* GridUtils::logfile;
int MpiManager::my_rank;
int MpiManager::num_ranks;

// ****************************************************************************
/// Default Constructor
GridObj::GridObj(void)
{
}

// ****************************************************************************
/// Default Destructor
GridObj::~GridObj(void)
{
}

// ****************************************************************************
/// \brief	Constructor for top level grid.
///
///			Coarse limits are set to zero and then L0-specific initialiser called.
///
/// \param level always should be zero astop level grid.
GridObj::GridObj(int level)
{
	// Defaults
	this->t = 0;

	// Assign level and region number
	this->level = level;
	this->region_number = 0;
	// Set limits of refinement to zero as top level
	for (int i = 0; i < 2; i++) {
		this->CoarseLimsX[i] = 0;
		this->CoarseLimsY[i] = 0;
		this->CoarseLimsZ[i] = 0;
	}

	// Assign refinement ratio
	this->refinement_ratio = (1.0 / pow(2.0, static_cast<double>(level)));

	// Reset timers
	this->timeav_mpi_overhead = 0.0;
	this->timeav_timestep = 0.0;

	*GridUtils::logfile << "Constructing Grid level " << level  << std::endl;

	// Call L0 non-MPI initialiser
	this->LBM_initGrid();

}

// ****************************************************************************
/// \brief	MPI constructor for top level grid.
///
///			When using MPI, this constructors a local grid which represents an 
///			appropriate portion of the top-level grid as dictated by the extent
///			of this rank.
///
/// \param level always should be zero astop level grid.
/// \param local_size vector indicating dimensions of local grid including halo.
/// \param GlobalLimsInd vector indicating the global indices of the edges of 
///			this local grid.
/// \param GlobalLimsPos vector indicating the global positions of the edges of 
///			this local grid.
GridObj::GridObj(int level, std::vector<int> local_size, 
				 std::vector< std::vector<int> > GlobalLimsInd, 
				 std::vector< std::vector<double> > GlobalLimsPos)
{
	// Assign
	this->level = level;
    this->region_number = 0;
	this->t = 0;
	// Set limits of refinement to zero as top level
	for (int i = 0; i < 2; i++) {
		this->CoarseLimsX[i] = 0;
		this->CoarseLimsY[i] = 0;
		this->CoarseLimsZ[i] = 0;
	}

	// Assign refinement ratio
	this->refinement_ratio = (1.0 / pow(2.0, static_cast<double>(level)));

	// Reset timers
	this->timeav_mpi_overhead = 0.0;
	this->timeav_timestep = 0.0;

	*GridUtils::logfile << "Constructing Grid level " << level << " on rank " << MpiManager::my_rank << std::endl;
	
	// Call MPI initialisation routine
	LBM_initGrid(local_size, GlobalLimsInd, GlobalLimsPos); 

}

// ****************************************************************************
/// \brief Constructor for a sub-grid.
/// \param RegionNumber ID indicating the region of nested refinement to which 
///			this sub-grid belongs.
/// \param pGrid pointer to parent grid.
GridObj::GridObj(int RegionNumber, GridObj& pGrid)
{

	// Assign
	this->level = pGrid.level + 1;
    this->region_number = RegionNumber;
	this->t = 0;
	this->parentGrid = &pGrid;

	// Assign refinement ratio
	this->refinement_ratio = (1.0 / pow(2.0, static_cast<double>(level)));

	// Reset timers
	this->timeav_mpi_overhead = 0.0;
	this->timeav_timestep = 0.0;
	
	*GridUtils::logfile << "Constructing Sub-Grid level " << level << ", region " << region_number << ", on rank " << MpiManager::my_rank << std::endl;

}

// ****************************************************************************
/// \brief Wrapper method to add sub-grid to this grid.
/// \param RegionNumber ID indicating the region of nested refinement to which 
///			this sub-grid belongs.
void GridObj::LBM_addSubGrid(int RegionNumber) {

	// Return if no subgrid required for the current grid
	if ( !GridUtils::hasThisSubGrid(*this, RegionNumber) ) return;
	
	// Ok to proceed and add the subgrid passing parent grid as a reference for initialisation
	subGrid.emplace_back( RegionNumber, *this);	
	
	// Initialise the subgrid passing position of corner of the refined region on parent grid
	this->subGrid.back().LBM_initSubGrid(*this);

	// Add another subgrid beneath the one just created if necessary
	if (this->subGrid.back().level < L_NUM_LEVELS) {
		this->subGrid.back().LBM_addSubGrid(this->subGrid.back().region_number);
	}

}

// ****************************************************************************
// ****************************************************************************
// Other member methods are in their own files prefixed GridObj_