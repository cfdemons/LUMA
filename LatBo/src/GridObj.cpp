/* This file contains the constructors for the Grid objects. */

#include "../inc/stdafx.h"
#include "../inc/definitions.h"
#include "../inc/GridObj.h"
#include "../inc/MpiManager.h"

// Static declarations
std::ofstream* GridUtils::logfile;
int MpiManager::my_rank;
int MpiManager::num_ranks;

// ***************************************************************************************************
// Default Constructor
GridObj::GridObj(void)
{
}

// ***************************************************************************************************
// Default Destructor
GridObj::~GridObj(void)
{
}

// ***************************************************************************************************
// Basic Constructor for top level grid
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

	// Reset timers
	this->timeav_mpi_overhead = 0.0;
	this->timeav_timestep = 0.0;

	*GridUtils::logfile << "Constructing Grid level " << level  << std::endl;

	// Call L0 non-MPI initialiser
	this->LBM_init_grid();

}

// ***************************************************************************************************
// MPI constructor for level 0 with grid level, rank, local grid size and its global edges
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

	// Reset timers
	this->timeav_mpi_overhead = 0.0;
	this->timeav_timestep = 0.0;

	*GridUtils::logfile << "Constructing Grid level " << level << " on rank " << MpiManager::my_rank << std::endl;
	
	// Call MPI initialisation routine
	LBM_init_grid(local_size, GlobalLimsInd, GlobalLimsPos); 

}

// ***************************************************************************************************
// Overloaded constructor for a sub grid
GridObj::GridObj(int RegionNumber, GridObj& pGrid)
{

	// Assign
	this->level = pGrid.level + 1;
    this->region_number = RegionNumber;
	this->t = 0;

	// Reset timers
	this->timeav_mpi_overhead = 0.0;
	this->timeav_timestep = 0.0;
	
	*GridUtils::logfile << "Constructing Sub-Grid level " << level << ", region " << region_number << ", on rank " << MpiManager::my_rank << std::endl;

}

// ***************************************************************************************************
// Method to generate a subgrid
void GridObj::LBM_addSubGrid(int RegionNumber) {

	// Return if no subgrid required for the current grid
	if ( !GridUtils::hasThisSubGrid(*this, RegionNumber) ) return;
	
	// Ok to proceed and add the subgrid passing parent grid as a reference for initialisation
	subGrid.emplace_back( RegionNumber, *this);	
	
	// Initialise the subgrid passing position of corner of the refined region on parent grid
	this->subGrid.back().LBM_init_subgrid(*this);

	// Add another subgrid beneath the one just created if necessary
	if (this->subGrid.back().level < NumLev) {
		this->subGrid.back().LBM_addSubGrid(this->subGrid.back().region_number);
	}

}

// ***************************************************************************************************
// ***************************************************************************************************
// Other member methods are in their own files prefixed GridObj_