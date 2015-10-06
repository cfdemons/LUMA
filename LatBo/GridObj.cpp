/* This file contains the constructors for the Grid object.
*/

#include "stdafx.h"
#include "definitions.h"
#include "GridObj.h"
#include <vector>

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
GridObj::GridObj(int level, std::ofstream* logfile)
{
	// Defaults
	this->my_rank = 0;
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

	// Create utility class
	this->gUtils.setLogFile(logfile);

	*logfile << "Constructing Grid level " << level  << std::endl;

	this->LBM_init_grid(); // Call L0 non-MPI initialiser

}

// ***************************************************************************************************
// MPI constructor for level 0 with grid level, rank, local grid size and its global edges
GridObj::GridObj(int level, int rank, std::vector<unsigned int> local_size, 
				 std::vector< std::vector<unsigned int> > GlobalLimsInd, 
				 std::vector< std::vector<double> > GlobalLimsPos,
				 std::ofstream* logfile)
{
	// Assign
	this->level = level;
    this->region_number = 0;
	this->my_rank = rank;
	this->t = 0;
	// Set limits of refinement to zero as top level
	for (int i = 0; i < 2; i++) {
		this->CoarseLimsX[i] = 0;
		this->CoarseLimsY[i] = 0;
		this->CoarseLimsZ[i] = 0;
	}

	// Logfile assign
	this->gUtils.setLogFile(logfile);

	*logfile << "Building Grid level " << level << " on rank " << my_rank << std::endl;

	LBM_init_grid( local_size, GlobalLimsInd, GlobalLimsPos ); // Run initialisation routine

}

// ***************************************************************************************************
// Overload constructor for MPI sub grid
GridObj::GridObj(int level, int RegionNumber, int rank, std::ofstream* logfile)
{
	// Assign
	this->level = level;
    this->region_number = RegionNumber;
	this->my_rank = rank;
	this->t = 0;

	// Logfile assign
	this->gUtils.setLogFile(logfile);

	*logfile << "Building Grid level " << level << ", region " << region_number << ", on rank " << my_rank << std::endl;

}

// ***************************************************************************************************
// Method to generate a subgrid
void GridObj::LBM_addSubGrid(int RegionNumber) {

	subGrid.emplace_back( this->level + 1, RegionNumber, this->my_rank, this->gUtils.logfile);
	
	// Update limits as cannot be done through the constructor since
	// we cannot pass more than 5 arguments to emplace_back()
	if (this->level == 0) {

		// First layer of sub-grids get limits from this top level grid -- convert to local indices
		this->subGrid[subGrid.size()-1].CoarseLimsX[0] = RefXstart[RegionNumber] - this->XInd[1] + 1;	this->subGrid[subGrid.size()-1].CoarseLimsX[1] = RefXend[RegionNumber] - this->XInd[1] + 1;
		this->subGrid[subGrid.size()-1].CoarseLimsY[0] = RefYstart[RegionNumber] - this->YInd[1] + 1;	this->subGrid[subGrid.size()-1].CoarseLimsY[1] = RefYend[RegionNumber] - this->YInd[1] + 1;
#if (dims == 3)
		this->subGrid[subGrid.size()-1].CoarseLimsZ[0] = RefZstart[RegionNumber] - this->ZInd[1] + 1;	this->subGrid[subGrid.size()-1].CoarseLimsZ[1] = RefZend[RegionNumber] - this->ZInd[1] + 1;
#else
		this->subGrid[subGrid.size()-1].CoarseLimsZ[0] = RefZstart[RegionNumber];	this->subGrid[subGrid.size()-1].CoarseLimsZ[1];
#endif
	
	} else {
		
		// Lower sub-grids get limits from higher sub-grid -- already local
		this->subGrid[subGrid.size()-1].CoarseLimsX[0] = 2; this->subGrid[subGrid.size()-1].CoarseLimsX[1] = this->XInd.size() - 2;
		this->subGrid[subGrid.size()-1].CoarseLimsY[0] = 2; this->subGrid[subGrid.size()-1].CoarseLimsY[1] = this->YInd.size() - 2;
		this->subGrid[subGrid.size()-1].CoarseLimsZ[0] = 2; this->subGrid[subGrid.size()-1].CoarseLimsZ[1] = this->ZInd.size() - 2;
	}

	// Initialise the subgrid passing position of corner of the refined region on parent grid
#if dims == 3
	this->subGrid[subGrid.size()-1].LBM_init_subgrid(	this->XPos[ subGrid[subGrid.size()-1].CoarseLimsX[0] ],
														this->YPos[ subGrid[subGrid.size()-1].CoarseLimsY[0] ],
														this->ZPos[ subGrid[subGrid.size()-1].CoarseLimsZ[0] ],
														this->dx, this->omega, this->mrt_omega	);
#else

    this->subGrid[subGrid.size()-1].LBM_init_subgrid(	this->XPos[ subGrid[subGrid.size()-1].CoarseLimsX[0] ],
														this->YPos[ subGrid[subGrid.size()-1].CoarseLimsY[0] ],
														0.0,
														this->dx, this->omega, this->mrt_omega	);
#endif

	// Add another subgrid beneath the one just created if necessary
	if (this->subGrid[subGrid.size()-1].level < NumLev) {
		this->subGrid[subGrid.size()-1].LBM_addSubGrid(this->subGrid[subGrid.size()-1].region_number);
	}

}

// ***************************************************************************************************
// ***************************************************************************************************
// Other member methods are in their own files prefixed GridObj_