#include "stdafx.h"
#include "GridObj.h"
#include "definitions.h"

// ***************************************************************************************************
// Default Constructor
GridObj::GridObj(void)
{
}

// ***************************************************************************************************
// Deafult Destructor
GridObj::~GridObj(void)
{
}

// ***************************************************************************************************
// Constructor for top level grid
GridObj::GridObj(int level)
{
	// Assign level and region number
	this->level = level;
	this->region_number = 0;
	// Set limits of refinement to zero as top level
	for (int i = 0; i < 2; i++) {
		this->CoarseLimsX[i] = 0;
		this->CoarseLimsY[i] = 0;
		this->CoarseLimsZ[i] = 0;
	}

	std::cout << "Constructing Grid level " << level  << std::endl;
	this->LBM_init_grid(); // Call L0 initialiser

}

// ***************************************************************************************************

// Overload constructor for sub grid
GridObj::GridObj(int level, int RegionNumber)
{
	// Assign level and region number
	this->level = level;
	this->region_number = RegionNumber;
	
	std::cout << "Constructing Grid level " << level << ", region #" << RegionNumber << std::endl;

}

// ***************************************************************************************************
// Method to generate a subgrid
void GridObj::LBM_addSubGrid(int RegionNumber) {

	this->subGrid.emplace_back( this->level + 1, RegionNumber); //  Generate subgrid object

	// Assign limits of refinement immediately
	if (level == 0) { // First subgrid

		this->subGrid[subGrid.size()-1].CoarseLimsX[0] = RefXstart[RegionNumber]; this->subGrid[subGrid.size()-1].CoarseLimsX[1] = RefXend[RegionNumber];
		this->subGrid[subGrid.size()-1].CoarseLimsY[0] = RefYstart[RegionNumber]; this->subGrid[subGrid.size()-1].CoarseLimsY[1] = RefYend[RegionNumber];
		this->subGrid[subGrid.size()-1].CoarseLimsZ[0] = RefZstart[RegionNumber]; this->subGrid[subGrid.size()-1].CoarseLimsZ[1] = RefZend[RegionNumber];

	} else { // Lower subgrids

		this->subGrid[subGrid.size()-1].CoarseLimsX[0] = 2; this->subGrid[subGrid.size()-1].CoarseLimsX[1] = this->XPos.size() - 3;
		this->subGrid[subGrid.size()-1].CoarseLimsY[0] = 2; this->subGrid[subGrid.size()-1].CoarseLimsY[1] = this->YPos.size() - 3;
		this->subGrid[subGrid.size()-1].CoarseLimsZ[0] = 2; this->subGrid[subGrid.size()-1].CoarseLimsZ[1] = this->ZPos.size() - 3;

	}

	// Initialise the subgrid passing necessary data from parent
#if dims == 3
	this->subGrid[subGrid.size()-1].LBM_init_subgrid(	this->XPos[ subGrid[subGrid.size()-1].CoarseLimsX[0] ],
														this->YPos[ subGrid[subGrid.size()-1].CoarseLimsY[0] ],
														this->ZPos[ subGrid[subGrid.size()-1].CoarseLimsZ[0] ],
														this->dx, this->omega	);
#else

    this->subGrid[subGrid.size()-1].LBM_init_subgrid(	this->XPos[ subGrid[subGrid.size()-1].CoarseLimsX[0] ],
														this->YPos[ subGrid[subGrid.size()-1].CoarseLimsY[0] ],
														0.0,
														this->dx, this->omega	);
#endif

	// Add another subgrid beneath the one just created if necessary
	if (this->subGrid[subGrid.size()-1].level < NumLev) {
		this->subGrid[subGrid.size()-1].LBM_addSubGrid(this->subGrid[subGrid.size()-1].region_number);
	}

}

// ***************************************************************************************************
bool GridObj::isEdge (size_t i, size_t j, size_t k, int RegionNumber) {

	bool edgetest;
	
	// Test for edge of grid
	if (i == subGrid[RegionNumber].CoarseLimsX[1] && 
		j == subGrid[RegionNumber].CoarseLimsY[1] && 
		k == subGrid[RegionNumber].CoarseLimsZ[1])
		edgetest = true;

	else edgetest = false;

	return edgetest;
}

// ***************************************************************************************************
bool GridObj::isWithin (size_t i, size_t j, size_t k, int RegionNumber) {

	bool intest;
	
	// Test for edge of grid
	if (	( i >= subGrid[RegionNumber].CoarseLimsX[0] && i <= subGrid[RegionNumber].CoarseLimsX[1] ) &&
			( j >= subGrid[RegionNumber].CoarseLimsY[0] && j <= subGrid[RegionNumber].CoarseLimsX[1] ) &&
			( k >= subGrid[RegionNumber].CoarseLimsZ[0] && k <= subGrid[RegionNumber].CoarseLimsX[1] )
			)
		intest = true;

	else intest = false;

	return intest;
}


// ***************************************************************************************************
// Other methods are in their own files prefixed GridObj_
