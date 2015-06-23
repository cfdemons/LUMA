#include "stdafx.h"
#include "GridObj.h"
#include "definitions.h"
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
// Method to call body building routine for IBM
void GridObj::build_body(int type) {

	// Declarations
	std::vector<double> dimensions, centrepoint;

	// Increase the iBody array by single IB_body object
	iBody.emplace_back();

	if (type == 1) {	// Build a rectangle/cuboid

		// Dimensions
		dimensions.push_back(ibb_w);
		dimensions.push_back(ibb_l);
		dimensions.push_back(ibb_d);
		centrepoint.push_back(ibb_x);
		centrepoint.push_back(ibb_y);
		centrepoint.push_back(ibb_z);

		// Build
		iBody[iBody.size()-1].makeBody(dimensions,centrepoint,false);

	} else if (type == 2) { // Build a circle/sphere

		// Dimensions
		centrepoint.push_back(ibb_x);
		centrepoint.push_back(ibb_y);
		centrepoint.push_back(ibb_z);

		// Build
		iBody[iBody.size()-1].makeBody(ibb_r,centrepoint,false);

	} else if (type == 3) {

		// Build both with custom dimensions
		// Dimensions
		centrepoint.push_back( 2*(double)(b_x-a_x)/5 );
		centrepoint.push_back( 2*(double)(b_y - a_y)/5 );
		centrepoint.push_back(ibb_z);
		dimensions.push_back(ibb_w);
		dimensions.push_back(ibb_l);
		dimensions.push_back(ibb_d);

		// Build circle first
		iBody[iBody.size()-1].makeBody(ibb_r,centrepoint,false);
		
		// Grow vector
		iBody.emplace_back();

		// Build rectangle
		centrepoint[0] = 3*(double)(b_x - a_x)/5;
		centrepoint[1] = 3*(double)(b_y - a_y)/5;
		centrepoint[2] = ibb_z;
		iBody[iBody.size()-1].makeBody(dimensions,centrepoint,false);


	}

}

// ***************************************************************************************************
// ***************************************************************************************************
// Other member methods are in their own files prefixed GridObj_
