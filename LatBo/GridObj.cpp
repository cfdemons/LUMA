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
void GridObj::build_body(int body_type) {

	// Declarations
	std::vector<double> dimensions, centrepoint, start_position, end_position;

	// Increase the iBody array by single IB_body object
	iBody.emplace_back();

	if (body_type == 1) {
		// =========================== Build a rectangle/cuboid ===========================  //

		// Dimensions
		dimensions.push_back(ibb_w);
		dimensions.push_back(ibb_l);
		dimensions.push_back(ibb_d);
		centrepoint.push_back(ibb_x);
		centrepoint.push_back(ibb_y);
		centrepoint.push_back(ibb_z);

		// Build
		iBody[iBody.size()-1].makeBody(dimensions,centrepoint,false,ibb_deform,iBody.size()-1);

	} else if (body_type == 2) {
		// =========================== Build a circle/sphere ===========================  //

		// Dimensions
		centrepoint.push_back(ibb_x);
		centrepoint.push_back(ibb_y);
		centrepoint.push_back(ibb_z);

		// Build
		iBody[iBody.size()-1].makeBody(ibb_r,centrepoint,false,ibb_deform,iBody.size()-1);

	} else if (body_type == 3) {
		// =========================== Build both with custom dimensions ===========================  //

		// Dimensions
		centrepoint.push_back( 2*(double)(b_x-a_x)/5 );
		centrepoint.push_back( 2*(double)(b_y - a_y)/5 );
		centrepoint.push_back(ibb_z);
		dimensions.push_back(ibb_w);
		dimensions.push_back(ibb_l);
		dimensions.push_back(ibb_d);

		// Build circle first
		iBody[iBody.size()-1].makeBody(ibb_r,centrepoint,false,ibb_deform,iBody.size()-1);
		
		// Grow vector
		iBody.emplace_back();

		// Build rectangle
		centrepoint[0] = 3*(double)(b_x - a_x)/5;
		centrepoint[1] = 3*(double)(b_y - a_y)/5;
		centrepoint[2] = ibb_z;
		iBody[iBody.size()-1].makeBody(dimensions,centrepoint,false,ibb_deform,iBody.size()-1);


	} else if (body_type == 4) {
		// =========================== Build flexible filament ===========================  //

		// End points
		start_position.push_back(ibb_start_x);
		start_position.push_back(ibb_start_y);
		start_position.push_back(ibb_start_z);

		// Angle vector
		std::vector<double> angles;
		angles.push_back(ibb_angle_vert);
#if (dims == 3)
		angles.push_back(ibb_angle_horz);
#else
		angles.push_back(0.0);
#endif

		// Pack BC info
		std::vector<int> BCs;
		BCs.push_back(start_BC);
		BCs.push_back(end_BC);

		// Build body
		iBody[iBody.size()-1].makeBody(start_position, ibb_length, angles, BCs, true, true, iBody.size()-1);

		
	} else if (body_type == 5) {
		// =========================== Build array of flexible filaments ===========================  //

		std::vector<double> width_length, centrepoint, angles;
		unsigned int group = 999;


		// If 2D, just build a filament
#if (dims == 2)
		std::cout << "Building plate as a filament as only 2D" << std::endl;
		iBody.resize(iBody.size()-1); // Reduce size of iBody vector
		build_body(4);	// Call add again but as a filament
#else
		// In 3D, build an array of filaments

		// Centre point
		centrepoint.push_back(ibb_x);
		centrepoint.push_back(ibb_y);
		centrepoint.push_back(ibb_z);		

		// Dimensions
		width_length.push_back(ibb_w);
		width_length.push_back(ibb_d);

		// Angles
		angles.push_back(ibb_angle_vert);
		angles.push_back(ibb_angle_horz);

		// Pack BC info
		std::vector<int> BCs;
		BCs.push_back(start_BC);
		BCs.push_back(end_BC);

		// Declarations
		std::vector<double> start_point; start_point.resize(3);

		// Get length projected onto horizontal plane
		double width_h = width_length[0]*cos(angles[0]);
	
		// Define spacing for filaments
		double  lengthwise_spacing;
		// Must use odd number of filaments to ensure ones on edges and in centre
		if (num_markers % 2 == 0) {
			lengthwise_spacing = (width_length[1]) / (num_markers-1);
		} else {
			lengthwise_spacing = (width_length[1]) / num_markers;
		}

		// Loop
		int counter = (int)(width_length[1] / lengthwise_spacing);
		for (int i = 0; i < counter; i++ ) {

			// Get starting point
			start_point[0] = (centrepoint[0] - (width_h*cos(angles[1]))/2 ) + i*( lengthwise_spacing*sin(angles[1]) );
			start_point[1] = centrepoint[1] - (width_length[0] * sin(angles[0]))/2 ;
			start_point[2] = (centrepoint[2] - (width_length[1]*cos(angles[1]))/2 ) + i*( lengthwise_spacing*cos(angles[1]) );

			if (i != 0) {
				// Increase the iBody array by single IB_body object when not on first loop
				iBody.emplace_back();
			}

			// Build filament with fixed group ID making centre filament flexible
			if (i == (counter-1)/2) {
				iBody[iBody.size()-1].makeBody(start_point, width_length[0], angles, BCs, true, true, group);
			} else {
				iBody[iBody.size()-1].makeBody(start_point, width_length[0], angles, BCs, false, true, group);
			}


		}

#endif


	} else if (body_type == 6) {
		// =========================== Build conference case ===========================  //


	}

}

// ***************************************************************************************************
// ***************************************************************************************************
// Other member methods are in their own files prefixed GridObj_
