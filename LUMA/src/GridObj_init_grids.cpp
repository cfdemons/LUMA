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

/* This file contains the routines necessary to initialise the macroscopic quantities and the grids and labels.
*/

#include "../inc/stdafx.h"
#include "../inc/GridObj.h"
#include "../inc/MpiManager.h"
#include "../inc/GridUtils.h"

using namespace std;

// ****************************************************************************
/// \brief	Method to import an input profile from a file.
///
///			Input data may be over- or under-sampled but it must span the 
///			physical dimensions of the inlet otherwise the software does
///			not known how to scale the data to fit. Inlet profile is always
///			assumed to be oriented vertically (y-direction).
void GridObj::LBM_init_getInletProfile() {

	size_t i, j;
	double y, tmp;
	std::vector<double> ybuffer, uxbuffer, uybuffer, uzbuffer;

	// Indicate to log
	*GridUtils::logfile << "Loading inlet profile..." << std::endl;

	// Buffer information from file
	std::ifstream inletfile;
	inletfile.open("./input/inlet_profile.in", std::ios::in);
	if (!inletfile.is_open()) {
		// Error opening file
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "Cannot open inlet profile file named \"inlet_profile.in\". Exiting." << std::endl;
		exit(LUMA_FAILED);

	} else {

		std::string line_in;	// String to store line
		std::istringstream iss;	// Buffer stream

		while( !inletfile.eof() ) {

			// Get line and put in buffer
			std::getline(inletfile,line_in,'\n');
			iss.str(line_in);
			iss.seekg(0); // Reset buffer position to start of buffer

			// Get y position
			iss >> tmp;
			ybuffer.push_back(tmp);

			// Get x velocity
			iss >> tmp;
			uxbuffer.push_back(tmp);

			// Get y velocity
			iss >> tmp;
			uybuffer.push_back(tmp);

			// Get z velocity
			iss >> tmp;
			uzbuffer.push_back(tmp);

		}

	}


	// Resize vectors
	ux_in.resize(M_lim);
	uy_in.resize(M_lim);
	uz_in.resize(M_lim);

	// Loop over site positions (for left hand inlet, y positions)
	for (j = 0; j < M_lim; j++) {

		// Get Y-position
		y = YPos[j];

		// Find values either side of desired
		for (i = 0; i < ybuffer.size(); i++) {

			// Bottom point
			if (i == 0 && ybuffer[i] > y) {
				
				// Extrapolation
				ux_in[j] = -(uxbuffer[i+1] - uxbuffer[i]) + uxbuffer[i];
				uy_in[j] = -(uybuffer[i+1] - uybuffer[i]) + uybuffer[i];
				uz_in[j] = -(uzbuffer[i+1] - uzbuffer[i]) + uzbuffer[i];
				break;

			// Top point
			} else if (i == ybuffer.size()-1 && ybuffer[i] < y) {

				// Extrapolation
				ux_in[j] = (uxbuffer[i] - uxbuffer[i-1]) + uxbuffer[i];
				uy_in[j] = (uybuffer[i] - uybuffer[i-1]) + uybuffer[i];
				uz_in[j] = (uzbuffer[i] - uzbuffer[i-1]) + uzbuffer[i];
				break;


			// Any middle point
			} else if (ybuffer[i] < y && ybuffer[i+1] > y) {

				// Interpolation
				ux_in[j] = ((uxbuffer[i+1] - uxbuffer[i])/(ybuffer[i+1] - ybuffer[i])) * (y-ybuffer[i]) + uxbuffer[i];
				uy_in[j] = ((uybuffer[i+1] - uybuffer[i])/(ybuffer[i+1] - ybuffer[i])) * (y-ybuffer[i]) + uybuffer[i];
				uz_in[j] = ((uzbuffer[i+1] - uzbuffer[i])/(ybuffer[i+1] - ybuffer[i])) * (y-ybuffer[i]) + uzbuffer[i];
				break;

			} else if (ybuffer[i] == y ) {

				// Copy
				ux_in[j] = uxbuffer[i];
				uy_in[j] = uybuffer[i];
				uz_in[j] = uzbuffer[i];
				break;			
			
			} else {

				continue;

			}
				
				
		}

	}

	if (ux_in.size() == 0 || uy_in.size() == 0 || uz_in.size() == 0) {
		// No data read in
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "Failed to read in inlet profile data. Exiting." << std::endl;
		exit(LUMA_FAILED);
	}


}

// ****************************************************************************
/// \brief	Method to initialise the lattice velocity.
///
///			Unless the L_NO_FLOW macro is defined, the initial velocity 
///			everywhere will be set to the values specified in the definitions 
///			file.
void GridObj::LBM_initVelocity ( ) {


	// Loop over grid
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {
				
#ifdef L_NO_FLOW				
				for (size_t d = 0; d < L_dims; d++) {

					// No flow case
					u(i,j,k,d,M_lim,K_lim,L_dims) = 0.0;
				}

#else
				/* Input velocity is specified by individual vectors for x, y and z which 
				 * have either been read in from an input file or defined by an expression 
				 * given in the definitions.
				 */
				u(i,j,k,0,M_lim,K_lim,L_dims) = L_u_0x;
				u(i,j,k,1,M_lim,K_lim,L_dims) = L_u_0y;
#if (L_dims == 3)
				u(i,j,k,2,M_lim,K_lim,L_dims) = L_u_0z;
#endif


#endif
			}
		}
	}

	// Perform solid site reset of velocity for any solids
	bc_solidSiteReset();

}


// ****************************************************************************
/// \brief	Method to initialise the lattice density.
void GridObj::LBM_initRho ( ) {

	// Loop over grid
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {		
			for (int k = 0; k < K_lim; k++) {

				// Uniform Density
				rho(i,j,k,M_lim,K_lim) = L_rho_in;
			}
		}
	}

}

// ****************************************************************************
/// \brief	Wrapper to initialise all L0 lattice quantities.
///
///			This method wraps the MPI-specific version. It is called by the 
///			serial build and sets the MPI-specific arguments to default values
///			before calling the full initialiser.
void GridObj::LBM_initGrid( ) {

	// Set default value for the following MPI-specific settings
	std::vector<int> local_size;
	std::vector< std::vector<int> > global_edge_ind;
	std::vector< std::vector<double> > global_edge_pos;

	// Set local size to total grid size
	local_size.push_back(L_N);
	local_size.push_back(L_M);
	local_size.push_back(L_K);

	// Set global edge indices and position indices arbitrarily to zero as they won't be accessed anyway
	global_edge_ind.resize(1, std::vector<int>(1) );
	global_edge_ind[0][0] = 0;	
	global_edge_pos.resize(1, std::vector<double>(1) );
	global_edge_pos[0][0] = 0.0;

	// Call MPI initialiser wiht these default options
	LBM_initGrid(local_size, global_edge_ind, global_edge_pos);

}


// ****************************************************************************
/// \brief	Method to initialise all L0 lattice quantities.
/// \param	local_size	local grid size on this rank including halo.
/// \param	global_edge_ind	global indices of the rank edges.
/// \param	global_edge_pos	global positions of the rank edges.
void GridObj::LBM_initGrid( std::vector<int> local_size, 
							std::vector< std::vector<int> > global_edge_ind, 
							std::vector< std::vector<double> > global_edge_pos ) {
	
	// Store physical spacing
	// Global dimensions
	double Lx = L_b_x - L_a_x;
	double Ly = L_b_y - L_a_y;
	double Lz = L_b_z - L_a_z;
	dx = 2 * (Lx / (2 * static_cast<double>(L_N)));
	dy = 2 * (Ly / (2 * static_cast<double>(L_M)));
	dz = 2 * (Lz / (2 * static_cast<double>(L_K)));

	// Physical time step = physical grid spacing
	dt = dx;
	


	////////////////////////////
	// Check input parameters //
	////////////////////////////

#if (L_dims == 3)
	// Check that lattice volumes are cubes in 3D
	if ( (Lx/L_N) != (Ly/L_M) || (Lx/L_N) != (Lz/L_K) ) {
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "Need to have lattice volumes which are cubes -- either change L_N/L_M/L_K or change domain dimensions. Exiting." << std::endl;
		exit(LUMA_FAILED);
	}
	
#else
	// 2D so need square lattice cells
	if ( (Lx/L_N) != (Ly/L_M) ) {
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "Need to have lattice cells which are squares -- either change L_N/L_M or change domain dimensions. Exiting." << std::endl;
		exit(LUMA_FAILED);
	}

#endif

    // Checks to make sure grid size is suitable for refinement
	if (L_NumLev != 0) {

#if (L_dims == 3)
		for (int reg = 0; reg < L_NumReg; reg++) {
			// Check grid is big enough to allow embedded refinement of factor 2
			if (	(
					RefXend[level][reg]-RefXstart[level][reg]+1 < 3 || 
					RefYend[level][reg]-RefYstart[level][reg]+1 < 3 || 
					RefZend[level][reg]-RefZstart[level][reg]+1 < 3
					) || (
					(
					RefXend[level][reg]-RefXstart[level][reg]+1 == 3 ||
					RefYend[level][reg]-RefYstart[level][reg]+1 == 3 ||
					RefZend[level][reg]-RefZstart[level][reg]+1 == 3
					) && L_NumLev > 1 )
				) {
					std::cout << "Error: See Log File" << std::endl;
					*GridUtils::logfile << "Refined region is too small to support refinement. Exiting." << std::endl;
					exit(LUMA_FAILED);
			}
		}
#else
		for (int reg = 0; reg < L_NumReg; reg++) {
			// Check grid is big enough to allow embedded refinement of factor 2
			if (	(
					RefXend[level][reg]-RefXstart[level][reg]+1 < 3 || 
					RefYend[level][reg]-RefYstart[level][reg]+1 < 3
					) || (
					(
					RefXend[level][reg]-RefXstart[level][reg]+1 == 3 ||
					RefYend[level][reg]-RefYstart[level][reg]+1 == 3
					) && L_NumLev > 1 )
				) {
					std::cout << "Error: See Log File" << std::endl;
					*GridUtils::logfile << "Refined region is too small to support refinement. Exiting." << std::endl;
					exit(LUMA_FAILED);
			}
		}
#endif
	}


    ///////////////////
	// Checks passed //
	///////////////////

	// Get local grid sizes (includes overlap)
	N_lim = local_size[0];
	M_lim = local_size[1];
#if (L_dims == 3)
	K_lim = local_size[2];
#else
	K_lim = 1;
#endif
	

	// NODE NUMBERS on L0
	/* When using MPI:
	 * Node numbers should be specified in the global system.
	 * The overlapping sites have the same index as an index
	 * on its overlapping grid. The setup assumes periodicity
	 * on all edges, even where the edge of the domain is a 
	 * boundary and so halo cells exist on all edges.
	 */
#ifdef L_BUILD_FOR_MPI	

	// Build index vectors
	XInd = GridUtils::onespace( (int)global_edge_ind[0][MpiManager::my_rank], (int)global_edge_ind[1][MpiManager::my_rank] - 1 );
	YInd = GridUtils::onespace( (int)global_edge_ind[2][MpiManager::my_rank], (int)global_edge_ind[3][MpiManager::my_rank] - 1 );
	ZInd = GridUtils::onespace( (int)global_edge_ind[4][MpiManager::my_rank], (int)global_edge_ind[5][MpiManager::my_rank] - 1 );

	// Add overlap indices to both ends of the vector taking into account periodicity
	XInd.insert( XInd.begin(), (XInd[0]-1 + L_N) % L_N ); XInd.insert( XInd.end(), (XInd[XInd.size()-1]+1 + L_N) % L_N );
	YInd.insert( YInd.begin(), (YInd[0]-1 + L_M) % L_M ); YInd.insert( YInd.end(), (YInd[YInd.size()-1]+1 + L_M) % L_M );
#if (L_dims == 3)
	ZInd.insert( ZInd.begin(), (ZInd[0]-1 + L_K) % L_K ); ZInd.insert( ZInd.end(), (ZInd[ZInd.size()-1]+1 + L_K) % L_K );
#endif

#else
	// When not builiding for MPI indices are straightforward
	XInd = GridUtils::onespace( 0, L_N-1 );
	YInd = GridUtils::onespace( 0, L_M-1 );
	ZInd = GridUtils::onespace( 0, L_K-1 );
#endif



	// L0 lattice site POSITION VECTORS
	/* When using MPI:
	 * As with the indices, the overlap is assume periodic
	 * in all directions.
	 */

#ifdef L_BUILD_FOR_MPI

	// Create position vectors excluding overlap
	XPos = GridUtils::linspace( global_edge_pos[0][MpiManager::my_rank] + dx/2, global_edge_pos[1][MpiManager::my_rank] - dx/2, N_lim-2 );
	YPos = GridUtils::linspace( global_edge_pos[2][MpiManager::my_rank] + dy/2, global_edge_pos[3][MpiManager::my_rank] - dy/2, M_lim-2 );
	ZPos = GridUtils::linspace( global_edge_pos[4][MpiManager::my_rank] + dz/2, global_edge_pos[5][MpiManager::my_rank] - dz/2, K_lim-2 );

	// Add overlap sites taking into account periodicity
	XPos.insert( XPos.begin(), fmod(XPos[0]-dx + Lx, Lx) ); XPos.insert( XPos.end(), fmod(XPos[XPos.size()-1]+dx + Lx, Lx) );
	YPos.insert( YPos.begin(), fmod(YPos[0]-dy + Ly, Ly) ); YPos.insert( YPos.end(), fmod(YPos[YPos.size()-1]+dy + Ly, Ly) );
#if (L_dims == 3)
	ZPos.insert( ZPos.begin(), fmod(ZPos[0]-dz + Lz, Lz) ); ZPos.insert( ZPos.end(), fmod(ZPos[ZPos.size()-1]+dz + Lz, Lz) );
#endif

	// Update the sender/recv layer positions in the MpiManager
	MpiManager* mpim = MpiManager::getInstance();

	// X
	mpim->sender_layer_pos.X[0] = XPos[1] - dx/2;					mpim->sender_layer_pos.X[1] = XPos[1] + dx/2;
	mpim->sender_layer_pos.X[2] = XPos[local_size[0] - 2] - dx/2;	mpim->sender_layer_pos.X[3] = XPos[local_size[0] - 2] + dx/2;
	mpim->recv_layer_pos.X[0]	= XPos[0] - dx/2;					mpim->recv_layer_pos.X[1]	= XPos[0] + dx/2;
	mpim->recv_layer_pos.X[2]	= XPos[local_size[0] - 1] - dx/2;	mpim->recv_layer_pos.X[3]	= XPos[local_size[0] - 1] + dx/2;
	// Y
	mpim->sender_layer_pos.Y[0] = YPos[1] - dy/2;					mpim->sender_layer_pos.Y[1] = YPos[1] + dy/2;
	mpim->sender_layer_pos.Y[2] = YPos[local_size[1] - 2] - dy/2;	mpim->sender_layer_pos.Y[3] = YPos[local_size[1] - 2] + dy/2;
	mpim->recv_layer_pos.Y[0]	= YPos[0] - dy/2;					mpim->recv_layer_pos.Y[1]	= YPos[0] + dy/2;
	mpim->recv_layer_pos.Y[2]	= YPos[local_size[1] - 1] - dy/2;	mpim->recv_layer_pos.Y[3]	= YPos[local_size[1] - 1] + dy/2;
	// Z
#if (L_dims == 3)
	mpim->sender_layer_pos.Z[0] = ZPos[1] - dz/2;					mpim->sender_layer_pos.Z[1] = ZPos[1] + dz/2;
	mpim->sender_layer_pos.Z[2] = ZPos[local_size[2] - 2] - dz/2;	mpim->sender_layer_pos.Z[3] = ZPos[local_size[2] - 2] + dz/2;
	mpim->recv_layer_pos.Z[0]	= ZPos[0] - dz/2;					mpim->recv_layer_pos.Z[1]	= ZPos[0] + dz/2;
	mpim->recv_layer_pos.Z[2]	= ZPos[local_size[2] - 1] - dz/2;	mpim->recv_layer_pos.Z[3]	= ZPos[local_size[2] - 1] + dz/2;
#endif

#ifdef L_MPI_VERBOSE
	*mpim->logout << "X sender layers are: " << mpim->sender_layer_pos.X[0] << " -- " << mpim->sender_layer_pos.X[1] << " (min edge) , " << mpim->sender_layer_pos.X[2] << " -- " << mpim->sender_layer_pos.X[3] << " (max edge)" << std::endl;
	*mpim->logout << "X recv layers are: " << mpim->recv_layer_pos.X[0] << " -- " << mpim->recv_layer_pos.X[1] << " (min edge) , " << mpim->recv_layer_pos.X[2] << " -- " << mpim->recv_layer_pos.X[3] << " (max edge)" << std::endl;

	*mpim->logout << "Y sender layers are: " << mpim->sender_layer_pos.Y[0] << " -- " << mpim->sender_layer_pos.Y[1] << " (min edge) , " << mpim->sender_layer_pos.Y[2] << " -- " << mpim->sender_layer_pos.Y[3] << " (max edge)" << std::endl;
	*mpim->logout << "Y recv layers are: " << mpim->recv_layer_pos.Y[0] << " -- " << mpim->recv_layer_pos.Y[1] << " (min edge) , " << mpim->recv_layer_pos.Y[2] << " -- " << mpim->recv_layer_pos.Y[3] << " (max edge)" << std::endl;

#if (L_dims == 3)
	*mpim->logout << "Z sender layers are: " << mpim->sender_layer_pos.Z[0] << " -- " << mpim->sender_layer_pos.Z[1] << " (min edge) , " << mpim->sender_layer_pos.Z[2] << " -- " << mpim->sender_layer_pos.Z[3] << " (max edge)" << std::endl;
	*mpim->logout << "Z recv layers are: " << mpim->recv_layer_pos.Z[0] << " -- " << mpim->recv_layer_pos.Z[1] << " (min edge) , " << mpim->recv_layer_pos.Z[2] << " -- " << mpim->recv_layer_pos.Z[3] << " (max edge)" << std::endl;
#endif
#endif

#else
	// When not builiding for MPI positions are straightforward
	XPos = GridUtils::linspace( L_a_x + dx/2, L_b_x - dx/2, L_N );
	YPos = GridUtils::linspace( L_a_y + dy/2, L_b_y - dy/2, L_M );
	ZPos = GridUtils::linspace( L_a_z + dz/2, L_b_z - dz/2, L_K );
#endif

	

	// Define TYPING MATRICES
	LatTyp.resize( N_lim*M_lim*K_lim );

	// Label as coarse site
	std::fill(LatTyp.begin(), LatTyp.end(), eFluid);

	// Add boundary-specific labels
	LBM_initBoundLab();



	// Initialise L0 MACROSCOPIC quantities

	// Get the inlet profile data
#ifdef L_USE_INLET_PROFILE

	LBM_init_getInletProfile();
	
#endif

	// Velocity field
	u.resize( N_lim*M_lim*K_lim*L_dims );
	LBM_initVelocity();
	
	// Density field
	rho.resize( N_lim*M_lim*K_lim );
	LBM_initRho();

	// Cartesian force vector
	force_xyz.resize(N_lim*M_lim*K_lim*L_dims, 0.0);

	// Lattice force vector
	force_i.resize(N_lim*M_lim*K_lim*L_nVels, 0.0);

	// Time averaged quantities
	rho_timeav.resize(N_lim*M_lim*K_lim, 0.0);
	ui_timeav.resize(N_lim*M_lim*K_lim*L_dims, 0.0);
	uiuj_timeav.resize(N_lim*M_lim*K_lim*(3*L_dims-3), 0.0);


	// Initialise L0 POPULATION matrices (f, feq)
	f.resize( N_lim*M_lim*K_lim*L_nVels );
	feq.resize( N_lim*M_lim*K_lim*L_nVels );



	// Loop over grid
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {
				for (int v = 0; v < L_nVels; v++) {

					// Initialise f to feq
					f(i,j,k,v,M_lim,K_lim,L_nVels) = LBM_collide(i, j, k, v);

				}
			}
		}
	}
	feq = f; // Make feq = feq too


	// Initialise OTHER parameters
	// Compute kinematic viscosity based on target Reynolds number
#if defined L_IBM_ON && defined L_INSERT_CIRCLE_SPHERE
	// If IBM circle use diameter (in lattice units i.e. rescale wrt to physical spacing)
	nu = (L_ibb_r*2 / dx) * L_u_ref / L_Re;
#elif defined L_IBM_ON && defined L_INSERT_RECTANGLE_CUBOID
	// If IBM rectangle use y-dimension (in lattice units)
	nu = (L_ibb_l / dx) * L_u_ref / L_Re;
#elif defined L_SOLID_BLOCK_ON
	// Use block length (scaled back to L0 units)
	nu = ((L_block_x_max - L_block_x_min) / pow(2,L_block_on_grid_lev)) * L_u_ref / L_Re;
#elif defined L_IBM_ON && defined L_IBB_FROM_FILE
	// If IBM object read from file then use scale length as reference
	nu = (L_ibb_length_ref / dx) * L_u_ref / L_Re;
#elif defined L_SOLID_FROM_FILE
	// Use object length (scaled back to L0 units)
	nu = (L_object_length_ref / pow(2,L_object_on_grid_lev)) * L_u_ref / L_Re;
#elif defined L_BFL_ON
	// Use bfl body length (scaled back to L0 units)
	nu = (L_bfl_length_ref / pow(2,L_bfl_on_grid_lev)) * L_u_ref / L_Re;
#else
	// If no object then use domain height (in lattice units)
	nu = (L_M - L_wall_thickness_bottom - L_wall_thickness_top) * L_u_ref / L_Re;	// Based on actual width of channel (in lattice units)
#endif

	// Relaxation frequency on L0
	// Assign relaxation frequency using lattice viscosity
	omega = 1 / ( (nu / pow(cs,2)) + .5 );

	/* Above is valid for L0 only when dx = 1 -- general expression is:
	 * omega = 1 / ( ( (nu * dt) / (pow(cs,2)*pow(dx,2)) ) + .5 );
	 */

}
	


// ****************************************************************************
/// \brief	Method to initialise all sub-grid quantities.
/// \param	pGrid	reference to parent grid.
void GridObj::LBM_initSubGrid (GridObj& pGrid) {
	
	// Declarations
	int IndXstart, IndYstart, IndZstart = 0;
	int offset = (int)pow(2,pGrid.level);	// How many sites thick the recv data region is

	/* MPI specific setup:
	 * 1. Store coarse grid refinement limits;
	 * 2. Get node numbering ends in global system
	 *
	 * These are stored as local indices as they are used to map between the 
	 * fine and coarse grid cells during multi-grid operations. Therefore, we 
	 * must only store local values relevant to the grid on the rank and not the
	 * refined region as a whole or mapping will not be correct.
	 *
	 * When not using MPI, these can be read straight from the definitions file and 
	 * converted to local coordinates corresponding to the parent grid.
	 * However, when using MPI, the edges of the refined grid might not be on this 
	 * rank so we must round the coarse limits to the edge of the parent grid so 
	 * the correct offset is supplied to the mapping routine.
	 *
	 * When dealing with sub-grids embedded in walls, it is more complicated. If the 
	 * sub-grid starts on a max receiver layer which is also a periodic overlap then
	 * we need to make sure the limits are set properly. Likewise if the sub-grid ends
	 * on a min receiver layer which is also periodic.
	 */

	// If region is not contained on a single rank adjust limits accordingly:
	
	// X //

	// Start Limit
	// Find the local index of the refinement limits if they are on the rank at all
	auto found_x = std::find(pGrid.XInd.begin(), pGrid.XInd.end(), RefXstart[pGrid.level][region_number]);
	if (found_x != pGrid.XInd.end()) {	// Starts on this rank
		CoarseLimsX[0] = static_cast<int>(found_x - pGrid.XInd.begin());	// Store local index as the limit
		// Start index is simply zero
		IndXstart = 0;

	// Starts on some rank to the left of this one
	} else if ( (int)RefXstart[pGrid.level][region_number] < pGrid.XInd[offset] - offset ) {
		// Set limit to start edge of the rank
		CoarseLimsX[0] = 0;
		// Compute starting index of sub-grid on this rank (take into account rest of grid somewhere to left)
		IndXstart = (pGrid.XInd[CoarseLimsX[0]] - RefXstart[pGrid.level][region_number]) * 2;

	}

	// End Limit
	found_x = std::find(pGrid.XInd.begin(), pGrid.XInd.end(), RefXend[pGrid.level][region_number]);
	if (found_x != pGrid.XInd.end()) {	// Ends on this rank
		CoarseLimsX[1] = static_cast<int>(found_x - pGrid.XInd.begin());

	// End on some rank to the right of this one
	} else if ( (int)RefXend[pGrid.level][region_number] > pGrid.XInd[pGrid.XInd.size() - 1 - offset] + offset ) {
		// Set grid limits to end edge of the rank
		CoarseLimsX[1] = static_cast<int>(pGrid.XInd.size()) - 1;

	// Else the end must be on a rank to the left and hence grid wraps periodically so set end to right-hand edge
	} else if ( (int)RefXend[pGrid.level][region_number] < pGrid.XInd[offset] - offset ) {
		// Set grid limits to end edge of the rank
		CoarseLimsX[1] = static_cast<int>(pGrid.XInd.size()) - 1;

	}


	// Y //
	auto found_y = std::find(pGrid.YInd.begin(), pGrid.YInd.end(), RefYstart[pGrid.level][region_number]);
	if (found_y != pGrid.YInd.end()) {
		CoarseLimsY[0] = static_cast<int>(found_y - pGrid.YInd.begin());
		IndYstart = 0;

	} else if ( (int)RefYstart[pGrid.level][region_number] < pGrid.YInd[offset] - offset ) {
		CoarseLimsY[0] = 0;
		IndYstart = (pGrid.YInd[CoarseLimsY[0]] - RefYstart[pGrid.level][region_number]) * 2;
	}

	found_y = std::find(pGrid.YInd.begin(), pGrid.YInd.end(), RefYend[pGrid.level][region_number]);
	if (found_y != pGrid.YInd.end()) {
		CoarseLimsY[1] = static_cast<int>(found_y - pGrid.YInd.begin());

	} else if ( (int)RefYend[pGrid.level][region_number] > pGrid.YInd[pGrid.YInd.size() - 1 - offset] + offset ) {
		CoarseLimsY[1] = static_cast<int>(pGrid.YInd.size()) - 1;

	} else if ( (int)RefYend[pGrid.level][region_number] < pGrid.YInd[offset] - offset ) {
		CoarseLimsY[1] = static_cast<int>(pGrid.YInd.size()) - 1;

	}


#if (L_dims == 3)
	// Z //
	auto found_z = std::find(pGrid.ZInd.begin(), pGrid.ZInd.end(), RefZstart[pGrid.level][region_number]);
	if (found_z != pGrid.ZInd.end()) {
		CoarseLimsZ[0] = static_cast<int>(found_z - pGrid.ZInd.begin());
		IndZstart = 0;

	} else if ( (int)RefZstart[pGrid.level][region_number] < pGrid.ZInd[offset] - offset ) {
		CoarseLimsZ[0] = 0;
		IndZstart = (pGrid.ZInd[CoarseLimsZ[0]] - RefZstart[pGrid.level][region_number]) * 2;
	}

	found_z = std::find(pGrid.ZInd.begin(), pGrid.ZInd.end(), RefZend[pGrid.level][region_number]);
	if (found_z != pGrid.ZInd.end()) {
		CoarseLimsZ[1] = static_cast<int>(found_z - pGrid.ZInd.begin());

	} else if ( (int)RefZend[pGrid.level][region_number] > pGrid.ZInd[pGrid.ZInd.size() - 1 - offset] + offset ) {
		CoarseLimsZ[1] = static_cast<int>(pGrid.ZInd.size()) - 1;

	} else if ( (int)RefZend[pGrid.level][region_number] < pGrid.ZInd[offset] - offset ) {
		CoarseLimsZ[1] = static_cast<int>(pGrid.ZInd.size()) - 1;

	}
#else
	// Reset the refined region z-limits if only 2D
	for (int i = 0; i < 2; i++) {
		CoarseLimsZ[i] = 0;
	}
#endif


	*GridUtils::logfile << "Local Coarse Lims are " << 
		CoarseLimsX[0] << "-" << CoarseLimsX[1] << ", " << 
		CoarseLimsY[0] << "-" << CoarseLimsY[1] << ", " <<
		CoarseLimsZ[0] << "-" << CoarseLimsZ[1] << std::endl;

	/* If sub-grid wraps periodically we do not support it wrapping back round to the same rank
	 * again as this will confuse the mapping function so we error here. */

	if (
		(CoarseLimsX[1] < CoarseLimsX[0]) ||
		(CoarseLimsY[1] < CoarseLimsY[0]) ||
		(CoarseLimsZ[1] < CoarseLimsZ[0])
		) {
		
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "Refined region wraps periodically which is not supported. Exiting." << std::endl;
		exit(LUMA_FAILED);
	}	


	
	/* Note that the CoarseLims are local values used to identify how much and which part of the 
	 * parent grid is covered by the child sub-grid and might not conincide with edges of global 
	 * refined patch defined by the input file if broken up across ranks. So we get the offset 
	 * in position from the local refined patch edges specific to this rank. */

	// Get positons of local refined patch edges
	double posOffsetX[2] = { pGrid.XPos[ CoarseLimsX[0] ], pGrid.XPos[ CoarseLimsX[1] ] };
	double posOffsetY[2] = { pGrid.YPos[ CoarseLimsY[0] ], pGrid.YPos[ CoarseLimsY[1] ] };
#if (L_dims == 3)
	double posOffsetZ[2] = { pGrid.ZPos[ CoarseLimsZ[0] ], pGrid.ZPos[ CoarseLimsZ[0] ] };
#else
	double posOffsetZ[2] = { 0.0, 0.0 };
#endif

	// Get local grid size of the sub grid based on limits
	int local_size[L_dims] = {
		(int)((CoarseLimsX[1] - CoarseLimsX[0] + .5)*2) + 1,
		(int)((CoarseLimsY[1] - CoarseLimsY[0] + .5)*2) + 1
#if (L_dims == 3)
		, (int)((CoarseLimsZ[1] - CoarseLimsZ[0] + .5)*2) + 1
#endif
	};

	// Generate NODE NUMBERS
	XInd = GridUtils::onespace( IndXstart, IndXstart + local_size[0] - 1 );
	YInd = GridUtils::onespace( IndYstart, IndYstart + local_size[1] - 1 );
#if (L_dims == 3)
	ZInd = GridUtils::onespace( IndZstart, IndZstart + local_size[2] - 1 );
#else
	ZInd.insert(ZInd.begin(), 0); // Default for 2D
#endif



	// Generate POSITION VECTORS of nodes
	// Define spacing
	dx = pGrid.dx/2;
	dy = dx;
	dz = dx;

	// Populate the position vectors
	XPos = GridUtils::linspace(posOffsetX[0] - dx / 2, (posOffsetX[0] - dx / 2) + (XInd.size() - 1) * dx, static_cast<int>(XInd.size()) );
	YPos = GridUtils::linspace(posOffsetY[0] - dy / 2, (posOffsetY[0] - dy / 2) + (YInd.size() - 1) * dy, static_cast<int>(YInd.size()) );
#if L_dims == 3
	ZPos = GridUtils::linspace(posOffsetZ[0] - dz / 2, (posOffsetZ[0] - dz / 2) + (ZInd.size() - 1) * dz, static_cast<int>(ZInd.size()) );
#else
	ZPos.insert( ZPos.begin(), 1 ); // 2D default
#endif

	
	// Generate TYPING MATRICES
	
	// Get local grid sizes
	N_lim = static_cast<int>(XInd.size());
	M_lim = static_cast<int>(YInd.size());
	K_lim = static_cast<int>(ZInd.size());
	
	// Resize
	LatTyp.resize( N_lim * M_lim * K_lim );


	// Default labelling of coarse
	std::fill(LatTyp.begin(), LatTyp.end(), eFluid);
	
	// Call refined labelling routine passing parent grid
	LBM_initRefinedLab(pGrid);

	
	
	// Assign MACROSCOPIC quantities

#ifdef L_USE_INLET_PROFILE

	LBM_init_getInletProfile();
	
#endif

	// Velocity
	u.resize(N_lim * M_lim * K_lim * L_dims);
	LBM_initVelocity( );

	// Density
	rho.resize(N_lim * M_lim * K_lim);
	LBM_initRho( );

	// Cartesian force vector
	force_xyz.resize(N_lim * M_lim * K_lim * L_dims, 0.0);

	// Lattice force vector
	force_i.resize(N_lim * M_lim * K_lim * L_nVels, 0.0);

	// Time averaged quantities
	rho_timeav.resize(N_lim*M_lim*K_lim, 0.0);
	ui_timeav.resize(N_lim*M_lim*K_lim*L_dims, 0.0);
	uiuj_timeav.resize(N_lim*M_lim*K_lim*(3*L_dims-3), 0.0);


	// Generate POPULATION MATRICES for lower levels
	// Resize
	f.resize(N_lim * M_lim * K_lim * L_nVels);
	feq.resize(N_lim * M_lim * K_lim * L_nVels);
	

	// Loop over grid
	for (int i = 0; i != N_lim; ++i) {
		for (int j = 0; j != M_lim; ++j) {
			for (int k = 0; k != K_lim; ++k) {
				for (int v = 0; v < L_nVels; v++) {
					
					// Initialise f to feq
					f(i,j,k,v,M_lim,K_lim,L_nVels) = LBM_collide(i, j, k, v);

				}
			}
		}
	}
	feq = f; // Set feq to feq

	// Compute relaxation time from coarser level assume refinement by factor of 2
	omega = 1 / ( ( (1/pGrid.omega - .5) *2) + .5);

	// Lattice viscosity is constant across subgrids
	nu = pGrid.nu;

}

// ****************************************************************************
/// \brief	Method to initialise label-based solids
void GridObj::LBM_initSolidLab() {

#ifdef L_SOLID_BLOCK_ON
	// Return if not to be put on the current grid
	if (L_block_on_grid_lev != level || L_block_on_grid_reg != region_number) return;

	// Declarations
	int i, j, k, local_i, local_j, local_k;

	/* Check solid block contained on the global grid specified in the definitions file.
	 * If not then exit as user has specifed block outside the grid on which they want it 
	 * placed. */

	// Get global grid sizes
	int Ng_lim, Mg_lim, Kg_lim;
	if (level != 0) {
		Ng_lim = (RefXend[level-1][region_number] - RefXstart[level-1][region_number] + 1) * 2;
		Mg_lim = (RefYend[level-1][region_number] - RefYstart[level-1][region_number] + 1) * 2;
		Kg_lim = (RefZend[level-1][region_number] - RefZstart[level-1][region_number] + 1) * 2;
	} else {
		Ng_lim = L_N; Mg_lim = L_M; Kg_lim = L_K;
	}


	// Check block placement -- must not be on TL (last two sites) if on a level other than 0
	if	(
		(
		
		(level == 0) && 
		
		(L_block_x_max > Ng_lim - 1 || L_block_x_min < 0 || L_block_y_max > Mg_lim - 1 || L_block_y_min < 0 
#if (L_dims == 3)
		|| L_block_z_max > Kg_lim - 1 || L_block_z_min < 0
#endif
		)

		) || (

		(level != 0) && 

		(L_block_x_max >= Ng_lim - 2 || L_block_x_min <= 1 || L_block_y_max >= Mg_lim - 2 || L_block_y_min <= 1 
#if (L_dims == 3)
		|| L_block_z_max >= Kg_lim - 2 || L_block_z_min <= 1
#endif
		)
	
		)
		){

		// Block outside grid
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "Block is placed outside or on the TL of the selected grid. Exiting." << std::endl;
		exit(LUMA_FAILED);
	}


	// Loop over object definition in global indices
	for (i = L_block_x_min; i <= L_block_x_max; i++) {
		for (j = L_block_y_min; j <= L_block_y_max; j++) {
			for (k = L_block_z_min; k <= L_block_z_max; k++)
			{

				// Only label if the site is on current rank
				if ( GridUtils::isOnThisRank(i,j,k,*this) ) {
						
						// Map global indices to local indices
						local_i = i - XInd[1] + 1;
						local_j = j - YInd[1] + 1;
#if (L_dims == 3)
						local_k = k - ZInd[1] + 1;
#else
						local_k = k;
#endif

						LatTyp(local_i,local_j,local_k,M_lim,K_lim) = eSolid;
				}

			}
		}
	}

#endif	// L_SOLID_BLOCK_ON

}

// ****************************************************************************
/// \brief	Method to initialise wall and object labels on L0.
///
///			The virtual wind tunnel definitions are implemented by this method.
void GridObj::LBM_initBoundLab ( ) {

#if defined L_WALLS_ON || defined L_INLET_ON || defined L_OUTLET_ON
	int i, j, k;
#endif

	// Try to add the solid block
	LBM_initSolidLab();

#ifdef L_INLET_ON
	// Left hand face only

	// Check for potential singularity in BC
	if (L_u_max == 1 || L_u_ref == 1) {
		// Singularity so exit
		std::cout << "Error: See Log File" << std::endl;
		*GridUtils::logfile << "Inlet BC fails with L_u_0x = 1, choose something else. Exiting." << std::endl;
		exit(LUMA_FAILED);
	}

	// Search index vector to see if left hand wall on this rank
	for (i = 0; i < N_lim; i++ ) {
		if (XInd[i] == 0) {		// Wall found

			// Label inlet
			for (j = 0; j < M_lim; j++) {
				for (k = 0; k < K_lim; k++) {

					// Inlet site
					LatTyp(i,j,k,M_lim,K_lim) = eInlet;

				}
			}
			break;

		}
	}
#endif

#ifdef L_OUTLET_ON
	// Right hand face only

	// Search index vector to see if right hand wall on this rank
	for (i = 0; i < N_lim; i++ ) {
		if (XInd[i] == L_N - 1) {		// Wall found

			// Label outlet
			for (j = 0; j < M_lim; j++) {
				for (k = 0; k < K_lim; k++) {

#ifdef L_FREESTREAM_TUNNEL
					LatTyp(i,j,k,M_lim,K_lim) = eInlet;
#else
					LatTyp(i,j,k,M_lim,K_lim) = eOutlet;
#endif

				}
			}
			break;

		}
	}
#endif

#if ((defined L_WALLS_ON && !defined L_WALLS_ON_2D) || defined L_FREESTREAM_TUNNEL) && (L_dims == 3)

	// Search index vector to see if FRONT wall on this rank
	for (k = 0; k < K_lim; k++ ) {
		if (ZInd[k] < L_wall_thickness_front) {		// Wall found

			// Label wall
			for (i = 0; i < N_lim; i++) {
				for (j = 0; j < M_lim; j++) {

#if (defined L_UPSTREAM_TUNNEL || defined L_FREESTREAM_TUNNEL)
					LatTyp(i,j,k,M_lim,K_lim) = eInlet;
#else
					LatTyp(i,j,k,M_lim,K_lim) = eSolid;
#endif

				}
			}

		}
	}

	// Search index vector to see if BACK wall on this rank
	for (k = 0; k < K_lim; k++ ) {
		if (ZInd[k] > L_K-1 - L_wall_thickness_back) {		// Wall found

			// Label wall
			for (i = 0; i < N_lim; i++) {
				for (j = 0; j < M_lim; j++) {

#if (defined L_UPSTREAM_TUNNEL || defined L_FREESTREAM_TUNNEL)
					LatTyp(i,j,k,M_lim,K_lim) = eInlet;
#else
					LatTyp(i,j,k,M_lim,K_lim) = eSolid;
#endif

				}
			}

		}
	}

#endif


	// Search index vector to see if BOTTOM wall on this rank
	for (j = 0; j < M_lim; j++ ) {
		if (YInd[j] < L_wall_thickness_bottom) {		// Wall found

			// Label wall
			for (i = 0; i < N_lim; i++) {
				for (k = 0; k < K_lim; k++) {

#ifdef L_WALLS_ON
					LatTyp(i,j,k,M_lim,K_lim) = eSolid;
#elif (defined L_UPSTREAM_TUNNEL || defined L_FREESTREAM_TUNNEL)
					LatTyp(i,j,k,M_lim,K_lim) = eInlet;	// Label as inlet (for rolling road -- velocity BC)
#endif

				}
			}
		}
	}



	// Search index vector to see if TOP wall on this rank
	for (j = 0; j < M_lim; j++ ) {
		if (YInd[j] > L_M-1 - L_wall_thickness_top) {		// Wall found

			// Label wall
			for (i = 0; i < N_lim; i++) {
				for (k = 0; k < K_lim; k++) {

#ifdef L_WALLS_ON
					LatTyp(i,j,k,M_lim,K_lim) = eSolid;
#elif (defined L_UPSTREAM_TUNNEL || defined L_FREESTREAM_TUNNEL)
					LatTyp(i,j,k,M_lim,K_lim) = eInlet;	// Label as free-stream
#endif

				}
			}
		}
	}

}

// ****************************************************************************
/// \brief	Method to initialise all labels on sub-grids.
///
///			Boundary labels are set by considering parent labels on overlapping
///			sites and then assigning child labels appropriately.
///
/// \param	pGrid	reference to parent grid.
void GridObj::LBM_initRefinedLab (GridObj& pGrid) {

	// Get parent local grid sizes
	size_t Np_lim = pGrid.N_lim;
	size_t Mp_lim = pGrid.M_lim;
	size_t Kp_lim = pGrid.K_lim;
	
	// Declare indices global and local
	int i, j, k, local_i, local_j, local_k;

	// Loop over global indices of refined patch and add "TL to lower" labels
    for (i = (int)RefXstart[pGrid.level][region_number]; i <= (int)RefXend[pGrid.level][region_number]; i++) {
        for (j = (int)RefYstart[pGrid.level][region_number]; j <= (int)RefYend[pGrid.level][region_number]; j++) {
#if (L_dims == 3)
			for (k = (int)RefZstart[pGrid.level][region_number]; k <= (int)RefZend[pGrid.level][region_number]; k++)
#else
			k = 0;
#endif
			{

			// Only act if the site is on parent rank (inc overlap) to avoid out of bounds errors
			if ( GridUtils::isOnThisRank(i,j,k,pGrid) ) {

#ifdef L_BUILD_FOR_MPI						
				// Compute local indices to access LatTyp array on parent
				std::vector<int> locals;
				GridUtils::global_to_local(i,j,k,&pGrid,locals);
				local_i = locals[0];
				local_j = locals[1];
				local_k = locals[2];
#else
				local_i = i;
				local_j = j;
				local_k = k;						
#endif

				// If on the edge of the global refined patch and it is simply a fluid then it is TL so label
				if	(
					(i == RefXstart[pGrid.level][region_number] || i == RefXend[pGrid.level][region_number]) ||
					(j == RefYstart[pGrid.level][region_number] || j == RefYend[pGrid.level][region_number])
#if (L_dims == 3)
					|| (k == RefZstart[pGrid.level][region_number] || k == RefZend[pGrid.level][region_number])
#endif
					) {

					// If parent site is fluid site then correct label
					if (pGrid.LatTyp(local_i,local_j,local_k,Mp_lim,Kp_lim) == eFluid) {
						// Change to "TL to lower" label
						pGrid.LatTyp(local_i,local_j,local_k,Mp_lim,Kp_lim) = eTransitionToFiner;
					}

				// Else it is not on the edges but in the middle of the global refined patch
				} else if (pGrid.LatTyp(local_i,local_j,local_k,Mp_lim,Kp_lim) == eFluid) {
					// Label it a "refined" site
					pGrid.LatTyp(local_i,local_j,local_k,Mp_lim,Kp_lim) = eRefined;

				}

			}

			}            
        }
    }
    

	// Generate grid type matrices for this level //

	// Declare array for parent site indices
	std::vector<int> p;
	eType par_label;
    
    // Loop over sub-grid local indices and add labels based on parent site labels
    for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
#if (L_dims == 3)
			for (int k = 0; k < K_lim; k++)
#else
			int k = 0;
#endif
			{
				// Get parent site indices as locals for the current sub-grid local site
				p = GridUtils::getCoarseIndices(
					i, CoarseLimsX[0], 
					j, CoarseLimsY[0], 
					k, CoarseLimsZ[0]
				);
				
				// Get parent site label using local indices
				par_label = pGrid.LatTyp(p[0],p[1],p[2],Mp_lim,Kp_lim);
								
				// If parent is a "TL to lower" then add "TL to upper" label
				if (par_label == eTransitionToFiner) { 
					LatTyp(i,j,k,M_lim,K_lim) = eTransitionToCoarser;
                
					// Else if parent is a "refined" label then label as coarse
				} else if (par_label == eRefined) { 
					LatTyp(i,j,k,M_lim,K_lim) = eFluid;
                
					// Else parent label is some other kind of boundary so copy the
					// label to retain the behaviour onto this grid
				} else { 
					LatTyp(i,j,k,M_lim,K_lim) = par_label;
					
					// If last site to be updated in fine block, change parent label
					// to ensure boundary values are pulled from fine grid
					if ((j % 2) != 0 && (i % 2) != 0) {
						if (pGrid.LatTyp(p[0],p[1],p[2],Mp_lim,Kp_lim) == eSolid || pGrid.LatTyp(p[0],p[1],p[2],Mp_lim,Kp_lim) == eRefinedSolid) {
							pGrid.LatTyp(p[0],p[1],p[2],Mp_lim,Kp_lim) = eRefinedSolid;

						} else if (pGrid.LatTyp(p[0],p[1],p[2],Mp_lim,Kp_lim) == eSymmetry || pGrid.LatTyp(p[0],p[1],p[2],Mp_lim,Kp_lim) == eRefinedSymmetry) {
							pGrid.LatTyp(p[0],p[1],p[2],Mp_lim,Kp_lim) = eRefinedSymmetry;
						
						} else if (pGrid.LatTyp(p[0],p[1],p[2],Mp_lim,Kp_lim) == eInlet || pGrid.LatTyp(p[0],p[1],p[2],Mp_lim,Kp_lim) == eRefinedInlet) {
							pGrid.LatTyp(p[0],p[1],p[2],Mp_lim,Kp_lim) = eRefinedInlet;

						}
					}
				}
            }
        }
    }

	
	// Try to add the solid block labels
	LBM_initSolidLab();
}

// ***************************************************************************************************
