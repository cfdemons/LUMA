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

using namespace std;

// ****************************************************************************
/// \brief	Method to import an input profile from a file.
///
///			Input data may be over- or under-sampled but it must span the 
///			physical dimensions of the inlet otherwise the software does
///			not known how to scale the data to fit. Inlet profile is always
///			assumed to be oriented vertically (y-direction).
void GridObj::LBM_init_getInletProfile() {

	size_t j;
	std::vector<double> ybuffer, uxbuffer, uybuffer, uzbuffer;

#ifdef L_PARABOLIC_INLET

	// Resize vectors
	ux_in.resize(M_lim);
	uy_in.resize(M_lim);
	uz_in.resize(M_lim);

	// Loop over site positions (for left hand inlet, y positions)
	for (j = 0; j < M_lim; j++) {

		// Set the inlet velocity profile values
		ux_in[j] = GridUnits::ud2ulbm(L_UMAX,this) * (1 - pow((YPos[j] - (L_BY - 2 * dh) / 2) / ((L_BY - 2 * dh) / 2), 2));
		uy_in[j] = 0.0;
		uz_in[j] = 0.0;
	}

#else

	size_t i;
	double y, tmp;

	// Indicate to log
	*GridUtils::logfile << "Loading inlet profile..." << std::endl;

	// Buffer information from file
	std::ifstream inletfile;
	inletfile.open("./input/inlet_profile.in", std::ios::in);
	if (!inletfile.is_open()) {
		// Error opening file
		L_ERROR("Cannot open inlet profile file named \"inlet_profile.in\". Exiting.", GridUtils::logfile);

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
		L_ERROR("Failed to read in inlet profile data. Exiting.", GridUtils::logfile);
	}
#endif // L_PARABOLIC_INLET


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
				for (size_t d = 0; d < L_DIMS; d++) {

					// No flow case
					u(i,j,k,d,M_lim,K_lim,L_DIMS) = 0.0;
				}

#else
				/* Input velocity is specified by individual vectors for x, y and z which 
				 * have either been read in from an input file or defined by an expression 
				 * given in the definitions.
				 */
				u(i,j,k,0,M_lim,K_lim,L_DIMS) = GridUnits::ud2ulbm(L_UX0,this);
				u(i,j,k,1,M_lim,K_lim,L_DIMS) = GridUnits::ud2ulbm(L_UY0,this);
#if (L_DIMS == 3)
				u(i,j,k,2,M_lim,K_lim,L_DIMS) = GridUnits::ud2ulbm(L_UZ0,this);
#endif


#endif
				if (LatTyp(i, j, k, M_lim, K_lim) == eSolid) {
					u(i, j, k, 0, M_lim, K_lim, L_DIMS) = 0.0;
					u(i, j, k, 1, M_lim, K_lim, L_DIMS) = 0.0;
#if (L_DIMS == 3)
					u(i, j, k, 2, M_lim, K_lim, L_DIMS) = 0.0;
#endif

				}


			}
		}
	}

}


// ****************************************************************************
/// \brief	Method to initialise the lattice density.
void GridObj::LBM_initRho ( ) {

	// Loop over grid
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {		
			for (int k = 0; k < K_lim; k++) {

				// Uniform Density
				rho(i,j,k,M_lim,K_lim) = L_RHOIN;
			}
		}
	}

}

// ****************************************************************************
/// \brief	Method to initialise all L0 lattice quantities.
void GridObj::LBM_initGrid() {

#ifdef L_INIT_VERBOSE
	*GridUtils::logfile << "Initialising grid level 0..." << std::endl;
#endif

	// Get GM instance
	GridManager *gm = GridManager::getInstance();

#ifdef L_BUILD_FOR_MPI

	// Get MpiManager instance
	MpiManager *mpim = MpiManager::getInstance();

#endif
	
	// Store physical spacing
	// Global dimensions
	double Lx = L_BX;
	double Ly = L_BY;
	double Lz = L_BZ;
	dh = Lx / static_cast<double>(L_N);

	// Store temporal spacing
	dt = L_TIMESTEP;

	//Gravity in LBM units
	g = GridUnits::fd2flbm(L_GRAVITY_FORCE, this);

	//Reference velocity in LBM units
	uref = GridUnits::ud2ulbm(1, this);

	


	////////////////////////////
	// Check input parameters //
	////////////////////////////

#if (L_DIMS == 3)
	// Check that lattice volumes are cubes in 3D
	if (abs((Lx / L_N) - (Ly / L_M)) > L_SMALL_NUMBER || abs((Lx / L_N) - (Lz / L_K)) > L_SMALL_NUMBER) {
		L_ERROR("Need to have lattice volumes which are cubes -- either change L_N/L_M/L_K or change domain dimensions. Exiting.", 
			GridUtils::logfile);
	}
	
#else
	// 2D so need square lattice cells
	if ( abs((Lx/L_N) - (Ly/L_M)) > L_SMALL_NUMBER ) {
		L_ERROR("Need to have lattice cells which are squares -- either change L_N/L_M or change domain dimensions. Exiting.",
			GridUtils::logfile);
	}

#endif


    ///////////////////
	// Checks passed //
	///////////////////

	// Get local grid sizes (includes halo)
	N_lim = gm->local_size[eXDirection];
	M_lim = gm->local_size[eYDirection];
#if (L_DIMS == 3)
	K_lim = gm->local_size[eZDirection];
#else
	K_lim = 1;
#endif
	

	// L0 lattice site POSITION VECTORS
	/* When using MPI:
	 * The overlap is assume periodic in all directions on L0.
	 */

#ifdef L_BUILD_FOR_MPI

	// Create position vectors excluding overlap
	XPos = GridUtils::linspace(mpim->rank_core_edge[eXMin][mpim->my_rank] + dh / 2, mpim->rank_core_edge[eXMax][mpim->my_rank] - dh / 2, N_lim - 2);
	YPos = GridUtils::linspace(mpim->rank_core_edge[eYMin][mpim->my_rank] + dh / 2, mpim->rank_core_edge[eYMax][mpim->my_rank] - dh / 2, M_lim - 2);
	ZPos = GridUtils::linspace(mpim->rank_core_edge[eZMin][mpim->my_rank] + dh / 2, mpim->rank_core_edge[eZMax][mpim->my_rank] - dh / 2, K_lim - 2);

	// Add overlap sites taking into account periodicity
	XPos.insert( XPos.begin(), fmod(XPos[0] - dh + Lx, Lx) );
	XPos.insert( XPos.end(), fmod(XPos[XPos.size() - 1] + dh + Lx, Lx) );

	YPos.insert( YPos.begin(), fmod(YPos[0] - dh + Ly, Ly) );
	YPos.insert( YPos.end(), fmod(YPos[YPos.size() - 1] + dh + Ly, Ly) );
#if (L_DIMS == 3)
	ZPos.insert( ZPos.begin(), fmod(ZPos[0] - dh + Lz, Lz) );
	ZPos.insert( ZPos.end(), fmod(ZPos[ZPos.size() - 1] + dh + Lz, Lz) );
#endif

	// Update the sender/recv layer positions in the MpiManager

	// X
	mpim->sender_layer_pos.X[eLeftMin]	= XPos[1] - dh/2;						mpim->sender_layer_pos.X[eLeftMax]	= XPos[1] + dh/2;
	mpim->sender_layer_pos.X[eRightMin] = XPos[gm->local_size[0] - 2] - dh / 2;	mpim->sender_layer_pos.X[eRightMax] = XPos[gm->local_size[0] - 2] + dh / 2;
	mpim->recv_layer_pos.X[eLeftMin]	= XPos[0] - dh/2;						mpim->recv_layer_pos.X[eLeftMax]	= XPos[0] + dh/2;
	mpim->recv_layer_pos.X[eRightMin]	= XPos[gm->local_size[0] - 1] - dh / 2;	mpim->recv_layer_pos.X[eRightMax]	= XPos[gm->local_size[0] - 1] + dh / 2;
	// Y
	mpim->sender_layer_pos.Y[eLeftMin]	= YPos[1] - dh/2;						mpim->sender_layer_pos.Y[eLeftMax]	= YPos[1] + dh/2;
	mpim->sender_layer_pos.Y[eRightMin] = YPos[gm->local_size[1] - 2] - dh / 2;	mpim->sender_layer_pos.Y[eRightMax] = YPos[gm->local_size[1] - 2] + dh / 2;
	mpim->recv_layer_pos.Y[eLeftMin]	= YPos[0] - dh/2;						mpim->recv_layer_pos.Y[eLeftMax]	= YPos[0] + dh/2;
	mpim->recv_layer_pos.Y[eRightMin]	= YPos[gm->local_size[1] - 1] - dh / 2;	mpim->recv_layer_pos.Y[eRightMax]	= YPos[gm->local_size[1] - 1] + dh / 2;
	// Z
#if (L_DIMS == 3)
	mpim->sender_layer_pos.Z[eLeftMin]	= ZPos[1] - dh/2;						mpim->sender_layer_pos.Z[eLeftMax]	= ZPos[1] + dh/2;
	mpim->sender_layer_pos.Z[eRightMin] = ZPos[gm->local_size[2] - 2] - dh/2;	mpim->sender_layer_pos.Z[eRightMax] = ZPos[gm->local_size[2] - 2] + dh/2;
	mpim->recv_layer_pos.Z[eLeftMin]	= ZPos[0] - dh/2;						mpim->recv_layer_pos.Z[eLeftMax]	= ZPos[0] + dh/2;
	mpim->recv_layer_pos.Z[eRightMin]	= ZPos[gm->local_size[2] - 1] - dh/2;	mpim->recv_layer_pos.Z[eRightMax]	= ZPos[gm->local_size[2] - 1] + dh/2;
#endif

#ifdef L_MPI_VERBOSE
	*mpim->logout << "X sender layers are: " << mpim->sender_layer_pos.X[eLeftMin] << " -- " << mpim->sender_layer_pos.X[eLeftMax] << " (min edge) , " << mpim->sender_layer_pos.X[eRightMin] << " -- " << mpim->sender_layer_pos.X[eRightMax] << " (max edge)" << std::endl;
	*mpim->logout << "X recv layers are: " << mpim->recv_layer_pos.X[eLeftMin] << " -- " << mpim->recv_layer_pos.X[eLeftMax] << " (min edge) , " << mpim->recv_layer_pos.X[eRightMin] << " -- " << mpim->recv_layer_pos.X[eRightMax] << " (max edge)" << std::endl;

	*mpim->logout << "Y sender layers are: " << mpim->sender_layer_pos.Y[eLeftMin] << " -- " << mpim->sender_layer_pos.Y[eLeftMax] << " (min edge) , " << mpim->sender_layer_pos.Y[eRightMin] << " -- " << mpim->sender_layer_pos.Y[eRightMax] << " (max edge)" << std::endl;
	*mpim->logout << "Y recv layers are: " << mpim->recv_layer_pos.Y[eLeftMin] << " -- " << mpim->recv_layer_pos.Y[eLeftMax] << " (min edge) , " << mpim->recv_layer_pos.Y[eRightMin] << " -- " << mpim->recv_layer_pos.Y[eRightMax] << " (max edge)" << std::endl;

#if (L_DIMS == 3)
	*mpim->logout << "Z sender layers are: " << mpim->sender_layer_pos.Z[eLeftMin] << " -- " << mpim->sender_layer_pos.Z[eLeftMax] << " (min edge) , " << mpim->sender_layer_pos.Z[eRightMin] << " -- " << mpim->sender_layer_pos.Z[eRightMax] << " (max edge)" << std::endl;
	*mpim->logout << "Z recv layers are: " << mpim->recv_layer_pos.Z[eLeftMin] << " -- " << mpim->recv_layer_pos.Z[eLeftMax] << " (min edge) , " << mpim->recv_layer_pos.Z[eRightMin] << " -- " << mpim->recv_layer_pos.Z[eRightMax] << " (max edge)" << std::endl;
#endif
#endif

#else
	// When not builiding for MPI positions are straightforward
	XPos = GridUtils::linspace(dh / 2, L_BX - dh / 2, static_cast<int>(L_N));
	YPos = GridUtils::linspace(dh / 2, L_BY - dh / 2, static_cast<int>(L_M));
	ZPos = GridUtils::linspace(dh / 2, L_BZ - dh / 2, static_cast<int>(L_K));

#endif	// L_BUILD_FOR_MPI

	// Absolute origins (position of first site on grid)
	XOrigin = dh / 2.0;
	YOrigin = dh / 2.0;
#if (L_DIMS == 3)
	ZOrigin = dh / 2.0;
#else
	ZOrigin = 0.0;
#endif

	

	// Define TYPING MATRICES
	LatTyp.resize(N_lim * M_lim * K_lim);

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
	u.resize( N_lim*M_lim*K_lim*L_DIMS );
	LBM_initVelocity();
	
	// Density field
	rho.resize( N_lim*M_lim*K_lim );
	LBM_initRho();

	// Cartesian force vector
	force_xyz.resize(N_lim*M_lim*K_lim*L_DIMS, 0.0);

	// Lattice force vector
	force_i.resize(N_lim*M_lim*K_lim*L_NUM_VELS, 0.0);

	// Time averaged quantities
	rho_timeav.resize(N_lim*M_lim*K_lim, 0.0);
	ui_timeav.resize(N_lim*M_lim*K_lim*L_DIMS, 0.0);
	uiuj_timeav.resize(N_lim*M_lim*K_lim*(3*L_DIMS-3), 0.0);


	// Initialise L0 POPULATION matrices (f, feq)
	f.resize( N_lim*M_lim*K_lim*L_NUM_VELS );
	feq.resize( N_lim*M_lim*K_lim*L_NUM_VELS );
	fNew.resize(N_lim * M_lim * K_lim * L_NUM_VELS);


	// Loop over grid
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {
				for (int v = 0; v < L_NUM_VELS; v++) {

					// Initialise f to feq
					f(i, j, k, v, M_lim, K_lim, L_NUM_VELS) = 
						_LBM_equilibrium_opt(k + j * K_lim + i * M_lim * K_lim, v);

				}
			}
		}
	}
	feq = f; // Make feq = feq too
	fNew = f;


	// Initialise OTHER parameters
//	// Compute kinematic viscosity based on target Reynolds number
//#if defined L_IBM_ON && defined L_INSERT_CIRCLE_SPHERE
//	// If IBM circle use diameter (in lattice units i.e. rescale wrt to physical spacing)
//	nu = (L_IBB_R*2 / dh) * L_UREF / L_RE;
//#elif defined L_IBM_ON && defined L_INSERT_RECTANGLE_CUBOID
//	// If IBM rectangle use y-dimension (in lattice units)
//	nu = (L_IBB_L / dh) * L_UREF / L_RE;
//#elif defined L_IBM_ON && defined L_IBB_FROM_FILE
//	// If IBM object read from file then use scale length as reference
//	nu = (L_IBB_REF_LENGTH / dh) * L_UREF / L_RE;
//#elif defined L_SOLID_FROM_FILE
//	// Use object length
//	nu = (L_OBJECT_REF_LENGTH / dh) * L_UREF / L_RE;
//#elif defined L_BFL_ON
//	// Use bfl body length
//	nu = (L_BFL_REF_LENGTH / dh) * L_UREF / L_RE;
//#elif defined WALLS_ON
//	// If no object then use domain height
//	nu = (L_BY / dh - std::round(L_WALL_THICKNESS_BOTTOM / dh) - std::round(L_WALL_THICKNESS_TOP / dh)) * L_UREF / L_RE;	// Based on actual width of channel
//#else
//	// Use reference length of 1.0
//	nu = (1.0 / dh) * L_UREF / L_RE;
//#endif

#ifdef L_NU
	nu = GridUnits::nud2nulbm(L_NU, this);
#else
	nu = GridUnits::nud2nulbm(1/static_cast<double>(L_RE), this);
#endif

	// Relaxation frequency on L0
	// Assign relaxation frequency using lattice viscosity
	omega = 1 / ( (nu / pow(cs,2)) + .5 );

	/* Above is valid for L0 only when dh = 1 -- general expression is:
	 * omega = 1 / ( ( (nu * dt) / (pow(cs,2)*pow(dh,2)) ) + .5 );
	 */

#ifndef L_BUILD_FOR_MPI
	// Update the writable data in the grid manager
	// When using MPI this is done when building the communicators.
	gm->createWritableDataStore(this);
#endif

#ifdef L_INIT_VERBOSE
	*GridUtils::logfile << "Initialisation Complete.";
#endif

}
	

// ****************************************************************************
/// \brief	Method to initialise all sub-grid quantities.
/// \param	pGrid	reference to parent grid.
void GridObj::LBM_initSubGrid (GridObj& pGrid) {

#ifdef L_INIT_VERBOSE
	*GridUtils::logfile << "Initialising sub-grid level " << level << ", region " << region_number << "..." << std::endl;
#endif

	// Get GM instance
	GridManager *gm = GridManager::getInstance();
	int gm_idx = level + region_number * L_NUM_LEVELS;

	// Define scales
	dh = pGrid.dh / 2.0;
	dt = pGrid.dt / 2.0;
	g = pGrid.g;
	uref = pGrid.uref;
	
	/* Get coarse grid refinement limits as indicies local to the parent grid
	 * on this rank. */
	LBM_initGridToGridMappings(pGrid);
	
	/* Get local grid size of the sub grid based on local ijk limits.
	 * i.e. how much of the parent grid it covers. Volumetric formulation with refinement
	 * factor of 2 makes this easy. */
	int local_size[L_DIMS] = {
		static_cast<int>((CoarseLimsX[eMaximum] - CoarseLimsX[eMinimum] + .5) * 2) + 1,
		static_cast<int>((CoarseLimsY[eMaximum] - CoarseLimsY[eMinimum] + .5) * 2) + 1
#if (L_DIMS == 3)
		,
		static_cast<int>((CoarseLimsZ[eMaximum] - CoarseLimsZ[eMinimum] + .5) * 2) + 1
#endif
	};

	// Set grid sizes
	N_lim = local_size[eXDirection];
	M_lim = local_size[eYDirection];
#if (L_DIMS == 3)
	K_lim = local_size[eZDirection];
#else
	K_lim = 1;
#endif
	
	
	// Generate POSITION VECTORS of nodes
	
	// Populate the position vectors
	LBM_initPositionVector(pGrid.XPos[CoarseLimsX[eMinimum]] - dh / 2.0, pGrid.XPos[CoarseLimsX[eMaximum]] - dh / 2.0, eXDirection);
	LBM_initPositionVector(pGrid.YPos[CoarseLimsY[eMinimum]] - dh / 2.0, pGrid.YPos[CoarseLimsY[eMaximum]] - dh / 2.0, eYDirection);
#if L_DIMS == 3
	LBM_initPositionVector(pGrid.ZPos[CoarseLimsZ[eMinimum]] - dh / 2.0, pGrid.ZPos[CoarseLimsZ[eMaximum]] - dh / 2.0, eZDirection);
#else
	ZPos.insert( ZPos.begin(), 1 ); // 2D default
#endif
		

	// Global edge origins (voxel centre position of first cell on grid)
	XOrigin = gm->global_edges[eXMin][gm_idx] + dh / 2;
	YOrigin = gm->global_edges[eYMin][gm_idx] + dh / 2;
#if (L_DIMS == 3)
	ZOrigin = gm->global_edges[eZMin][gm_idx] + dh / 2;
#else
	ZOrigin = 0.0;
#endif

	
	// Generate TYPING MATRICES

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
	u.resize(N_lim * M_lim * K_lim * L_DIMS);
	LBM_initVelocity( );

	// Density
	rho.resize(N_lim * M_lim * K_lim);
	LBM_initRho( );


	// Cartesian force vector
	force_xyz.resize(N_lim * M_lim * K_lim * L_DIMS, 0.0);

	// Lattice force vector
	force_i.resize(N_lim * M_lim * K_lim * L_NUM_VELS, 0.0);

	// Time averaged quantities
	rho_timeav.resize(N_lim*M_lim*K_lim, 0.0);
	ui_timeav.resize(N_lim*M_lim*K_lim*L_DIMS, 0.0);
	uiuj_timeav.resize(N_lim*M_lim*K_lim*(3*L_DIMS-3), 0.0);


	// Generate POPULATION MATRICES for lower levels
	// Resize
	f.resize(N_lim * M_lim * K_lim * L_NUM_VELS);
	feq.resize(N_lim * M_lim * K_lim * L_NUM_VELS);
	fNew.resize(N_lim * M_lim * K_lim * L_NUM_VELS);


	// Loop over grid
	for (int i = 0; i != N_lim; ++i) {
		for (int j = 0; j != M_lim; ++j) {
			for (int k = 0; k != K_lim; ++k) {
				for (int v = 0; v < L_NUM_VELS; v++) {
					
					// Initialise f to feq
					f(i, j, k, v, M_lim, K_lim, L_NUM_VELS) = 
						_LBM_equilibrium_opt(k + j * K_lim + i * M_lim * K_lim, v);

				}
			}
		}
	}
	feq = f; // Set feq to feq
	fNew = f;

	// Compute relaxation time from coarser level assume refinement by factor of 2
	omega = 1.0 / ( ( (1.0 / pGrid.omega - 0.5) * 2.0) + 0.5);

	// Lattice viscosity is constant across subgrids
	nu = pGrid.nu;

#ifdef L_INIT_VERBOSE
	*GridUtils::logfile << "Initialisation Complete." << std::endl;
#endif

}


// ****************************************************************************
/// \brief	Method to initialise the mapping parameters between this grid and 
///			the supplied parent.
///
///			The mappings computed by this method are local ijk references to allow
///			correct coupling during multi-grid operations.
///
/// \param	pGrid	reference to parent grid.
void GridObj::LBM_initGridToGridMappings(GridObj& pGrid)
{

	// If not using MPI then it is easy
	/* Coarse Limits stored as local indices as they are used to map between the
	 * fine and coarse grid cells during multi - grid operations.Therefore, we
	 * must only store local values relevant to the grid on the rank to ensure
	 * mapping is correct.
	 *
	 * When not using MPI, these can be computed from the definitions file values
	 * which specify the extent of the refinement at that level.
	 * However, when using MPI, the edges of the refined grid might not be on this
	 * rank so we must round the coarse limits to the edge of the parent grid so
	 * the correct offset is supplied to the mapping routine.
	 *
	 * When dealing with sub - grids embedded in walls, it is more complicated.If the
	 * sub - grid starts on a max receiver layer which is also a periodic overlap then
	 * we need to make sure the limits are set properly.Likewise if the sub - grid ends
	 * on a min receiver layer which is also periodic. */

	int position = 0;
	int position2 = 0;
	eLocationOnRank loc = eNone;
	double subgrid_start_voxel_centre;
	double subgrid_end_voxel_centre;
	int gm_idx = level + region_number * L_NUM_LEVELS;
	GridManager *gm = GridManager::getInstance();

#ifdef L_BUILD_FOR_MPI
	MpiManager *mpim = MpiManager::getInstance();
#endif
	

#ifdef L_INIT_VERBOSE
	*GridUtils::logfile << "Cell width on parent = " << pGrid.dh << std::endl;
#endif

	// X //

	// Set voxel centre positions from knowledge of grid edges
	subgrid_start_voxel_centre = gm->global_edges[eXMin][gm_idx] + dh / 2.0;
	subgrid_end_voxel_centre = gm->global_edges[eXMax][gm_idx] - dh / 2.0;

#ifdef L_INIT_VERBOSE
	*GridUtils::logfile << "X: Start and end voxel centres are at: " << 
		subgrid_start_voxel_centre << " & " << subgrid_end_voxel_centre << std::endl;
#endif

	// Set TL to on by default
	gm->subgrid_tlayer_key[eXMin][gm_idx - 1] = true;
	gm->subgrid_tlayer_key[eXMax][gm_idx - 1] = true;

	// Start Limit: Find whether edge of refined region is on this grid at all and return its local index

	// Starts on this rank
	if (GridUtils::isOnThisRank(subgrid_start_voxel_centre, eXDirection, &loc, &pGrid, &position))
	{
		// Set limit to ijk of enclosing voxel
		CoarseLimsX[eMinimum] = position;

#ifdef L_INIT_VERBOSE
		*GridUtils::logfile << "X: Starts on this rank @ local " << position << " Loc = " << loc << std::endl;
#endif
	}

#ifdef L_BUILD_FOR_MPI

	/* Starts on some rank to the left of this one. Note: Using the core edge position
	* is fine as if it was in the halo the above call would have returned true. */
	else if (subgrid_start_voxel_centre < mpim->rank_core_edge[eXMin][mpim->my_rank])
	{
		// Set limit to start edge of the rank
		CoarseLimsX[eMinimum] = 0;

#ifdef L_INIT_VERBOSE
		*GridUtils::logfile << "X: Starts on earlier rank" << std::endl;
#endif
	}

#endif // L_BUILD_FOR_MPI


	// End Limit: Find whether last sub-grid site position is on this grid at all and return its local index

	// Ends on this rank
	if (GridUtils::isOnThisRank(subgrid_end_voxel_centre, eXDirection, &loc, &pGrid, &position))
	{
		// Set limit to ijk of enclosing voxel
		CoarseLimsX[eMaximum] = position;

#ifdef L_INIT_VERBOSE
		*GridUtils::logfile << "X: Ends on this rank @ local " << position << " Loc = " << loc << std::endl;
#endif
	}

#ifdef L_BUILD_FOR_MPI

	// Ends on some rank to the right of this one
	else if (subgrid_end_voxel_centre > mpim->rank_core_edge[eXMax][mpim->my_rank])
	{
		// Set limit to end edge of the rank
		CoarseLimsX[eMaximum] = pGrid.N_lim - 1;

#ifdef L_INIT_VERBOSE
		*GridUtils::logfile << "X: Ends on later rank" << std::endl;
#endif
	}

	// Else the end must be on a rank to the left and hence grid wraps periodically so set end to right-hand edge
	else if (subgrid_end_voxel_centre < mpim->rank_core_edge[eXMin][mpim->my_rank])
	{
		// Set grid limits to end edge of the rank
		CoarseLimsX[eMaximum] = pGrid.N_lim - 1;

#ifdef L_INIT_VERBOSE
		*GridUtils::logfile << "X: Ends on earlier rank" << std::endl;
#endif
	}

#endif

	/* If start and finish on same process then we test to see if sites are next
	* to each other and coalesce the TL edges. */

	if (GridUtils::isOnThisRank(subgrid_start_voxel_centre, eXDirection, &loc, &pGrid, &position) &&
		GridUtils::isOnThisRank(subgrid_end_voxel_centre, eXDirection, &loc, &pGrid, &position2))
	{

		// Check sites are next to each other in space
		if (abs(pGrid.XPos[position] - pGrid.XPos[position2]) - pGrid.dh < L_SMALL_NUMBER)
		{
			// Set as if middle of complete block
			CoarseLimsX[eMinimum] = 0;
			CoarseLimsX[eMaximum] = pGrid.N_lim - 1;

			// Tell mpim that TL doesn't exist so HDF5 writer does not try to exclude valid sites
			gm->subgrid_tlayer_key[eXMin][gm_idx - 1] = false;
			gm->subgrid_tlayer_key[eXMax][gm_idx - 1] = false;

#ifdef L_INIT_VERBOSE
			*GridUtils::logfile << "X: Periodically wrapped and joined" << std::endl;
#endif
		}
	}



	// Y //
	subgrid_start_voxel_centre = gm->global_edges[eYMin][gm_idx] + dh / 2;
	subgrid_end_voxel_centre = gm->global_edges[eYMax][gm_idx] - dh / 2;
	gm->subgrid_tlayer_key[eYMin][gm_idx - 1] = true;
	gm->subgrid_tlayer_key[eYMax][gm_idx - 1] = true;

#ifdef L_INIT_VERBOSE
	*GridUtils::logfile << "Y: Start and end voxel centres are at: " <<
		subgrid_start_voxel_centre << " & " << subgrid_end_voxel_centre << std::endl;
#endif

	if (GridUtils::isOnThisRank(subgrid_start_voxel_centre, eYDirection, &loc, &pGrid, &position))
	{
		CoarseLimsY[eMinimum] = position;

#ifdef L_INIT_VERBOSE
		*GridUtils::logfile << "Y: Starts on this rank @ local " << position << " Loc = " << loc << std::endl;
#endif
	}

#ifdef L_BUILD_FOR_MPI
	else if (subgrid_start_voxel_centre < mpim->rank_core_edge[eYMin][mpim->my_rank])
	{
		CoarseLimsY[eMinimum] = 0;

#ifdef L_INIT_VERBOSE
		*GridUtils::logfile << "Y: Starts on earlier rank" << std::endl;
#endif
	}
#endif

	if (GridUtils::isOnThisRank(subgrid_end_voxel_centre, eYDirection, &loc, &pGrid, &position))
	{
		CoarseLimsY[eMaximum] = position;

#ifdef L_INIT_VERBOSE
		*GridUtils::logfile << "Y: Ends on this rank @ local " << position << " Loc = " << loc << std::endl;
#endif
	}

#ifdef L_BUILD_FOR_MPI
	else if (subgrid_end_voxel_centre > mpim->rank_core_edge[eYMax][mpim->my_rank])
	{
		CoarseLimsY[eMaximum] = pGrid.M_lim - 1;

#ifdef L_INIT_VERBOSE
		*GridUtils::logfile << "Y: Ends on later rank" << std::endl;
#endif
	}
	else if (subgrid_end_voxel_centre < mpim->rank_core_edge[eYMin][mpim->my_rank])
	{
		CoarseLimsY[eMaximum] = pGrid.M_lim - 1;

#ifdef L_INIT_VERBOSE
		*GridUtils::logfile << "Y: Ends on earlier rank" << std::endl;
#endif
	}
#endif

	if (GridUtils::isOnThisRank(subgrid_start_voxel_centre, eYDirection, &loc, &pGrid, &position) &&
		GridUtils::isOnThisRank(subgrid_end_voxel_centre, eYDirection, &loc, &pGrid, &position2))
	{
		if (abs(pGrid.YPos[position] - pGrid.YPos[position2]) - pGrid.dh < L_SMALL_NUMBER)
		{
			CoarseLimsY[eMinimum] = 0;
			CoarseLimsY[eMaximum] = pGrid.M_lim - 1;
			gm->subgrid_tlayer_key[eYMin][gm_idx - 1] = false;
			gm->subgrid_tlayer_key[eYMax][gm_idx - 1] = false;

#ifdef L_INIT_VERBOSE
			*GridUtils::logfile << "Y: Periodically wrapped and joined" << std::endl;
#endif
		}
	}


#if (L_DIMS == 3)
	// Z //
	subgrid_start_voxel_centre = gm->global_edges[eZMin][gm_idx] + dh / 2;
	subgrid_end_voxel_centre = gm->global_edges[eZMax][gm_idx] - dh / 2;
	gm->subgrid_tlayer_key[eZMin][gm_idx - 1] = true;
	gm->subgrid_tlayer_key[eZMax][gm_idx - 1] = true;

#ifdef L_INIT_VERBOSE
	*GridUtils::logfile << "Z: Start and end voxel centres are at: " <<
		subgrid_start_voxel_centre << " & " << subgrid_end_voxel_centre << std::endl;
#endif

	if (GridUtils::isOnThisRank(subgrid_start_voxel_centre, eZDirection, &loc, &pGrid, &position))
	{
		CoarseLimsZ[eMinimum] = position;

#ifdef L_INIT_VERBOSE
		*GridUtils::logfile << "Z: Starts on this rank @ local " << position << " Loc = " << loc << std::endl;
#endif
	}

#ifdef L_BUILD_FOR_MPI
	else if (subgrid_start_voxel_centre < mpim->rank_core_edge[eZMin][mpim->my_rank])
	{
		CoarseLimsZ[eMinimum] = 0;

#ifdef L_INIT_VERBOSE
		*GridUtils::logfile << "Z: Starts on earlier rank" << std::endl;
#endif
	}
#endif

	if (GridUtils::isOnThisRank(subgrid_end_voxel_centre, eZDirection, &loc, &pGrid, &position))
	{
		CoarseLimsZ[eMaximum] = position;

#ifdef L_INIT_VERBOSE
		*GridUtils::logfile << "Z: Ends on this rank @ local " << position << " Loc = " << loc << std::endl;
#endif
	}

#ifdef L_BUILD_FOR_MPI
	else if (subgrid_end_voxel_centre > mpim->rank_core_edge[eZMax][mpim->my_rank])
	{
		CoarseLimsZ[eMaximum] = pGrid.K_lim - 1;

#ifdef L_INIT_VERBOSE
		*GridUtils::logfile << "Z: Ends on later rank" << std::endl;
#endif
	}
	else if (subgrid_end_voxel_centre < mpim->rank_core_edge[eZMin][mpim->my_rank])
	{
		CoarseLimsZ[eMaximum] = pGrid.K_lim - 1;

#ifdef L_INIT_VERBOSE
		*GridUtils::logfile << "Z: Ends on earlier rank" << std::endl;
#endif
	}
#endif

	if (GridUtils::isOnThisRank(subgrid_start_voxel_centre, eZDirection, &loc, &pGrid, &position) &&
		GridUtils::isOnThisRank(subgrid_end_voxel_centre, eZDirection, &loc, &pGrid, &position2))
	{
		if (abs(pGrid.ZPos[position] - pGrid.ZPos[position2]) - pGrid.dh < L_SMALL_NUMBER)
		{
			CoarseLimsZ[eMinimum] = 0;
			CoarseLimsZ[eMaximum] = pGrid.K_lim - 1;
			gm->subgrid_tlayer_key[eZMin][gm_idx - 1] = false;
			gm->subgrid_tlayer_key[eZMax][gm_idx - 1] = false;

#ifdef L_INIT_VERBOSE
			*GridUtils::logfile << "Z: Periodically wrapped and joined" << std::endl;
#endif
		}
	}
#else
	// Reset the refined region z-limits if only 2D
	for (int i = 0; i < 2; i++) {
		CoarseLimsZ[i] = 0;
	}
#endif

	/* If sub-grid wraps periodically we do not support it wrapping back round to the same rank
	* again as this will confuse the mapping function so we error here. */

	if (
		(CoarseLimsX[eMaximum] < CoarseLimsX[eMinimum] && CoarseLimsX[eMinimum] - 1 != CoarseLimsX[eMaximum]) ||
		(CoarseLimsY[eMaximum] < CoarseLimsY[eMinimum] && CoarseLimsY[eMinimum] - 1 != CoarseLimsY[eMaximum]) ||
		(CoarseLimsZ[eMaximum] < CoarseLimsZ[eMinimum] && CoarseLimsZ[eMinimum] - 1 != CoarseLimsZ[eMaximum])
		) {
		L_ERROR("Refined region wraps periodically but is not connected which is not supported. Exiting.", GridUtils::logfile);
	}

#ifdef L_INIT_VERBOSE
	*GridUtils::logfile << "Local Coarse Lims are " <<
		CoarseLimsX[eMinimum] << "-" << CoarseLimsX[eMaximum] << ", " <<
		CoarseLimsY[eMinimum] << "-" << CoarseLimsY[eMaximum] << ", " <<
		CoarseLimsZ[eMinimum] << "-" << CoarseLimsZ[eMaximum] << std::endl;
#endif
}

// ****************************************************************************
/// \brief	Method to initialise the position vector on the grid between the 
///			start and end positions supplied.
///
///			This method can be used for either serial or parallel initialisation
///			as the halo and any wrap around is automatically taken into consideration.
///			As such, the start position can be after the end position and the resulting
///			vector will wrap at the correct point.
///
/// \param	start_pos	position of the first voxel centre in the vector.
///	\param	end_pos		position of the last voxel centre in the vector.
///	\param	dir			direction of the vector.
void GridObj::LBM_initPositionVector(double start_pos, double end_pos, eCartesianDirection dir)
{

#ifdef L_INIT_VERBOSE
	*GridUtils::logfile << "Building position vector for grid level " << level << ", direction " << dir << "...";
#endif

	GridManager *gm = GridManager::getInstance();

	int gm_idx = level + region_number * L_NUM_LEVELS;
	int start_idx = 0, end_idx = 0;
	std::vector<double> *arr = nullptr;
	int req_size = 0;

	switch (dir)
	{
	case eXDirection:
		start_idx = eXMin;
		end_idx = eXMax;
		arr = &XPos;
		req_size = N_lim;
		break;

	case eYDirection:
		start_idx = eYMin;
		end_idx = eYMax;
		arr = &YPos;
		req_size = M_lim;
		break;

	case eZDirection:
		start_idx = eZMin;
		end_idx = eZMax;
		arr = &ZPos;
		req_size = K_lim;
		break;

	default:
		break;
	}

	// Get limits of the grid for wrap around of positions
	double start_lim = gm->global_edges[start_idx][gm_idx];
	double end_lim = gm->global_edges[end_idx][gm_idx];

	// Construct one at a time taking into account periodicity if necessary
	arr->push_back(start_pos);
	do
	{
		if (arr->back() + dh > end_lim) arr->push_back(start_lim + dh / 2);
		else arr->push_back(arr->back() + dh);
	}
	while (arr->size() < req_size);

#ifdef L_INIT_VERBOSE
	*GridUtils::logfile << "Complete." << std::endl;
#endif
}

// ****************************************************************************
/// \brief	Method to initialise label-based solids
void GridObj::LBM_initSolidLab() {

#ifdef L_SOLID_BLOCK_ON
	// Return if not to be put on the current grid
	if (L_BLOCK_ON_GRID_LEV != level || L_BLOCK_ON_GRID_REG != region_number) return;

	// Declarations
	int i, j, k;
	int idx = level + region_number * L_NUM_LEVELS;
	GridManager *gm = GridManager::getInstance();

	/* Check solid block contained on the grid specified in the definitions file.
	 * If not then exit as user has specifed block outside the grid on which they want it 
	 * placed. */

	// Check block placement -- must not be on TL if on a level other than 0
	if	(
		(

		(level == 0) && 

		(L_BLOCK_MAX_X > L_BX || L_BLOCK_MIN_X < 0.0 || L_BLOCK_MAX_Y > L_BY || L_BLOCK_MIN_Y < 0.0 
#if (L_DIMS == 3)
		|| L_BLOCK_MAX_Z > L_BZ || L_BLOCK_MIN_Z < 0.0
#endif
		)

		) || (

		(level != 0) && 

		(
		(
		L_BLOCK_MAX_X > gm->global_edges[eXMax][idx] || 
		L_BLOCK_MIN_X < gm->global_edges[eXMin][idx] || 
		L_BLOCK_MAX_Y > gm->global_edges[eYMax][idx] || 
		L_BLOCK_MIN_Y < gm->global_edges[eYMin][idx] 
#if (L_DIMS == 3)
		||
		L_BLOCK_MAX_Z > gm->global_edges[eZMax][idx] || 
		L_BLOCK_MIN_Z < gm->global_edges[eZMin][idx]
#endif
		)
		||
		(
		GridUtils::isOnTransitionLayer(L_BLOCK_MIN_X, L_BLOCK_MIN_Y, L_BLOCK_MIN_Z, this) ||
		GridUtils::isOnTransitionLayer(L_BLOCK_MAX_X, L_BLOCK_MAX_Y, L_BLOCK_MAX_Z, this)
		)
		)
	
		)
		) {

		// Block outside grid
		L_ERROR("Block is placed outside or on the TL of the selected grid. Exiting.", GridUtils::logfile);
	}


	// Loop over grid
	for (i = 0; i < N_lim; ++i)
	{
		for (j = 0; j < M_lim; ++j)
		{
			for (k = 0; k < K_lim; ++k)
			{

				if (
					XPos[i] <= L_BLOCK_MAX_X && XPos[i] > L_BLOCK_MIN_X &&
					YPos[j] <= L_BLOCK_MAX_Y && YPos[j] > L_BLOCK_MIN_Y
#if (L_DIMS == 3)
					&& ZPos[k] <= L_BLOCK_MAX_Z && ZPos[k] > L_BLOCK_MIN_Z
#endif
					)
				{
					LatTyp(i, j, k, M_lim, K_lim) = eSolid;
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

	// Try to add the solid block
	LBM_initSolidLab();

	// If we need to label the edges...
#if defined L_WALLS_ON || defined L_INLET_ON || defined L_OUTLET_ON || defined L_FREESTREAM_TUNNEL || defined L_UPSTREAM_TUNNEL

	int i, j, k;


#ifdef L_INLET_ON
	// Left hand face only

	// Check for potential singularity in BC
	if (GridUnits::ud2ulbm(L_UMAX,this) == 1 || uref == 1) {
		// Singularity so exit
		L_ERROR("Inlet BC fails with L_UX0 in LBM units = 1, choose something else. Exiting.", GridUtils::logfile);
	}

	// Search position vector to see if left hand wall on this rank
	for (i = 0; i < N_lim; i++ ) {
		if (XPos[i] == dh / 2) {		// Wall found

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
		if (XPos[i] == L_BX - dh / 2) {		// Wall found

			// Label outlet
			for (j = 0; j < M_lim; j++) {
				for (k = 0; k < K_lim; k++) {

#ifdef L_FREESTREAM_TUNNEL
					LatTyp(i,j,k,M_lim,K_lim) = eInlet;
#else
					//LatTyp(i,j,k,M_lim,K_lim) = eOutlet;
					LatTyp(i, j, k, M_lim, K_lim) = eInlet;		// Set to free stream for now
#endif

				}
			}
			break;

		}
	}
#endif

#if ((defined L_WALLS_ON && !defined L_WALLS_ON_2D) || defined L_FREESTREAM_TUNNEL) && (L_DIMS == 3)

	// Search index vector to see if FRONT wall on this rank
	for (k = 0; k < K_lim; k++ ) {
		if (ZPos[k] <= L_WALL_THICKNESS_FRONT) {		// Wall found

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
		if (ZPos[k] >= L_BZ - L_WALL_THICKNESS_BACK) {		// Wall found

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
		if (YPos[j] <= L_WALL_THICKNESS_BOTTOM) {		// Wall found

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
		if (YPos[j] >= L_BY - L_WALL_THICKNESS_TOP) {		// Wall found

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

#endif

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
	
	// Declare indices
	int i, j, k, d;

	// Get edges and indication of TL presence on the refined region from MPIM
	double edges[6];		// Use eCartMinMax to access
	bool TL_present[6];		// Use eCartMinMax to access
	int gm_idx = level + region_number * L_NUM_LEVELS;
	GridManager *gm = GridManager::getInstance();
	for (d = 0; d < 6; ++d)
	{
		edges[d] = gm->global_edges[d][gm_idx];
		TL_present[d] = gm->subgrid_tlayer_key[d][gm_idx - 1];
	}

	// Loop over parent lattice and add "refined" and "TL to lower" labels
	for (i = 0; i < Np_lim; ++i)
	{
		for (j = 0; j < Mp_lim; ++j)
		{
			for (k = 0; k < Kp_lim; ++k)
			{
				if (
					pGrid.XPos[i] > edges[eXMin] && pGrid.XPos[i] < edges[eXMax] &&
					pGrid.YPos[j] > edges[eYMin] && pGrid.YPos[j] < edges[eYMax]
#if (L_DIMS == 3)
					&&
					pGrid.ZPos[k] > edges[eZMin] && pGrid.ZPos[k] < edges[eZMax]
#endif
					)
				{

					// If within single cell width of refined region edge and TL is present then it is TL to lower
					if (
						(pGrid.XPos[i] > edges[eXMin] && pGrid.XPos[i] < edges[eXMin] + pGrid.dh && TL_present[eXMin]) ||
						(pGrid.XPos[i] < edges[eXMax] && pGrid.XPos[i] > edges[eXMax] - pGrid.dh && TL_present[eXMax]) ||
						(pGrid.YPos[j] > edges[eYMin] && pGrid.YPos[j] < edges[eYMin] + pGrid.dh && TL_present[eYMin]) ||
						(pGrid.YPos[j] < edges[eYMax] && pGrid.YPos[j] > edges[eYMax] - pGrid.dh && TL_present[eYMax])
#if (L_DIMS == 3)
						||
						(pGrid.ZPos[k] > edges[eZMin] && pGrid.ZPos[k] < edges[eZMin] + pGrid.dh && TL_present[eZMin]) ||
						(pGrid.ZPos[k] < edges[eZMax] && pGrid.ZPos[k] > edges[eZMax] - pGrid.dh && TL_present[eZMax])
#endif
						)
					{

						// If parent site is fluid site then correct label
						if (pGrid.LatTyp(i, j, k, Mp_lim, Kp_lim) == eFluid) {
							// Change to "TL to lower" label
							pGrid.LatTyp(i, j, k, Mp_lim, Kp_lim) = eTransitionToFiner;
						}

					}
					else if (pGrid.LatTyp(i, j, k, Mp_lim, Kp_lim) == eFluid) {
						// Label it a "refined" site
						pGrid.LatTyp(i, j, k, Mp_lim, Kp_lim) = eRefined;
					}

				}
			}
		}
	}
    

	// Generate grid type matrices for this level //

	// Declare array for parent site indices
	std::vector<int> p;
	eType par_label;
    
    // Loop over sub-grid and add labels based on parent site labels
    for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
#if (L_DIMS == 3)
			for (int k = 0; k < K_lim; k++)
#else
			int k = 0;
#endif
			{
				// Get parent site indices as locals for the current sub-grid local site
				p = GridUtils::getCoarseIndices(
					i, CoarseLimsX[eMinimum], 
					j, CoarseLimsY[eMinimum], 
					k, CoarseLimsZ[eMinimum]
				);
				
				// Get parent site label using local indices
				par_label = pGrid.LatTyp(p[eXDirection], p[eYDirection], p[eZDirection], Mp_lim, Kp_lim);
								
				// If parent is a "TL to lower" then add "TL to upper" label
				if (par_label == eTransitionToFiner) { 
					LatTyp(i, j, k, M_lim, K_lim) = eTransitionToCoarser;
                
					// Else if parent is a "refined" label then label as coarse
				} else if (par_label == eRefined) { 
					LatTyp(i, j, k, M_lim, K_lim) = eFluid;
                
					// Else parent label is some other kind of boundary so copy the
					// label to retain the behaviour onto this grid
				} else { 
					LatTyp(i, j, k, M_lim, K_lim) = par_label;
					
					// If last site to be updated in fine block, change parent label
					// to ensure boundary values are pulled from fine grid
					if (
						(j % 2) != 0 && (i % 2) != 0
#if (L_DIMS == 3)
						&& (k % 2 != 0)
#endif
						) {
						if (pGrid.LatTyp(p[eXDirection], p[eYDirection], p[eZDirection], Mp_lim, Kp_lim) == eSolid || 
							pGrid.LatTyp(p[eXDirection], p[eYDirection], p[eZDirection], Mp_lim, Kp_lim) == eRefinedSolid)
							pGrid.LatTyp(p[eXDirection], p[eYDirection], p[eZDirection], Mp_lim, Kp_lim) = eRefinedSolid;

						else if (pGrid.LatTyp(p[eXDirection], p[eYDirection], p[eZDirection], Mp_lim, Kp_lim) == eSymmetry || 
							pGrid.LatTyp(p[eXDirection], p[eYDirection], p[eZDirection], Mp_lim, Kp_lim) == eRefinedSymmetry)
							pGrid.LatTyp(p[eXDirection], p[eYDirection], p[eZDirection], Mp_lim, Kp_lim) = eRefinedSymmetry;
						
						else if (pGrid.LatTyp(p[eXDirection], p[eYDirection], p[eZDirection], Mp_lim, Kp_lim) == eInlet || 
							pGrid.LatTyp(p[eXDirection], p[eYDirection], p[eZDirection], Mp_lim, Kp_lim) == eRefinedInlet)
							pGrid.LatTyp(p[eXDirection], p[eYDirection], p[eZDirection], Mp_lim, Kp_lim) = eRefinedInlet;

					}
				}
            }
        }
    }

	
	// Try to add the solid block labels
	LBM_initSolidLab();
}

// ***************************************************************************************************
