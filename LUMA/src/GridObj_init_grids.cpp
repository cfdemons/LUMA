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

/* This file contains the routines necessary to initialise the macroscopic quantities and the grids and labels.
*/

#include "../inc/stdafx.h"
#include "../inc/GridObj.h"

using namespace std;

// ****************************************************************************
/// \brief	Method to initialise the lattice velocity.
///
///			If the L_NO_FLOW macro is defined, velocity set to zero everywhere.
///			If L_INIT_VELOCITY_FROM_FILE defined, velocity read from file.
///			Otherwise set from stored profile.
void GridObj::LBM_initVelocity()
{

	// Setup the inlet profile data on this grid
	_LBM_initSetInletProfile();

#ifdef L_INIT_VELOCITY_FROM_FILE

	*GridUtils::logfile << "Loading initial velocity..." << std::endl;

	std::vector<double> x_coord, y_coord, z_coord, ux, uy, uz;
	GridUtils::readVelocityFromFile("./input/initial_velocity.in", x_coord, y_coord, z_coord, ux, uy, uz);
	int gridSize = L_N*L_M*L_K;

	/* Check that the data in the file has the same number of points as the current grid. 
	 * The initial velocity data is copied directly to the cells, not interpolated, 
	 * so the number of points has to match the number of cells. */
	if (x_coord.size() < gridSize) {
		L_ERROR("The initial velocity file has less points than the LBM grid -- change L_BX/L_BY/L_BZ/L_RESOLUTION or change the file initial_velocity.in. Exiting.",
			GridUtils::logfile);
	}
	else if (x_coord.size() > gridSize) {
		L_INFO("WARNING: The initial velocity file has more points than the LBM grid. LUMA will only read as many points as the LBM grid has.",
			GridUtils::logfile);
	}

	// Loop over the data and assign the part that corresponds to the current processor. 
	for (int i = 0; i < gridSize; i++)
	{
		eLocationOnRank loc;
		std::vector<int> indices;
		if (GridUtils::isOnThisRank(x_coord[i], y_coord[i], z_coord[i], &loc, this, &indices))
		{
			u(indices[0], indices[1], indices[2], 0, M_lim, K_lim, L_DIMS) = GridUnits::ud2ulbm(ux[i], this);
			u(indices[0], indices[1], indices[2], 1, M_lim, K_lim, L_DIMS) = GridUnits::ud2ulbm(uy[i], this);
#if (L_DIMS == 3)
			u(indices[0], indices[1], indices[2], 2, M_lim, K_lim, L_DIMS) = GridUnits::ud2ulbm(uz[i], this);
#endif

		}

	}

#else    // L_INIT_VELOCITY_FROM_FILE is not defined

	// Loop over grid
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {

#ifdef L_NO_FLOW
				if (LatTyp(i, j, k, M_lim, K_lim) != eVelocity)
				{
					for (size_t d = 0; d < L_DIMS; d++)
					{
						// No flow case
						u(i,j,k,d,M_lim,K_lim,L_DIMS) = 0.0;
					}
				}

				// Even with no flow still set velocity BCs
				else
#endif
				{

					/* Input velocity is specified by individual vectors for x, y and z which
					* have either been read in from an input file or defined by an expression
					* given in the definitions. */

					// If doing a ramp velocity then initial velocity should be set to ramp at t = 0
					double rampCoefficient = GridUtils::getVelocityRampCoefficient(0.0);

					// Get velocity from profile store
					u(i, j, k, eXDirection, M_lim, K_lim, L_DIMS) = ux_in[j] * rampCoefficient;
					u(i, j, k, eYDirection, M_lim, K_lim, L_DIMS) = uy_in[j] * rampCoefficient;
#if (L_DIMS == 3)
					u(i, j, k, eZDirection, M_lim, K_lim, L_DIMS) = uz_in[j] * rampCoefficient;
#endif
				}

				// Wall sites set to zero
				if (LatTyp(i, j, k, M_lim, K_lim) == eSolid)
				{
					u(i, j, k, eXDirection, M_lim, K_lim, L_DIMS) = 0.0;
					u(i, j, k, eYDirection, M_lim, K_lim, L_DIMS) = 0.0;
#if (L_DIMS == 3)
					u(i, j, k, eZDirection, M_lim, K_lim, L_DIMS) = 0.0;
#endif

				}
			}
		}
	}

#endif	// L_INIT_VELOCITY_FROM_FILE

}


// ****************************************************************************
/// \brief	Method to initialise the lattice density.
void GridObj::LBM_initRho() {

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
	double Lx = gm->global_edges[eXMax][0];
	double Ly = gm->global_edges[eYMax][0];
	double Lz = gm->global_edges[eZMax][0];
	dh = L_COARSE_SITE_WIDTH;

	// Store temporal spacing
	dt = L_TIMESTEP;

	// Store mass conversion
	dm = (L_PHYSICAL_RHO / L_RHOIN) * dh * dh * dh;

	// Gravity in LBM units
	gravity = GridUnits::fd2flbm(L_GRAVITY_FORCE, this);

	// Reference velocity in LBM units -- always 1 in dimensionless units but varies in LBM units
	uref = GridUnits::ud2ulbm(1, this);


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
	

	// L0 lattice site POSITION VECTORS //

#ifdef L_BUILD_FOR_MPI
	
	// Create position vectors using receiver layers as reference positions for wrapping
	// X
	if (MpiManager::getInstance()->dimensions[eXDirection] == 1)
		LBM_initPositionVector(gm->global_edges[eXMin][0] + dh / 2.0, gm->global_edges[eXMax][0] - dh / 2.0, eXDirection);
	else
		LBM_initPositionVector(mpim->recv_layer_pos.X[eLeftMin] + dh / 2.0, mpim->recv_layer_pos.X[eRightMax] - dh / 2.0, eXDirection);

	// Y
	if (MpiManager::getInstance()->dimensions[eYDirection] == 1)
		LBM_initPositionVector(gm->global_edges[eYMin][0] + dh / 2.0, gm->global_edges[eYMax][0] - dh / 2.0, eYDirection);
	else
		LBM_initPositionVector(mpim->recv_layer_pos.Y[eLeftMin] + dh / 2.0, mpim->recv_layer_pos.Y[eRightMax] - dh / 2.0, eYDirection);

	// Z
#if (L_DIMS == 3)
	if (MpiManager::getInstance()->dimensions[eZDirection] == 1)
		LBM_initPositionVector(gm->global_edges[eZMin][0] + dh / 2.0, gm->global_edges[eZMax][0] - dh / 2.0, eZDirection);
	else
		LBM_initPositionVector(mpim->recv_layer_pos.Z[eLeftMin] + dh / 2.0, mpim->recv_layer_pos.Z[eRightMax] - dh / 2.0, eZDirection);
#else
	LBM_initPositionVector(0.0, gm->global_edges[eZMax][0], eZDirection);
#endif

#else
	// When not builiding for MPI positions can be built from global edges directly
	XPos = GridUtils::linspace(dh / 2.0, gm->global_edges[eXMax][0] - dh / 2.0, static_cast<int>(L_N));
	YPos = GridUtils::linspace(dh / 2.0, gm->global_edges[eYMax][0] - dh / 2.0, static_cast<int>(L_M));
#if (L_DIMS == 3)
	ZPos = GridUtils::linspace(dh / 2.0, gm->global_edges[eZMax][0] - dh / 2.0, static_cast<int>(L_K));
#else
	ZPos.push_back(0.0);
#endif

#endif	// L_BUILD_FOR_MPI

#ifdef L_INIT_VERBOSE
	// Write out the position vectors
	std::string msg;
	msg += "XPos: ";
	for (size_t i = 0; i < XPos.size(); ++i) msg += std::to_string(XPos[i]) + " ";
	L_INFO(msg, GridUtils::logfile); msg.clear();
	msg += "YPos: ";
	for (size_t i = 0; i < YPos.size(); ++i) msg += std::to_string(YPos[i]) + " ";
	L_INFO(msg, GridUtils::logfile); msg.clear();
#if (L_DIMS == 3)
	msg += "ZPos: ";
	for (size_t i = 0; i < ZPos.size(); ++i) msg += std::to_string(ZPos[i]) + " ";
	L_INFO(msg, GridUtils::logfile); msg.clear();
#endif
#endif
	

	// Define TYPING MATRICES
	LatTyp.resize(N_lim * M_lim * K_lim);

	// Label as coarse site
	std::fill(LatTyp.begin(), LatTyp.end(), eFluid);

	// Can't use regularised boundaries with D3Q27 because of the corners
#if (defined L_REGULARISED_BOUNDARIES && L_NUM_VELS == 27)
	L_ERROR("Cannot use regularised boundaries with D3Q27 because of the corner treatment. Exiting.", GridUtils::logfile);
#endif

	// Add boundary-specific labels
	LBM_initBoundLab();

	// Initialise L0 MACROSCOPIC quantities

	// Velocity field
	u.resize(N_lim * M_lim * K_lim * L_DIMS);
	LBM_initVelocity();
	
#ifdef L_IBM_ON
	// Set start-of-timestep-velocity
	u_n.resize(N_lim * M_lim * K_lim * L_DIMS);
	u_n = u;
#endif

	// Density field
	rho.resize(N_lim * M_lim * K_lim);
	LBM_initRho();

#if (defined L_GRAVITY_ON || defined L_IBM_ON)
	// Cartesian force vector
	force_xyz.resize(N_lim * M_lim * K_lim * L_DIMS, 0.0);

	// Initialise with gravity
	for (int id = 0; id < N_lim * M_lim * K_lim; ++id)
		force_xyz[L_GRAVITY_DIRECTION + id * L_DIMS] = rho[id] * gravity * refinement_ratio;

	// Lattice force vector
	force_i.resize(N_lim * M_lim * K_lim * L_NUM_VELS, 0.0);
#endif

	// Time averaged quantities
	rho_timeav.resize(N_lim * M_lim * K_lim, 0.0);
	ui_timeav.resize(N_lim * M_lim * K_lim * L_DIMS, 0.0);
	uiuj_timeav.resize(N_lim * M_lim * K_lim * (3 * L_DIMS - 3), 0.0);


	// Initialise L0 POPULATION matrices (f, feq)
	f.resize(N_lim * M_lim * K_lim * L_NUM_VELS);
	feq.resize(N_lim * M_lim * K_lim * L_NUM_VELS);
	fNew.resize(N_lim * M_lim * K_lim * L_NUM_VELS);


	// Loop over grid
	for (int i = 0; i < N_lim; i++)
	{
		for (int j = 0; j < M_lim; j++)
		{
			for (int k = 0; k < K_lim; k++)
			{
				for (int v = 0; v < L_NUM_VELS; v++)
				{
					// Initialise f to feq
					f(i, j, k, v, M_lim, K_lim, L_NUM_VELS) = 
						_LBM_equilibrium_opt(k + j * K_lim + i * M_lim * K_lim, v);

				}
			}
		}
	}
	feq = f; // Make feq = feq too
	fNew = f;


#ifdef L_NU
	nu = GridUnits::nud2nulbm(L_NU, this);
#else
	nu = GridUnits::nud2nulbm(1.0 / static_cast<double>(L_RE), this);
#endif

	// Relaxation frequency on L0
	// Assign relaxation frequency using lattice viscosity
	omega = 1.0 / ( (nu / SQ(cs)) + 0.5 );

	/* Above is valid for L0 only when dh = 1 -- general expression is:
	 * omega = 1 / ( ( (nu * dt) / (pow(cs,2)*pow(dh,2)) ) + .5 );
	 */

	/* Check that the relaxation frequency is within acceptable values. 
	 * Suggest a better value for dt to the user if omega is not within 
	 * acceptable limits. Note that the use of BGKSMAG allows for omega >=2. */
#ifndef L_USE_BGKSMAG
	if (omega >= 2.0)
		L_ERROR("LBM relaxation frequency omega too large. Change L_TIMESTEP or L_RESOLUTION. Exiting.", GridUtils::logfile);
#endif

	// Check if there are incompressibility issues and warn the user if so
	if (uref > (0.17 * cs))
	{
		std::string msg;
		msg += "Reference velocity in LBM units larger than 17% of the speed of sound. ";
		msg += "Compressibility effects may impair the quality of the results. ";
		msg += "Try L_TIMESTEP = " + std::to_string(0.17 * dh * cs) + " or smaller.";
		L_WARN(msg, GridUtils::logfile);
		if (uref >= cs)
		{
			msg.clear();
			msg += "Reference velocity in LBM units equal to or larger than the speed of sound cs. ";
			msg += "Results of the simulation are not valid. Exiting.";
			L_ERROR(msg, GridUtils::logfile);
		}
	}

#ifndef L_BUILD_FOR_MPI
	// Update the writable data in the grid manager
	// When using MPI this is done when building the communicators.
	gm->createWritableDataStore(this);
#endif

#ifdef L_INIT_VERBOSE
	L_INFO("Initialisation Complete.", GridUtils::logfile);
#endif
}
	

// ****************************************************************************
/// \brief	Method to initialise all sub-grid quantities.
/// \param	pGrid	reference to parent grid.
void GridObj::LBM_initSubGrid (GridObj& pGrid)
{

#ifdef L_INIT_VERBOSE
	*GridUtils::logfile << "Initialising sub-grid level " << level << ", region " << region_number << "..." << std::endl;
#endif

	// Get GM instance
	GridManager *gm = GridManager::getInstance();
	int gm_idx = level + region_number * L_NUM_LEVELS;

	// Define scales
	dh = pGrid.dh / 2.0;
	dt = pGrid.dt / 2.0;
	dm = pGrid.dm / 8.0;
	gravity = pGrid.gravity;
	uref = pGrid.uref;
	
	/* Get coarse grid refinement limits as indicies local to the parent grid
	 * on this rank. */
	LBM_initGridToGridMappings(pGrid);
	
	/* Get local grid size of the sub grid based on local ijk limits.
	 * i.e. how much of the parent grid it covers. Volumetric formulation with refinement
	 * factor of 2 makes this easy. */
	int local_size[L_DIMS] =
	{
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
	ZPos.insert( ZPos.begin(), 0.0 ); // 2D default
#endif

	
	// Generate TYPING MATRICES

	// Resize
	LatTyp.resize(N_lim * M_lim * K_lim);

	// Default labelling of coarse
	std::fill(LatTyp.begin(), LatTyp.end(), eFluid);
	
	// Call refined labelling routine passing parent grid
	LBM_initRefinedLab(pGrid);

	
	// Assign MACROSCOPIC quantities

	// Velocity
	u.resize(N_lim * M_lim * K_lim * L_DIMS);
	LBM_initVelocity();

	// Set start-of-timestep-velocity
#ifdef L_IBM_ON
	u_n.resize(N_lim * M_lim * K_lim * L_DIMS);
	u_n = u;
#endif

	// Density
	rho.resize(N_lim * M_lim * K_lim);
	LBM_initRho();


#if (defined L_GRAVITY_ON || defined L_IBM_ON)

	// Cartesian force vector
	force_xyz.resize(N_lim * M_lim * K_lim * L_DIMS, 0.0);

	// Initialise with gravity
	for (int id = 0; id < N_lim * M_lim * K_lim; ++id)
		force_xyz[L_GRAVITY_DIRECTION + id * L_DIMS] = rho[id] * gravity * refinement_ratio;

	// Lattice force vector
	force_i.resize(N_lim * M_lim * K_lim * L_NUM_VELS, 0.0);

#endif

	// Time averaged quantities
#ifdef L_COMPUTE_TIME_AVERAGED_QUANTITIES
	rho_timeav.resize(N_lim * M_lim * K_lim, 0.0);
	ui_timeav.resize(N_lim * M_lim * K_lim * L_DIMS, 0.0);
	uiuj_timeav.resize(N_lim * M_lim * K_lim * (3 * L_DIMS - 3), 0.0);
#endif


	// Generate POPULATION MATRICES for lower levels
	// Resize
	f.resize(N_lim * M_lim * K_lim * L_NUM_VELS);
	feq.resize(N_lim * M_lim * K_lim * L_NUM_VELS);
	fNew.resize(N_lim * M_lim * K_lim * L_NUM_VELS);


	// Loop over grid
	for (int i = 0; i < N_lim; ++i)
	{
		for (int j = 0; j < M_lim; ++j)
		{
			for (int k = 0; k < K_lim; ++k)
			{
				for (int v = 0; v < L_NUM_VELS; v++)
				{
					
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
	 * fine and coarse grid cells during multi-grid operations. Therefore, we
	 * must only store local values relevant to the grid on the rank to ensure
	 * mapping is correct.
	 *
	 * When not using MPI, these can be computed from the definitions file values
	 * which specify the extent of the refinement at that level.
	 * However, when using MPI, the edges of the refined grid might not be on this
	 * rank so we must round the coarse limits to the edge of the parent grid so
	 * the correct offset is supplied to the mapping routine.
	 *
	 * When dealing with sub-grids embedded in walls, it is more complicated. If the
	 * sub-grid starts on a max receiver layer which is also a periodic overlap then
	 * we need to make sure the limits are set properly. Likewise if the sub-grid ends
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
		*GridUtils::logfile << "X: Starts on this rank @ " << pGrid.XPos[position] << ", local idx " << position << ", Loc = " << loc << std::endl;
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
		*GridUtils::logfile << "X: Ends on this rank @ " << pGrid.XPos[position] << ", local idx " << position << ", Loc = " << loc << std::endl;
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

	// TODO: Possible that this can be set incorrectly if the grid has no core and grid is coalesced incorrectly
//	if (GridUtils::isOnThisRank(subgrid_start_voxel_centre, eXDirection, &loc, &pGrid, &position) &&
//		GridUtils::isOnThisRank(subgrid_end_voxel_centre, eXDirection, &loc, &pGrid, &position2))
//	{
//
//		// Check sites are next to each other in space
//		if (abs(pGrid.XPos[position] - pGrid.XPos[position2]) - pGrid.dh < L_SMALL_NUMBER
//			|| gm->periodic_flags[eXDirection][gm_idx])
//		{
//			// Set as if middle of complete block
//			CoarseLimsX[eMinimum] = 0;
//			CoarseLimsX[eMaximum] = pGrid.N_lim - 1;
//
//			// Tell mpim that TL doesn't exist so HDF5 writer does not try to exclude valid sites
//			gm->subgrid_tlayer_key[eXMin][gm_idx - 1] = false;
//			gm->subgrid_tlayer_key[eXMax][gm_idx - 1] = false;
//
//#ifdef L_INIT_VERBOSE
//			*GridUtils::logfile << "X: Periodically wrapped and joined" << std::endl;
//#endif
//		}
//	}

	// If they do not start and end on the same process we still need to update the TL keys 
	if (gm->periodic_flags[eXDirection][gm_idx])
	{
		gm->subgrid_tlayer_key[eXMin][gm_idx - 1] = false;
		gm->subgrid_tlayer_key[eXMax][gm_idx - 1] = false;
	}


	// Y //
	subgrid_start_voxel_centre = gm->global_edges[eYMin][gm_idx] + dh / 2.0;
	subgrid_end_voxel_centre = gm->global_edges[eYMax][gm_idx] - dh / 2.0;
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
		*GridUtils::logfile << "Y: Starts on this rank @ " << pGrid.YPos[position] << ", local idx " << position << ", Loc = " << loc << std::endl;
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
		*GridUtils::logfile << "Y: Ends on this rank @ " << pGrid.YPos[position] << ", local idx " << position << ", Loc = " << loc << std::endl;
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

	/*if (GridUtils::isOnThisRank(subgrid_start_voxel_centre, eYDirection, &loc, &pGrid, &position) &&
		GridUtils::isOnThisRank(subgrid_end_voxel_centre, eYDirection, &loc, &pGrid, &position2))
	{
		if (abs(pGrid.YPos[position] - pGrid.YPos[position2]) - pGrid.dh < L_SMALL_NUMBER
			|| gm->periodic_flags[eYDirection][gm_idx])
		{
			CoarseLimsY[eMinimum] = 0;
			CoarseLimsY[eMaximum] = pGrid.M_lim - 1;
			gm->subgrid_tlayer_key[eYMin][gm_idx - 1] = false;
			gm->subgrid_tlayer_key[eYMax][gm_idx - 1] = false;

#ifdef L_INIT_VERBOSE
			*GridUtils::logfile << "Y: Periodically wrapped and joined" << std::endl;
#endif
		}
	}*/
	if (gm->periodic_flags[eYDirection][gm_idx])
	{
		gm->subgrid_tlayer_key[eYMin][gm_idx - 1] = false;
		gm->subgrid_tlayer_key[eYMax][gm_idx - 1] = false;
	}


#if (L_DIMS == 3)
	// Z //
	subgrid_start_voxel_centre = gm->global_edges[eZMin][gm_idx] + dh / 2.0;
	subgrid_end_voxel_centre = gm->global_edges[eZMax][gm_idx] - dh / 2.0;
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
		*GridUtils::logfile << "Z: Starts on this rank @ " << pGrid.ZPos[position] << ", local idx " << position << ", Loc = " << loc << std::endl;
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
		*GridUtils::logfile << "Z: Ends on this rank @ " << pGrid.ZPos[position] << ", local idx " << position << ", Loc = " << loc << std::endl;
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

	/*if (GridUtils::isOnThisRank(subgrid_start_voxel_centre, eZDirection, &loc, &pGrid, &position) &&
		GridUtils::isOnThisRank(subgrid_end_voxel_centre, eZDirection, &loc, &pGrid, &position2))
	{
		if (abs(pGrid.ZPos[position] - pGrid.ZPos[position2]) - pGrid.dh < L_SMALL_NUMBER
			|| gm->periodic_flags[eZDirection][gm_idx])
		{
			CoarseLimsZ[eMinimum] = 0;
			CoarseLimsZ[eMaximum] = pGrid.K_lim - 1;
			gm->subgrid_tlayer_key[eZMin][gm_idx - 1] = false;
			gm->subgrid_tlayer_key[eZMax][gm_idx - 1] = false;

#ifdef L_INIT_VERBOSE
			*GridUtils::logfile << "Z: Periodically wrapped and joined" << std::endl;
#endif
		}
	}*/
	if (gm->periodic_flags[eZDirection][gm_idx])
	{
		gm->subgrid_tlayer_key[eZMin][gm_idx - 1] = false;
		gm->subgrid_tlayer_key[eZMax][gm_idx - 1] = false;
	}
#else
	// Reset the refined region z-limits if only 2D
	for (int i = 0; i < 2; i++) {
		CoarseLimsZ[i] = 0;
	}
#endif

	/* If sub-grid wraps periodically we do not support it wrapping back round to the same rank
	* again as this will confuse the mapping function so we error here. */

	std::string dir;
	if (CoarseLimsX[eMaximum] < CoarseLimsX[eMinimum] && CoarseLimsX[eMinimum] - 1 != CoarseLimsX[eMaximum]) dir = "X";
	if (CoarseLimsY[eMaximum] < CoarseLimsY[eMinimum] && CoarseLimsY[eMinimum] - 1 != CoarseLimsY[eMaximum]) dir = "Y";
	if (CoarseLimsZ[eMaximum] < CoarseLimsZ[eMinimum] && CoarseLimsZ[eMinimum] - 1 != CoarseLimsZ[eMaximum]) dir = "Z";
	if (dir != "")
		L_ERROR("Refined region wraps periodically in " + dir + "-direction but is not connected which is not supported. Exiting.", GridUtils::logfile);

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
	*GridUtils::logfile << "Building position vector for grid level " << level << ", region " << region_number << ", direction " << dir << "..." << std::endl;
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

#ifdef L_INIT_VERBOSE
	L_INFO("Required Size = " + std::to_string(req_size) + 
		", Start limit = " + std::to_string(start_lim) + 
		", Start position = " + std::to_string(start_pos) +
		", End limit = " + std::to_string(end_lim) +
		", End position = " + std::to_string(end_pos),
		GridUtils::logfile);
#endif

	// Construct one at a time taking into account periodicity if necessary
	arr->push_back(start_pos);
	do
	{
		if (arr->back() + dh > end_lim) arr->push_back(start_lim + dh / 2.0);
		else arr->push_back(arr->back() + dh);
	}
	while (arr->size() < req_size);

#ifdef L_INIT_VERBOSE
	*GridUtils::logfile << "Complete." << std::endl;

	// Write out the position vectors
	if (dir == eXDirection) *GridUtils::logfile << "XPos: ";
	else if (dir == eYDirection) *GridUtils::logfile << "YPos: ";
	else if (dir == eZDirection) *GridUtils::logfile << "ZPos: ";
	for (size_t i = 0; i < arr->size(); ++i)
		*GridUtils::logfile << std::to_string((*arr)[i]) << " ";
	*GridUtils::logfile << "...Complete." << std::endl;
#endif

}

// ****************************************************************************
/// \brief	Method to initialise wall and object labels on L0.
///
///			The virtual wind tunnel definitions are implemented by this method.
void GridObj::LBM_initBoundLab ( )
{
	// Declarations
	int i, j, k;

	// LEFT WALL //

	// Search position vector to see if left hand wall on this rank
	for (i = 0; i < N_lim; i++)
	{
		// Wall found
		if (XPos[i] <= L_WALL_THICKNESS_LEFT)
		{
			// Label boundary
			for (j = 0; j < M_lim; j++)
			{
				for (k = 0; k < K_lim; k++)
				{
					LatTyp(i, j, k, M_lim, K_lim) = 
						LBM_setBCPrecedence(LatTyp(i, j, k, M_lim, K_lim), L_WALL_LEFT);
				}
			}
		}
	}

	// RIGHT WALL //

	// Search index vector to see if right hand wall on this rank
	for (i = 0; i < N_lim; i++)
	{
		if (XPos[i] >= GridManager::getInstance()->global_edges[eXMax][0] - L_WALL_THICKNESS_RIGHT)
		{
			// Label boundary
			for (j = 0; j < M_lim; j++)
			{
				for (k = 0; k < K_lim; k++)
				{
					LatTyp(i, j, k, M_lim, K_lim) = 
						LBM_setBCPrecedence(LatTyp(i, j, k, M_lim, K_lim), L_WALL_RIGHT);
				}
			}
		}
	}

#if (L_DIMS == 3)

	// FRONT WALL //

	for (k = 0; k < K_lim; k++)
	{
		if (ZPos[k] <= L_WALL_THICKNESS_FRONT)
		{
			for (i = 0; i < N_lim; i++)
			{
				for (j = 0; j < M_lim; j++)
				{
					LatTyp(i, j, k, M_lim, K_lim) = 
						LBM_setBCPrecedence(LatTyp(i, j, k, M_lim, K_lim), L_WALL_FRONT);
				}
			}
		}
	}

	// BACK WALL //

	for (k = 0; k < K_lim; k++)
	{
		if (ZPos[k] >= GridManager::getInstance()->global_edges[eZMax][0] - L_WALL_THICKNESS_BACK)
		{
			for (i = 0; i < N_lim; i++)
			{
				for (j = 0; j < M_lim; j++)
				{
					LatTyp(i, j, k, M_lim, K_lim) = 
						LBM_setBCPrecedence(LatTyp(i, j, k, M_lim, K_lim), L_WALL_BACK);
				}
			}
		}
	}
#endif

	// BOTTOM WALL //

	for (j = 0; j < M_lim; j++)
	{
		if (YPos[j] <= L_WALL_THICKNESS_BOTTOM)
		{
			for (i = 0; i < N_lim; i++)
			{
				for (k = 0; k < K_lim; k++)
				{
					LatTyp(i, j, k, M_lim, K_lim) = 
						LBM_setBCPrecedence(LatTyp(i, j, k, M_lim, K_lim), L_WALL_BOTTOM);
				}
			}
		}
	}

	// TOP WALL //

	for (j = 0; j < M_lim; j++)
	{
		if (YPos[j] >= GridManager::getInstance()->global_edges[eYMax][0] - L_WALL_THICKNESS_TOP)
		{
			for (i = 0; i < N_lim; i++)
			{
				for (k = 0; k < K_lim; k++)
				{
					LatTyp(i, j, k, M_lim, K_lim) = 
						LBM_setBCPrecedence(LatTyp(i, j, k, M_lim, K_lim), L_WALL_TOP);
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
				// If this parent site is within the bounds of the sub-grid
				if (
					pGrid.XPos[i] > edges[eXMin] && pGrid.XPos[i] < edges[eXMax] &&
					pGrid.YPos[j] > edges[eYMin] && pGrid.YPos[j] < edges[eYMax]
#if (L_DIMS == 3)
					&&
					pGrid.ZPos[k] > edges[eZMin] && pGrid.ZPos[k] < edges[eZMax]
#endif
					)
				{

					// If within single cell width of sub-grid edge and TL is present then it is TL to lower
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
				if (par_label == eTransitionToFiner)
				{ 
					LatTyp(i, j, k, M_lim, K_lim) = eTransitionToCoarser;
				}

				// Else if parent is a "refined" label then label as coarse
				else if (par_label == eRefined)
				{ 
					LatTyp(i, j, k, M_lim, K_lim) = eFluid;
				}

				// Else parent label is some other kind of boundary so copy the
				// label to retain the behaviour onto this grid
				else
				{ 
					LatTyp(i, j, k, M_lim, K_lim) = par_label;
				}
            }
        }
    }

}

// ****************************************************************************
/// \brief	Method to import an input profile from a file.
///
///			Input data may be over- or under-sampled but it must span the 
///			physical dimensions of the inlet otherwise the software does
///			not known how to scale the data to fit. Inlet profile is always
///			assumed to be oriented vertically (y-direction) and data should be
///			specified in dimensionless units.
void GridObj::_LBM_initGetInletProfileFromFile()
{

	// Declare index
	size_t j;
	size_t i;
	double y;

	// Indicate to log
	*GridUtils::logfile << "Loading inlet profile..." << std::endl;

	IVector<double> xbuffer, ybuffer, zbuffer, uxbuffer, uybuffer, uzbuffer;
	GridUtils::readVelocityFromFile("./input/inlet_profile.in", xbuffer, ybuffer, zbuffer, uxbuffer, uybuffer, uzbuffer);

	// Loop over site positions (for left hand inlet, y positions)
	for (j = 0; j < M_lim; j++)
	{

		// Get Y-position
		y = YPos[j];

		// Find values either side of desired
		for (i = 0; i < ybuffer.size(); i++)
		{

			// Bottom point
			if (i == 0 && ybuffer[i] > y)
			{
				// Extrapolation
				ux_in[j] = -(uxbuffer[i + 1] - uxbuffer[i]) + uxbuffer[i];
				uy_in[j] = -(uybuffer[i + 1] - uybuffer[i]) + uybuffer[i];
				uz_in[j] = -(uzbuffer[i + 1] - uzbuffer[i]) + uzbuffer[i];
				break;

			}
			// Top point
			else if (i == ybuffer.size() - 1 && ybuffer[i] < y)
			{
				// Extrapolation
				ux_in[j] = (uxbuffer[i] - uxbuffer[i - 1]) + uxbuffer[i];
				uy_in[j] = (uybuffer[i] - uybuffer[i - 1]) + uybuffer[i];
				uz_in[j] = (uzbuffer[i] - uzbuffer[i - 1]) + uzbuffer[i];
				break;

			}
			// Any middle point
			else if (ybuffer[i] < y && ybuffer[i + 1] > y)
			{
				// Interpolation
				ux_in[j] = ((uxbuffer[i + 1] - uxbuffer[i]) / (ybuffer[i + 1] - ybuffer[i])) * (y - ybuffer[i]) + uxbuffer[i];
				uy_in[j] = ((uybuffer[i + 1] - uybuffer[i]) / (ybuffer[i + 1] - ybuffer[i])) * (y - ybuffer[i]) + uybuffer[i];
				uz_in[j] = ((uzbuffer[i + 1] - uzbuffer[i]) / (ybuffer[i + 1] - ybuffer[i])) * (y - ybuffer[i]) + uzbuffer[i];
				break;

			}
			// Exactly on a point
			else if (ybuffer[i] == y)
			{
				// Copy
				ux_in[j] = uxbuffer[i];
				uy_in[j] = uybuffer[i];
				uz_in[j] = uzbuffer[i];
				break;

			}
			// Not near enough to this point
			else
			{
				continue;
			}
		}
	}

	// Check for data read fail
	if (ux_in.size() == 0 || uy_in.size() == 0 || uz_in.size() == 0)
	{
		// No data read in
		L_ERROR("Failed to read in inlet profile data. Exiting.", GridUtils::logfile);
	}

}

// ****************************************************************************
/// \brief	Method to initialise the velocity profile used for velocity BCs.
///
///			By default sets the profile to be based on the macro defined velocity 
///			components. If required will initialise from a file or as an
///			analytical parabola.
void GridObj::_LBM_initSetInletProfile()
{
	// Resize vectors
	ux_in.resize(M_lim);
	uy_in.resize(M_lim);
	uz_in.resize(M_lim);
	
#ifdef L_USE_INLET_PROFILE

	// Get from file
	_LBM_initGetInletProfileFromFile();

#elif defined L_PARABOLIC_INLET
	// Specify as an anlytical parabola
	double b = GridManager::getInstance()->global_edges[eYMax][0] - L_WALL_THICKNESS_TOP;
	double p = (b + L_WALL_THICKNESS_BOTTOM) / 2.0;
	double q = b - p;

	// Loop over site positions (for left hand inlet, y positions)
	for (int j = 0; j < M_lim; j++)
	{
		// Set the inlet velocity profile values
		ux_in[j] = GridUnits::ud2ulbm(1.5 * L_UX0, this) * (1.0 - pow((YPos[j] - p) / q, 2.0));
		uy_in[j] = 0.0;
		uz_in[j] = 0.0;
	}

#else

	// Otherwise, set from basic definitions
	for (int j = 0; j < M_lim; j++)
	{
		ux_in[j] = GridUnits::ud2ulbm(L_UX0, this);
		uy_in[j] = GridUnits::ud2ulbm(L_UY0, this);
		uz_in[j] = GridUnits::ud2ulbm(L_UZ0, this);
	}

#endif
}

// *****************************************************************************
/// \brief	Used to preserve the precedence of BCs during labelling.
///
///			This is necessary in particular for corners where types of BCs meet.
///			This method ensures that velocity BCs are always preserved over slip 
///			conditions for example.
///
///	\param	currentBC	BC type already assigned to given site.
/// \param	desiredBC	new BC type for a given site.
/// \returns			permitted BC type for given site.
eType GridObj::LBM_setBCPrecedence(eType currentBC, eType desiredBC)
{
	if (currentBC == eSolid || desiredBC == eSolid) return eSolid;
	else if (currentBC == eVelocity) return eVelocity;
	else return desiredBC;
}

// ***************************************************************************************************
