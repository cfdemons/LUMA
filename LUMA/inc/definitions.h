/*
 * --------------------------------------------------------------
 *
 * ------ Lattice Boltzmann @ The University of Manchester ------
 *
 * -------------------------- L-U-M-A ---------------------------
 *
 *  Copyright (C) The University of Manchester 2017
 *  E-mail contact: info@luma.manchester.ac.uk
 *
 * This software is for academic use only and not available for
 * further distribution commericially or otherwise without written consent.
 *
 */

/// LUMA version
#define LUMA_VERSION "1.5.0"


// Header guard
#ifndef LBM_DEFINITIONS_H
#define LBM_DEFINITIONS_H

// Declarations here
#include <time.h>			// Timing functionality
#include <iostream>			// IO functionality
#include <fstream>			// File functionality
#include <vector>			// Vector template access
#include <iomanip>			// Output precision control
#include <math.h>			// Mathematics
#include <string>			// String template access
#include <mpi.h>			// Enable MPI


/*
*******************************************************************************
**************************** Debugging Options ********************************
*******************************************************************************
*/


//#define L_MEGA_DEBUG				///< Debug F, Feq, Macroscopic all in one file -- Warning: Heavy IO which kills performance
//#define L_INC_RECV_LAYER			///< Flag to include writing out receiver layer sites in MPI builds
//#define L_INIT_VERBOSE			///< Write out initialisation information such as refinement mappings
//#define L_MPI_VERBOSE				///< Write out the buffers used by MPI plus more setup data
//#define L_MPI_WRITE_LOAD_BALANCE	///< Write out the load balancing information based on active cell count
//#define L_IBM_DEBUG				///< Write IBM body and matrix data out to text files
//#define L_IBBODY_TRACER			///< Write out IBBody positions
//#define L_BFL_DEBUG				///< Write out BFL marker positions and Q values out to files
//#define L_CLOUD_DEBUG				///< Write out to a file the cloud that has been read in
//#define L_LOG_TIMINGS				///< Write out the initialisation, time step and mpi timings to an output file
//#define L_HDF_DEBUG				///< Write some HDF5 debugging information
//#define L_TEXTOUT					///< Verbose ASCII output of grid information
//#define L_MOMEX_DEBUG				///< Debug momentum exchange by writing out F contributions verbosely
#define L_SHOW_TIME_TO_COMPLETE		        ///< Write the estimated time to completion to the terminal


/*
*******************************************************************************
************************* Global configuration data ***************************
*******************************************************************************
*/

// Using MPI?
#define L_BUILD_FOR_MPI				///< Enable MPI features in build

// Output Options
#define L_OUT_EVERY 100			///< How many timesteps before whole grid output
#define L_OUT_EVERY_FORCES 20		///< Specific output frequency of body forces
#define L_OUTPUT_PRECISION 8		///< Precision of output (for text writers)

// Types of output
//#define L_IO_LITE					///< ASCII dump on output
#define L_HDF5_OUTPUT				///< HDF5 dump on output
#define L_LD_OUT					///< Write out lift and drag (all bodies)
//#define L_IO_FGA                  ///< Write the components of the macroscopic velocity in a .fga file. (To be used in Unreal Engine 4).
//#define L_COMPUTE_TIME_AVERAGED_QUANTITIES

// High frequency output options
//#define L_PROBE_OUTPUT						///< Turn on probe output
#define L_PROBE_OUT_FREQ 200					///< Write out frequency of probe output
const static int cNumProbes[3] = {3, 3, 3};		///< Number of probes in each direction (x, y, z)
const static double cProbeLimsX[2] = {0.1, 0.2};	///< Limits of X plane for array of probes
const static double cProbeLimsY[2] = {0.1, 0.2};	///< Limits of Y plane for array of probes
const static double cProbeLimsZ[2] = {0.1, 0.2};	///< Limits of Z plane for array of probes


// Gravity
//#define L_GRAVITY_ON						///< Turn on gravity force
/// Expression for the gravity force in dimensionless units
#define L_GRAVITY_FORCE (12.0 / L_RE)
#define L_GRAVITY_DIRECTION eXDirection		///< Gravity direction (specify using enumeration)

// Initialisation
#define L_NO_FLOW							///< Initialise the domain with no flow
//#define L_INIT_VELOCITY_FROM_FILE			///< Read initial velocity from file
//#define L_RESTARTING						///< Initialise the GridObj with quantities read from a restart file
#define L_RESTART_OUT_FREQ 100000			///< Frequency of write out of restart file

// LBM configuration
//#define L_USE_KBC_COLLISION					///< Use KBC collision operator instead of LBGK by default
//#define L_USE_BGKSMAG
#define L_CSMAG 0.07


/*
*******************************************************************************
******************************** Time data ************************************
*******************************************************************************
*/

#define L_TOTAL_TIMESTEPS 60000		///< Number of time steps to run simulation for


/*
*******************************************************************************
**************************** Domain Dimensions ********************************
*******************************************************************************
*/

// MPI Data
#define L_MPI_XCORES 2		///< Number of MPI ranks to divide domain into in X direction
#define L_MPI_YCORES 2		///< Number of MPI ranks to divide domain into in Y direction
/// Number of MPI ranks to divide domain into in Z direction.
#define L_MPI_ZCORES 2

// Balanced decomposition
#define L_MPI_SMART_DECOMPOSE
#define L_MPI_SD_MAX_ITER 1000

/*
*******************************************************************************
****************************** Physical Data **********************************
*******************************************************************************
*/

// Lattice properties
#define L_DIMS 2				///< Number of dimensions to the problem
#define L_RESOLUTION 200			///< Number of coarse lattice sites per unit length
#define L_TIMESTEP 0.000462963			///< The timestep in non-dimensional units

// Non-dimensional domain dimensions
#define L_BX 2.2								///< End of domain in X (non-dimensional units)
#define L_BY (0.41 + 2.0 / L_RESOLUTION)		///< End of domain in Y (non-dimensional units)
#define L_BZ 1.0								///< End of domain in Z (non-dimensional units)

// Physical velocity
#define L_PHYSICAL_U 1.0		///< Reference velocity of the real fluid to model [m/s]

// Reference density
#define L_RHO_REF 1.0


/*
*******************************************************************************
******************************** Fluid Data ***********************************
*******************************************************************************
*/

// Fluid data in lattice units
//#define L_USE_INLET_PROFILE	   	///< Use an inlet profile
#define L_PARABOLIC_INLET	   		///< Use analytic expression for inlet profile - if not then ASCII file is read (requires L_USE_INLET_PROFILE)

// If not using an inlet profile, specify values or expressions here
#define L_UX0 1.0			///< Initial/inlet x-velocity
#define L_UY0 0.0			///< Initial/inlet y-velocity
#define L_UZ0 0.0			///< Initial/inlet z-velocity

#define L_RHOIN 1			///< Initial density. In lattice units. 
//#define L_NU 0            ///< Dimensionless kinematic viscosity L_NU = 1/Re. Comment it to use L_RE instead.
#define L_RE 1000			///< Desired Reynolds number


/*
*******************************************************************************
****************************** Object Management ******************************
*******************************************************************************
*/

// General //
#define L_GEOMETRY_FILE					///< If defined LUMA will read for geometry config file
#define L_VTK_BODY_WRITE				///< Write out the bodies to a VTK file

// IBM //
#define L_IBM_ON						///< Turn on IBM
//#define L_UNIVERSAL_EPSILON_CALC		///< Do universal epsilon calculation (should be used if supports from different bodies overlap)

/*
*******************************************************************************
********************************** Wall Data **********************************
*******************************************************************************
*/

// Virtual Wind Tunnels
//#define L_FREESTREAM_TUNNEL		///< Adds a velocity BC to all faces

// Type of Inlet/Outlet BC (default Forced Equilibrium)
#define L_VELOCITY_REGULARISED	///< Specify the inlet/outlet BC to be a regularised velocity condition (Latt & Chopard)

// Inlet (left-hand wall)
#define L_INLET_ON				///< Turn on inlet boundary

// Outlet (right-hand wall)
#define L_OUTLET_ON				///< Turn on outlet boundary
#define EXTRAPOLATED_OUTLET		///< Extrapolate the velocity from the outlet

// Solids
#define L_WALLS_ON			///< Turn on no-slip walls (default is top, bottom, front, back unless L_WALLS_ON_2D is used)
//#define L_WALLS_ON_2D			///< Limit no-slip walls to top and bottom no-slip walls only
//#define L_WALL_FLOOR_ONLY
#define L_WALL_THICKNESS_BOTTOM L_COARSE_SITE_THICKNESS			///< Thickness of wall
#define L_WALL_THICKNESS_TOP L_COARSE_SITE_THICKNESS			///< Thickness of top wall
#define L_WALL_THICKNESS_FRONT L_COARSE_SITE_THICKNESS			///< Thickness of front (3D) wall
#define L_WALL_THICKNESS_BACK L_COARSE_SITE_THICKNESS			///< Thickness of back (3D) wall


/*
*******************************************************************************
****************************** Multi-grid Data ********************************
*******************************************************************************
*/

#define L_NUM_LEVELS 0		///< Levels of refinement (0 = coarse grid only)
#define L_NUM_REGIONS 0		///< Number of refined regions (can be arbitrary if L_NUM_LEVELS = 0)
#define L_AUTO_SUBGRIDS		///< Activate auto sub-grid generation using the padding parameters below

// If you want coincident edges then set to (-2.0 * dh)
#define L_PADDING_X_MIN 0.1		///< Padding between X start of each sub-grid and its child edge
#define L_PADDING_X_MAX 0.1		///< Padding between X end of each sub-grid and its child edge
#define L_PADDING_Y_MIN 0.1		///< Padding between Y start of each sub-grid and its child edge
#define L_PADDING_Y_MAX 0.1		///< Padding between Y end of each sub-grid and its child edge
#define L_PADDING_Z_MIN 0.1		///< Padding between Z start of each sub-grid and its child edge
#define L_PADDING_Z_MAX 0.1		///< Padding between Z end of each sub-grid and its child edge

#if L_NUM_LEVELS != 0
// Position of each refined region

static double cRefStartX[L_NUM_LEVELS][L_NUM_REGIONS] = {
	{ 2.0 }
};
static double cRefEndX[L_NUM_LEVELS][L_NUM_REGIONS] = {
	{ 7.0 }
};
static double cRefStartY[L_NUM_LEVELS][L_NUM_REGIONS] = {
	{ 1.0 }
};
static double cRefEndY[L_NUM_LEVELS][L_NUM_REGIONS] = {
	{ 9.0 }
};
static double cRefStartZ[L_NUM_LEVELS][L_NUM_REGIONS] = {
	{ 2.5 }
};
static double cRefEndZ[L_NUM_LEVELS][L_NUM_REGIONS] = {
	{ 4.4 }
};

#endif


/*
*******************************************************************************
************************* Clean-up: NOT FOR EDITING ***************************
*******************************************************************************
*/

#define L_N static_cast<int>(L_BX * L_RESOLUTION)
#define L_M static_cast<int>(L_BY * L_RESOLUTION)
#define L_K static_cast<int>(L_BZ * L_RESOLUTION)
#define L_COARSE_SITE_THICKNESS (static_cast<double>(L_BX)/static_cast<double>(L_N))

// Set dependent options
#if (L_DIMS == 3)

	#ifdef L_USE_KBC_COLLISION
		#define L_NUM_VELS 27		///< Number of lattice velocities
	#else
		#define L_NUM_VELS 19		///< Number of lattice velocities
	#endif

	#define L_MPI_DIRS 26	///< Number of MPI directions

#else
	#define L_NUM_VELS 9		// Use D2Q9

	// MPI config to 2D
	#define L_MPI_DIRS 8

	// Set Z limits for 2D
	#undef L_BZ
	#define L_BZ 0

	#undef L_K
	#define L_K 1

	#undef L_MPI_ZCORES
	#define L_MPI_ZCORES 1

	// Set object limits for 2D
	#undef L_BLOCK_MIN_Z
	#define L_BLOCK_MIN_Z 0.0

	#undef L_BLOCK_MAX_Z
	#define L_BLOCK_MAX_Z 0.0

	// Set z inlet velocity
	#undef L_UZ0
	#define L_UZ0 0.0

#endif

#if L_NUM_LEVELS == 0
	// Set region info to default as no refinement
	static double cRefStartX[1][1]	= {0.0};
	static double cRefEndX[1][1]	= {0.0};
	static double cRefStartY[1][1]	= {0.0};
	static double cRefEndY[1][1]	= {0.0};
	static double cRefStartZ[1][1]	= {0.0};
	static double cRefEndZ[1][1]	= {0.0};

	#undef L_NUM_REGIONS
	#define L_NUM_REGIONS 1
#endif

// Clean up for using profiled inlet
#ifdef L_USE_INLET_PROFILE
	#undef L_UX0
	#define L_UX0 ux_in[j]
	#undef L_UY0
	#define L_UY0 uy_in[j]
	#undef L_UZ0
	#define L_UZ0 uz_in[j]
#endif

#endif
