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
#define LUMA_VERSION "1.6.0-alpha"


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
#define L_SHOW_TIME_TO_COMPLETE		///< Write the estimated time to completion to the terminal


/*
*******************************************************************************
************************* Global configuration data ***************************
*******************************************************************************
*/

// Using MPI?
//#define L_BUILD_FOR_MPI			///< Enable MPI features in build

// Output Options
#define L_OUT_EVERY 500000			///< How many timesteps before whole grid output
#define L_OUT_EVERY_FORCES 1000		///< Specific output frequency of body forces
#define L_OUTPUT_PRECISION 6		///< Precision of output (for text writers)

// Types of output
//#define L_IO_LITE					///< ASCII dump on output
//#define L_HDF5_OUTPUT				///< HDF5 dump on output
//#define L_LD_OUT					///< Write out lift and drag (all bodies)
//#define L_IO_FGA                  ///< Write the components of the macroscopic velocity in a .fga file. (To be used in Unreal Engine 4).
//#define L_COMPUTE_TIME_AVERAGED_QUANTITIES

// High frequency output options
#define L_PROBE_OUTPUT							///< Turn on probe output
#define L_PROBE_OUT_FREQ 100						///< Write out frequency of probe output
#define L_PROBE_NUM_X 2							///< Number of probes in X direction
#define L_PROBE_NUM_Y 5							///< Number of probes in Y direction
#define L_PROBE_NUM_Z 1							///< Number of probes in Z direction
#define L_PROBE_MIN_X 0.5						///< Start position of probe array in X direction
#define L_PROBE_MIN_Y (0.4 + (1.0 / L_RESOLUTION))	///< Start position of probe array in Y direction
#define L_PROBE_MIN_Z 0.0						///< Start position of probe array in Z direction
#define L_PROBE_MAX_X 1.5						///< End position of probe array in X direction
#define L_PROBE_MAX_Y (1.6 + (1.0 / L_RESOLUTION))	///< End position of probe array in Y direction
#define L_PROBE_MAX_Z 0.0						///< End position of probe array in Z direction

// Gravity
#define L_GRAVITY_ON						///< Turn on gravity force
/// Expression for the gravity force in dimensionless units
#define L_GRAVITY_FORCE 1.59421e-2
#define L_GRAVITY_DIRECTION eXDirection		///< Gravity direction (specify using enumeration)

// Initialisation
#define L_NO_FLOW							///< Initialise the domain with no flow
//#define L_INIT_VELOCITY_FROM_FILE			///< Read initial velocity from file
//#define L_RESTARTING						///< Initialise the GridObj with quantities read from a restart file
#define L_RESTART_OUT_FREQ 100000				///< Frequency of write out of restart file

// LBM configuration
//#define L_USE_KBC_COLLISION					///< Use KBC collision operator instead of LBGK by default
//#define L_USE_BGKSMAG
#define L_CSMAG 0.1


/*
*******************************************************************************
******************************** Time data ************************************
*******************************************************************************
*/

#define L_TOTAL_TIMESTEPS 1000		///< Number of time steps to run simulation for


/*
*******************************************************************************
**************************** Domain Dimensions ********************************
*******************************************************************************
*/

// MPI Data
#define L_MPI_XCORES 2		///< Number of MPI ranks to divide domain into in X direction
#define L_MPI_YCORES 2		///< Number of MPI ranks to divide domain into in Y direction
#define L_MPI_ZCORES 2		///< Number of MPI ranks to divide domain into in Z direction.

// Decomposition strategy
#define L_MPI_SMART_DECOMPOSE		///< Use smart decomposition to improve load balancing
#define L_MPI_SD_MAX_ITER 1000		///< Max number of iterations to be used for samrt decomposition algorithm

// Topology report
//#define L_MPI_TOPOLOGY_REPORT		///< Have the MPI Manager report on different combinations of X Y Z cores
#define L_MPI_TOP_XCORES 12			///< Max number of X MPI ranks to use for the topology report
#define L_MPI_TOP_YCORES 12			///< Max number of Y MPI ranks to use for the topology report
#define L_MPI_TOP_ZCORES 12			///< Max number of Z MPI ranks to use for the topology report

/*
*******************************************************************************
****************************** Physical Data **********************************
*******************************************************************************
*/

// Lattice properties
#define L_DIMS 3			///< Number of dimensions to the problem
#define L_RESOLUTION 31		///< Number of coarse lattice sites per unit length
#define L_TIMESTEP 9.32259e-4	///< The timestep in non-dimensional units

// Non-dimensional domain dimensions
#define L_BX 2.0		///< End of domain in X (non-dimensional units)
#define L_BY (2.0 + (2.0 * L_COARSE_SITE_THICKNESS))		///< End of domain in Y (non-dimensional units)
#define L_BZ 2.0		///< End of domain in Z (non-dimensional units)

// Physical velocity
#define L_PHYSICAL_U 0.2		///< Reference velocity of the real fluid to model [m/s]


/*
*******************************************************************************
******************************** Fluid Data ***********************************
*******************************************************************************
*/

// Fluid data in lattice units
//#define L_USE_INLET_PROFILE	///< Use an inlet profile
//#define L_PARABOLIC_INLET		///< Use analytic expression for inlet profile - if not then ASCII file is read (requires L_USE_INLET_PROFILE)

// If not using an inlet profile, specify values or expressions here
#define L_UX0 1.0			///< Initial/inlet x-velocity
#define L_UY0 0.0			///< Initial/inlet y-velocity
#define L_UZ0 0.0			///< Initial/inlet z-velocity

#define L_RHOIN 1			///< Initial density. In lattice units. 
//#define L_NU 0            ///< Dimensionless kinematic viscosity L_NU = 1/Re. Comment it to use L_RE instead.  
#define L_RE 3294.01			///< Desired Reynolds number
//#define L_REYNOLDS_RAMP 1000	///< Defines over how many time steps to ramp the Reynolds number


/*
*******************************************************************************
****************************** Object Management ******************************
*******************************************************************************
*/

// General //
//#define L_GEOMETRY_FILE					///< If defined LUMA will read for geometry config file
//#define L_VTK_BODY_WRITE				///< Write out the bodies to a VTK file

// IBM //
//#define L_IBM_ON						///< Turn on IBM
//#define L_STOP_EPSILON_RECOMPUTE		///< Prevent recomputing of epsilon in an attempt to save time

/*
*******************************************************************************
********************************** Wall Data **********************************
*******************************************************************************
*/

// BC types (unspecified is periodic)
#define L_WALL_LEFT		eFluid		///< BC used on the left of the domain
#define L_WALL_RIGHT	eFluid		///< BC used on the right of the domain
#define L_WALL_BOTTOM	eSolid			///< BC used on the bottom of the domain
#define L_WALL_TOP		eSolid			///< BC used on the top of the domain
#define L_WALL_FRONT	eFluid			///< BC used on the front of the domain
#define L_WALL_BACK		eFluid			///< BC used on the bottom of the domain

// BC qualifiers
//#define L_REGULARISED_BOUNDARIES	///< Specify the velocity and pressure BCs to be regularised (Latt & Chopard)
//#define L_OUTLET_EXTRAPOLATED		///< Specifies that the outlet BC extrapolates information from the domain
//#define L_VELOCITY_RAMP 100		///< Defines the number of timesteps over which to ramp up the inlet velocity

// General
#define L_WALL_THICKNESS_BOTTOM (1.0 * L_COARSE_SITE_THICKNESS)	///< Thickness of wall
#define L_WALL_THICKNESS_TOP L_COARSE_SITE_THICKNESS			///< Thickness of top wall
#define L_WALL_THICKNESS_LEFT (1.0 * L_COARSE_SITE_THICKNESS)	///< Thickness of left wall
#define L_WALL_THICKNESS_RIGHT L_COARSE_SITE_THICKNESS			///< Thickness of right wall
#define L_WALL_THICKNESS_FRONT L_COARSE_SITE_THICKNESS			///< Thickness of front (3D) wall
#define L_WALL_THICKNESS_BACK L_COARSE_SITE_THICKNESS			///< Thickness of back (3D) wall


/*
*******************************************************************************
****************************** Multi-grid Data ********************************
*******************************************************************************
*/

#define L_NUM_LEVELS 0		///< Levels of refinement (0 = coarse grid only)
#define L_NUM_REGIONS 1		///< Number of refined regions (can be arbitrary if L_NUM_LEVELS = 0)
//#define L_AUTO_SUBGRIDS		///< Activate auto sub-grid generation using the padding parameters below

// If you want coincident edges then set to (-2.0 * dh)
#define L_PADDING_X_MIN 0.0		///< Padding between X start of each sub-grid and its child edge
#define L_PADDING_X_MAX 0.0		///< Padding between X end of each sub-grid and its child edge
#define L_PADDING_Y_MIN (-2.0 * dh)		///< Padding between Y start of each sub-grid and its child edge
#define L_PADDING_Y_MAX 0.0		///< Padding between Y end of each sub-grid and its child edge
#define L_PADDING_Z_MIN 0.0		///< Padding between Z start of each sub-grid and its child edge
#define L_PADDING_Z_MAX 0.0		///< Padding between Z end of each sub-grid and its child edge

#if L_NUM_LEVELS != 0
// Position of each refined region

static double cRefStartX[L_NUM_LEVELS][L_NUM_REGIONS] =
{
	{ 0.1 }, { 0.1 }, { 0.5 }, { 0.6 }, { 1.4 }
};
static double cRefEndX[L_NUM_LEVELS][L_NUM_REGIONS] =
{
	{ 3.0 }, { 2.8 }, { 2.7 }, { 2.5 }, { 2.0 }
};
static double cRefStartY[L_NUM_LEVELS][L_NUM_REGIONS] = 
{
	{ 1.0 / static_cast<double>(L_RESOLUTION) }, { 1.0 / static_cast<double>(L_RESOLUTION) },
	{ 1.0 / static_cast<double>(L_RESOLUTION) }, { 1.0 / static_cast<double>(L_RESOLUTION) },	{ 0.4 }
};
static double cRefEndY[L_NUM_LEVELS][L_NUM_REGIONS] = 
{
	{ 1.0 }, { 0.8 }, { 0.65 }, { 0.6 }, { 0.55 }
};
static double cRefStartZ[L_NUM_LEVELS][L_NUM_REGIONS] = 
{
	{ 0.2 }, { 0.4 }, { 0.5 }, { 0.55 }, { 0.575 }
};
static double cRefEndZ[L_NUM_LEVELS][L_NUM_REGIONS] = 
{
	{ 1.4 }, { 1.2 }, { 1.1 }, { 1.05 }, { 1.075 }
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

// Set probes
const static int cNumProbes[3] = { L_PROBE_NUM_X, L_PROBE_NUM_Y, L_PROBE_NUM_Z };
const static double cProbeLimsX[2] = { L_PROBE_MIN_X, L_PROBE_MAX_X };	///< Limits of X plane for array of probes
const static double cProbeLimsY[2] = { L_PROBE_MIN_Y, L_PROBE_MAX_Y };	///< Limits of Y plane for array of probes
const static double cProbeLimsZ[2] = { L_PROBE_MIN_Z, L_PROBE_MAX_Z };	///< Limits of Z plane for array of probes

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
