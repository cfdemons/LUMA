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

#pragma once

#include <vector>
#include "IVector.h"
#include "IBBody.h"
#include <iostream>
#include <fstream>
#include "hdf5luma.h"

/// \enum  eType
/// \brief Lattice typing labels
enum eType
{
	eSolid,					///< Rigid, solid site
	eFluid,					///< Fluid site
	eRefined,				///< Fluid site which is represented on a finer grid
	eTransitionToCoarser,	///< Fluid site coupled to a coarser grid
	eTransitionToFiner,		///< Fluid site coupled to a finer grid
	eBFL,					///< Site containing a BFL marker
	eSymmetry,				///< Symmetry boundary
	eInlet,					///< Inlet boundary
	eOutlet,				///< Outlet boundary
	eRefinedSolid,			///< Rigid, solid site represented on a finer grid
	eRefinedSymmetry,		///< Symmtery boundary represented on a finer grid
	eRefinedInlet			///< Inlet site represented on a finer grid
};

/// \enum  eBCType
/// \brief Flag for indicating which BCs to apply
enum eBCType
{
	eBCAll,				///< Apply all BCs
	eBCSolidSymmetry,	///< Apply just solid and symmetry BCs
	eBCInlet,			///< Apply just inlet BCs
	eBCOutlet,			///< Apply just outlet BCs
	eBCInletOutlet,		///< Apply inlet and outlet BCs
	eBCBFL				///< Apply just BFL BCs
};

/// \brief	Grid class.
///
///			This class represents a grid (lattice) and is capable of owning a 
///			nested hierarchy of child grids.
class GridObj
{

	// Allow MpiManager and ObjectManager objects to access private Grid information as required
	friend class MpiManager;
	friend class ObjectManager;
	friend class GridUtils;

public:

	// Constructors & Destructor //
	GridObj( ); // Default constructor
	GridObj(int level); // Basic grid constructor
	GridObj(int RegionNumber, GridObj& pGrid); // Sub grid constructor with region and reference to parent grid for initialisation
	// MPI L0 constructor with local size and global edges
	GridObj(int level, std::vector<int> local_size, 
		std::vector< std::vector<int> > GlobalLimsInd, 
		std::vector< std::vector<double> > GlobalLimsPos);
	~GridObj( ); // Default destructor


	/*
	***************************************************************************************************************
	********************************************* Member Data *****************************************************
	***************************************************************************************************************
	*/

private :

	/// 1D subgrid array (size = L_NumReg)
	std::vector<GridObj> subGrid;

	// Start and end indices of corresponding coarse level
	// When using MPI these values are local to a particular coarse grid
	int CoarseLimsX[2];		///< Local X indices corresponding to where this grid is locate on parent grid
	int CoarseLimsY[2];		///< Local Y indices corresponding to where this grid is locate on parent grid
	int CoarseLimsZ[2];		///< Local Z indices corresponding to where this grid is locate on parent grid

	// 1D arrays
public :
	std::vector<int> XInd;		///< Vector of global X indices of each site
	std::vector<int> YInd;		///< Vector of global Y indices of each site
	std::vector<int> ZInd;		///< Vector of global Z indices of each site
	std::vector<double> XPos;	///< Vector of global X positions of each site
	std::vector<double> YPos;	///< Vector of global Y positions of each site
	std::vector<double> ZPos;	///< Vector of global Z positions of each site

private :
	// Inlet velocity profile
	std::vector<double> ux_in;	///< Vector of x-component of inlet velocity read from file
	std::vector<double> uy_in;	///< Vector of y-component of inlet velocity read from file
	std::vector<double> uz_in;	///< Vector of z-component of inlet velocity read from file

	// Vector nodal properties
	// Flattened 4D arrays (i,j,k,vel)
	IVector<double> f;				///< Distribution functions
	IVector<double> feq;			///< Equilibrium distribution functions
	IVector<double> u;				///< Macropscopic velocity components
	IVector<double> force_xyz;		///< Macroscopic body force components
	IVector<double> force_i;		///< Mesoscopic body force components

	// Scalar nodal properties
	// Flattened 3D arrays (i,j,k)
	IVector<double> rho;			///< Macroscopic density

	// Grid scalars
	double dx;			///< Physical lattice X spacing
	double dy;			///< Physical lattice Y spacing
	double dz;			///< Physical lattice Z spacing
	int region_number;	///< Region number

	// Time averaged statistics
	IVector<double> rho_timeav;		///< Time-averaged density at each grid point (i,j,k)
	IVector<double> ui_timeav;		///< Time-averaged velocity at each grid point (i,j,k,L_dims)
	IVector<double> uiuj_timeav;	///< Time-averaged velocity products at each grid point (i,j,k,3*L_dims-3)


	// Public data members
public :

	IVector<eType> LatTyp;			///< Flattened 3D array of site labels
	int level;						///< Level in embedded grid hierarchy
	double dt;						///< Physical time step size
	int t;							///< Number of completed iterations on this level
	double nu;						///< Kinematic viscosity (in lattice units)
	double omega;					///< Relaxation frequency

	// Timing variables
	double timeav_mpi_overhead;		///< Time-averaged time of MPI communication
	double timeav_timestep;			///< Time-averaged time of a timestep

	// Local grid sizes
	int N_lim;		///< Local size of grid in X-direction
	int M_lim;		///< Local size of grid in Y-direction
	int K_lim;		///< Local size of grid in Z-direction


	/*
	***************************************************************************************************************
	********************************************* Member Methods **************************************************
	***************************************************************************************************************
	*/

public :

	// Initialisation functions
	void LBM_initVelocity();		// Initialise the velocity field
	void LBM_initRho();				// Initialise the density field
	void LBM_initGrid();			// Non-MPI wrapper for initialiser
	void LBM_initGrid(std::vector<int> local_size,
		std::vector< std::vector<int> > GlobalLimsInd,
		std::vector< std::vector<double> > GlobalLimsPos);		// Initialise top level grid with fields and labels
	void LBM_initSubGrid(GridObj& pGrid);		// Initialise subgrid with all quantities
	void LBM_initBoundLab();					// Initialise labels for walls
	void LBM_initSolidLab();					// Initialise labels for solid objects
	void LBM_initRefinedLab(GridObj& pGrid);	// Initialise labels for refined regions
	void LBM_init_getInletProfile();			// Initialise the store for inlet profile data from file

	// LBM operations
	void LBM_multi(bool IBM_flag);		// Launch the multi-grid kernel
	void LBM_collide();					// Apply collision + 1 overload for equilibrium calculation
	double LBM_collide(int i, int j, int k, int v);
	void LBM_kbcCollide(int i, int j, int k, IVector<double>& f_new);		// KBC collision operator
	void LBM_stream();							// Stream populations
	void LBM_macro();							// Compute macroscopic quantities + 1 overload for single site
	void LBM_macro(int i, int j, int k);
	void LBM_boundary(int bc_type_flag);		// Apply boundary conditions
	void LBM_forcegrid(bool reset_flag);		// Apply a force to the grid points (or simply reset force vectors if flag is true)

	// Boundary operations
	void bc_applyBounceBack(int label, int i, int j, int k);	// Application of HWBB BC
	void bc_applySpecReflect(int label, int i, int j, int k);	// Application of HWSR BC
	void bc_applyRegularised(int label, int i, int j, int k);	// Application of Regaulrised BC
	void bc_applyExtrapolation(int label, int i, int j, int k);			// Application of Extrapolation BC
	void bc_applyBfl(int i, int j, int k);														// Application of BFL BC
	void bc_applyNrbc(int i, int j, int k);														// Application of characteristic NRBC
	void bc_solidSiteReset();																	// Reset all the solid site velocities to zero
	double bc_getWallDensityForRBC(std::vector<double>& ftmp, int normal,
		int i, int j, int k);		// Gets wall density for generalised, regularised velocity BC

	// Multi-grid operations
	void LBM_explode(int RegionNumber);			// Explode populations from coarse to fine
	void LBM_coalesce(int RegionNumber);		// Coalesce populations from fine to coarse
	void LBM_addSubGrid(int RegionNumber);		// Add and initialise subgrid structure for a given region number

	// IO methods
	void io_textout(std::string output_tag);	// Writes out the contents of the class as well as any subgrids to a text file
	void io_restart(bool IO_flag);				// Reads/writes data from/to the global restart file
	void io_probeOutput();						// Output routine for point probes
	void io_lite(double tval, std::string Tag);	// Generic writer to individual files with Tag
	int io_hdf5(double tval);					// HDF5 writer returning integer to indicate success or failure

};

