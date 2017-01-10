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

#ifndef GRIDOBJ_H
#define GRIDOBJ_H

#include "stdafx.h"
#include "IVector.h"

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

/// \enum  eIOFlag
/// \brief Flag for indicating write or read action for IO methods
enum eIOFlag
{
	eWrite,				///< Write to file
	eRead,				///< Read from file
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

	/// 1D subgrid array (size = L_NUM_REGIONS)
	std::vector<GridObj> subGrid;

	/// Pointer to parent grid
	GridObj *parentGrid = nullptr;

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
	IVector<double> fNew;			///< Copy of distribution functions
	IVector<double> u;				///< Macropscopic velocity components
	IVector<double> force_xyz;		///< Macroscopic body force components
	IVector<double> force_i;		///< Mesoscopic body force components

	// Scalar nodal properties
	// Flattened 3D arrays (i,j,k)
	IVector<double> rho;			///< Macroscopic density

	// Time averaged statistics
	IVector<double> rho_timeav;		///< Time-averaged density at each grid point (i,j,k)
	IVector<double> ui_timeav;		///< Time-averaged velocity at each grid point (i,j,k,L_DIMS)
	IVector<double> uiuj_timeav;	///< Time-averaged velocity products at each grid point (i,j,k,3*L_DIMS-3)

	// Grid scale parameter
	double refinement_ratio;	///< Equivalent to (1 / pow(2, level))

	// Public data members
public :

	IVector<eType> LatTyp;			///< Flattened 3D array of site labels

	// Grid Scalars
	double dx;						///< Physical lattice X spacing
	double dy;						///< Physical lattice Y spacing
	double dz;						///< Physical lattice Z spacing
	int region_number;				///< Region number
	int level;						///< Level in embedded grid hierarchy
	double dt;						///< Physical time step size
	int t;							///< Number of completed iterations on this level
	double nu;						///< Kinematic viscosity (in lattice units)
	double omega;					///< Relaxation frequency

	// Timing variables
	double timeav_mpi_overhead;		///< Time-averaged time of MPI communication
	double timeav_timestep;			///< Time-averaged time of a timestep

	// Local grid sizes
	int N_lim;			///< Local size of grid in X-direction
	int M_lim;			///< Local size of grid in Y-direction
	int K_lim;			///< Local size of grid in Z-direction
	double XOrigin;		///< Global position of grid left edge
	double YOrigin;		///< Global position of grid bottom edge
	double ZOrigin;		///< Global position of grid front edge


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
	DEPRECATED void LBM_multi(bool ibmFlag);		// Launch the multi-grid kernel (DEPRECATED VERSION)
	DEPRECATED void LBM_multi();					// Launch the multi-grid kernel
	DEPRECATED void LBM_collide();					// Apply collision + 1 overload for equilibrium calculation
	double LBM_collide(int i, int j, int k, int v);
	void LBM_kbcCollide(int i, int j, int k, IVector<double>& f_new);		// KBC collision operator
	DEPRECATED void LBM_stream();					// Stream populations
	DEPRECATED void LBM_macro();					// Compute macroscopic quantities + 1 overload for single site
	void LBM_macro(int i, int j, int k);
	DEPRECATED void LBM_boundary(int bc_type_flag);	// Apply boundary conditions
	DEPRECATED void LBM_forceGrid();				// Apply a force to the grid points
	void LBM_resetForces();							// Resets the force vectors on the grid

	// Boundary operations
	void bc_applyBounceBack(int label, int i, int j, int k);	// Application of HWBB BC
	void bc_applySpecReflect(int label, int i, int j, int k);	// Application of HWSR BC
	void bc_applyRegularised(int label, int i, int j, int k);	// Application of Regaulrised BC
	void bc_applyExtrapolation(int label, int i, int j, int k);	// Application of Extrapolation BC
	void bc_applyBfl(int i, int j, int k);						// Application of BFL BC
	void bc_applyNrbc(int i, int j, int k);						// Application of characteristic NRBC
	DEPRECATED void bc_solidSiteReset();						// Reset all the solid site velocities to zero

	// Multi-grid operations
	DEPRECATED void LBM_explode(int RegionNumber);		// Explode populations from coarse to fine
	DEPRECATED void LBM_coalesce(int RegionNumber);		// Coalesce populations from fine to coarse
	void LBM_addSubGrid(int RegionNumber);				// Add and initialise subgrid structure for a given region number

	// IO methods
	void io_textout(std::string output_tag);	// Writes out the contents of the class as well as any subgrids to a text file
	void io_fgaout();							// Wrapper for _io_fgaout with 2/3D checking 
	void io_restart(eIOFlag IO_flag);			// Reads/writes data from/to the global restart file
	void io_probeOutput();						// Output routine for point probes
	void io_lite(double tval, std::string Tag);	// Generic writer to individual files with Tag
	int io_hdf5(double tval);					// HDF5 writer returning integer to indicate success or failure

private :
	void _io_fgaout(int timeStepL0);		// Writes out the macroscopic velocity components for the class as well as any subgrids 
											// to a different .fga file for each subgrid. .fga format is the one used for Unreal 
											// Engine 4 VectorField object.
	// Private optimised functions
	void _LBM_stream_opt(int i, int j, int k, int id, eType type_local, int subcycle);
	void _LBM_coalesce_opt(int i, int j, int k, int id, int v);
	void _LBM_explode_opt(int id, int v, int src_x, int src_y, int src_z);
	void _LBM_collide_opt(int id, double omega_s);
	void _LBM_macro_opt(int i, int j, int k, int id, eType type_local);
	void _LBM_forceGrid_opt(int id);
	double _LBM_equilibrium_opt(int id, int v);

public :
	void LBM_multi_opt(int subcycle = 0);


};

#endif
