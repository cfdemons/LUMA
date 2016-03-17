#pragma once

/** GridObj class which represents a lattice */

#include <vector>
#include "IVector.h"
#include "IBBody.h"
#include <iostream>
#include <fstream>

// Base class
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

	// 1D subgrid array (size = NumReg)
	std::vector<GridObj> subGrid;

	// Start and end indices of corresponding coarse level
	// When using MPI these values are local to a particular coarse grid
	size_t CoarseLimsX[2];
	size_t CoarseLimsY[2];
	size_t CoarseLimsZ[2];

	// 1D arrays
public :
	std::vector<int> XInd; // Vectors of indices
	std::vector<int> YInd;
	std::vector<int> ZInd;
	std::vector<double> XPos; // Vectors of positions of sites
	std::vector<double> YPos;
	std::vector<double> ZPos;

private :
	// Inlet velocity profile
	std::vector<double> ux_in, uy_in, uz_in;

	// Vector nodal properties
	// Flattened 4D arrays (i,j,k,vel)
	IVector<double> f;
	IVector<double> feq;
	IVector<double> u;
	IVector<double> force_xyz;
	IVector<double> force_i;

	// Scalar nodal properties
	// Flattened 3D arrays (i,j,k)
	IVector<double> rho;

	// Grid scalars
	double dx, dy, dz;	// Physical spacing
	int region_number;	// ID of region at a particular level in the embedded grid hierarchy

	// Time averaged statistics
	IVector<double> rho_timeav;		// Time-averaged density at each grid point (i,j,k)
	IVector<double> ui_timeav;		// Time-averaged velocity at each grid point (i,j,k,dims)
	IVector<double> uiuj_timeav;	// Time-averaged velocity products at each grid point (i,j,k,2*dims)


	// Public data members
public :

	IVector<int> LatTyp;			// Flattened 3D array of site labels
	int level;						// Level in embedded grid hierarchy
	double dt;						// Physical time step size
	int t;					// Number of completed iterations
	double nu;						// Kinematic viscosity (in lattice units)
	double omega;					// Relaxation frequency
	std::vector<double> mrt_omega;	// Relaxation frequencies in moment space (for MRT)

	// Timing variables
	double timeav_mpi_overhead;		// Time of MPI communication
	double timeav_timestep;			// Time of a timestep


	/*
	***************************************************************************************************************
	********************************************* Member Methods **************************************************
	***************************************************************************************************************
	*/

public :

	// Initialisation functions
	void LBM_init_vel();		// Initialise the velocity field
	void LBM_init_rho();		// Initialise the density field
	void LBM_init_grid();		// Non-MPI wrapper for initialiser
	void LBM_init_grid(std::vector<int> local_size,
		std::vector< std::vector<int> > GlobalLimsInd,
		std::vector< std::vector<double> > GlobalLimsPos);		// Initialise top level grid with fields and labels
	void LBM_init_subgrid(GridObj& pGrid);	// Initialise subgrid with all quantities
	void LBM_init_bound_lab();				// Initialise labels for walls
	void LBM_init_solid_lab();				// Initialise labels for solid objects
	void LBM_init_refined_lab(GridObj& pGrid);	// Initialise labels for refined regions
	void LBM_init_getInletProfile();		// Initialise the store for inlet profile data from file

	// LBM operations
	void LBM_multi(bool IBM_flag);		// Launch the multi-grid kernel
	void LBM_collide();					// Apply collision + 1 overload for equilibrium calculation
	double LBM_collide(int i, int j, int k, int v);
	void LBM_mrt_collide(IVector<double>& f_new, int i, int j, int k);	// MRT collision operation
	void LBM_stream();							// Stream populations
	void LBM_macro();							// Compute macroscopic quantities + 1 overload for single site
	void LBM_macro(int i, int j, int k);
	void LBM_boundary(int bc_type_flag);		// Apply boundary conditions
	void LBM_forcegrid(bool reset_flag);		// Apply a force to the grid points (or simply reset force vectors if flag is true)

	// Boundary operations
	void bc_applyZouHe(int label, int i, int j, int k, int M_lim, int K_lim);			// Application of Zou-He BC
	void bc_applyRegularised(int label, int i, int j, int k, int M_lim, int K_lim);		// Application of Regaulrised BC
	void bc_applyExtrapolation(int label, int i, int j, int k, int M_lim, int K_lim);	// Application of Extrapolation BC
	void bc_applyBfl(int i, int j, int k);												// Application of BFL BC
	void bc_solid_site_reset();	// Reset all the solid site velocities to zero

	// Multi-grid operations
	void LBM_explode(int RegionNumber);			// Explode populations from coarse to fine
	void LBM_coalesce(int RegionNumber);		// Coalesce populations from fine to coarse
	void LBM_addSubGrid(int RegionNumber);		// Add and initialise subgrid structure for a given region number

	// IO methods
	void io_textout(std::string output_tag);	// Writes out the contents of the class as well as any subgrids to a text file
	void io_restart(bool IO_flag);				// Reads/writes data from/to the global restart file
	void io_probe_output();						// Output routine for point probes
	void io_vtkwriter(double tval);				// VTK writer
	void io_tecplot(double tval);				// TecPlot write out
	void io_tecplot_debug(double tval, std::string tag);	// Special debugging writer to help debug problems with MPI


};

