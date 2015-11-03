/* Class defintion for the GridObj class
*/

#pragma once

#include <vector>
#include "ivector.h"
#include "IB_body.h"
#include "GridUtils.h"
#include <iostream>
#include <fstream>

// Base class
class GridObj
{

	// Allow MPI_manager objects to access Grid information
	friend class MPI_manager;

public:

	// Constructors & Destructor //
	GridObj( ); // Default constructor
	GridObj(int level, std::ofstream* logfile); // Basic grid constructor
	GridObj(int level, int RegionNumber, int rank, int max_ranks, std::ofstream* logfile); // MPI sub grid constructor with level, region and rank
	// MPI L0 constructor with level, rank, local size and global edges
	GridObj(int level, int rank, int max_ranks, std::vector<unsigned int> local_size, 
		std::vector< std::vector<unsigned int> > GlobalLimsInd, 
		std::vector< std::vector<double> > GlobalLimsPos,
		int my_coords[],
		std::ofstream* logfile);
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
private :
	std::vector<double> XPos; // Vectors of positions of sites
	std::vector<double> YPos;
	std::vector<double> ZPos;

	// Inlet velocity profile
	std::vector<double> ux_in, uy_in, uz_in;

	// Vector nodal properties
	// Flattened 4D arrays (i,j,k,vel)
	ivector<double> f;
	ivector<double> feq;
	ivector<double> u;
	ivector<double> force_xyz;
	ivector<double> force_i;

	// Scalar nodal properties
	// Flattened 3D arrays (i,j,k)
	ivector<double> rho;
	ivector<int> LatTyp;

	// Grid scalars
	double dx, dy, dz;	// Physical spacing
	int region_number;	// ID of region at a particular level in the embedded grid hierarchy

	// Time averaged statistics
	ivector<double> rho_timeav;		// Time-averaged density at each grid point (i,j,k)
	ivector<double> ui_timeav;		// Time-averaged velocity at each grid point (i,j,k,dims)
	ivector<double> uiuj_timeav;	// Time-averaged velocity products at each grid point (i,j,k,2*dims)

	// IBM objects
	std::vector<IB_body> iBody;		// Array of immersed boundary bodies


	// Public data members
public :

	int level;						// Level in embedded grid hierarchy
	double dt;						// Physical time step size
	unsigned int t;					// Number of completed iterations
	double nu;						// Kinematic viscosity (in lattice units)
	double omega;					// Relaxation frequency
	std::vector<double> mrt_omega;	// Relaxation frequencies in moment space (for MRT)
	int my_rank, max_ranks;			// MPI rank and totla number of ranks
	GridUtils gUtils;				// Utility class
	
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
	void LBM_init_grid(std::vector<unsigned int> local_size, 
		std::vector< std::vector<unsigned int> > GlobalLimsInd, 
		std::vector< std::vector<double> > GlobalLimsPos);		// Initialise top level grid with fields and labels
	void LBM_init_subgrid(double offsetX, double offsetY, double offsetZ, 
		double dx0, double omega_coarse, std::vector<double> mrt_omega_coarse);	// Initialise subgrid with all quantities
	void LBM_init_bound_lab();		// Initialise labels for objects and walls
	void LBM_init_refined_lab();	// Initialise labels for refined regions
	void LBM_init_getInletProfile();	// Initialise the store for inlet profile data from file

	// LBM operations
	void LBM_multi(bool IBM_flag);				// Launch the multi-grid kernel
	void LBM_collide(bool core_flag);			// Apply collision + 1 overload for equilibrium calculation
	double LBM_collide(int i, int j, int k, int v);
	void LBM_mrt_collide(ivector<double>& f_new, int i, int j, int k);	// MRT collision operation
	void LBM_stream();							// Stream populations
	void LBM_macro();							// Compute macroscopic quantities + 1 overload for single site
	void LBM_macro(int i, int j, int k);
	void LBM_boundary(int bc_type_flag);		// Apply boundary conditions
	void LBM_forcegrid(bool reset_flag);		// Apply a force to the grid points (or simply reset force vectors if flag is true)
	
	// Boundary operations
	void bc_applyZouHe(int label, int i, int j, int k, int M_lim, int K_lim);			// Application of Zou-He BC
	void bc_applyRegularised(int label, int i, int j, int k, int M_lim, int K_lim);		// Application of Regaulrised BC
	void bc_applyExtrapolation(int label, int i, int j, int k, int M_lim, int K_lim);	// Application of Extrapolation BC
	void bc_solid_site_reset();	// Reset all the solid site velocities to zero

	// Multi-grid operations
	void LBM_explode(int RegionNumber);			// Explode populations from coarse to fine
	void LBM_coalesce(int RegionNumber);		// Coalesce populations from fine to coarse
	void LBM_addSubGrid(int RegionNumber);		// Add and initialise subgrid structure for a given region number

	// IO methods
	void io_write_body_pos();					// Write out IB_body positions to text files
	void io_write_lift_drag();					// Write out IB_body lift and drag
	void io_textout(std::string output_tag);	// Writes out the contents of the class as well as any subgrids to a text file
	void io_restart(bool IO_flag);				// Reads/writes data from/to the global restart file
	void io_probe_output();						// Output routine for point probes
	void io_vtkwriter(double tval);				// VTK writer
	void io_tecplot(double tval);				// TecPlot write out

	// IBM methods
	void ibm_build_body(int body_type);						// Build a new pre-fab body
	void ibm_initialise();									// Initialise a built immersed body
	double ibm_deltakernel(double rad, double dilation);	// Evaluate kernel (delta function approximation)
	void ibm_interpol(unsigned int ib);						// Interpolation of velocity field onto markers of ib-th body
	void ibm_spread(unsigned int ib);						// Spreading of restoring force from ib-th body to grid
	void ibm_findsupport(unsigned int ib, unsigned int m);	// Populates support information for the m-th marker of ib-th body
	void ibm_computeforce(unsigned int ib);					// Compute restorative force at each marker in ib-th body
	double ibm_findepsilon(unsigned int ib);				// Method to find epsilon weighting parameter for ib-th body
	// Biconjugate gradient stablised method for solving asymmetric linear system required by finding epsilon
	double ibm_bicgstab(std::vector< std::vector<double> >& Amatrix, std::vector<double>& bVector, std::vector<double>& epsilon, 
						   double tolerance, unsigned int maxiterations);

	// Flexible methods
	void ibm_jacowire(unsigned int ib);				// Computes the tension and position of a 2D inextensible, flexible iBody filament
	void ibm_position_update( unsigned int ib );	// Updates the position of deformable body markers
	void ibm_position_update_grp( unsigned int group );	// Updates the positions of deformable bodies in a group using on the flexible group member
	// Methods to solve the Jacobian system
	void ibm_banbks(double **a, unsigned long n, unsigned int m1, unsigned int m2, double **al,
	unsigned long indx[], double b[]);
	void ibm_bandec(double **a, unsigned long n, unsigned int m1, unsigned int m2, double **al,
	unsigned long indx[], double *d);

	
};

