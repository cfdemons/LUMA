#pragma once

#include <vector>

// Base class
class GridObj
{

public:

	GridObj( ); // Default constructor
	GridObj(int level); // Constructor with level
	GridObj(int level, int RegionNumber); // Sub grid constructor
	~GridObj( ); // Default destructor


	/*	
	***************************************************************************************************************
	********************************************* Member Data *****************************************************
	***************************************************************************************************************
	*/
	
private :

	// Subgrid array (size = NumReg)
	std::vector<GridObj> subGrid;
	
	// Start and end indices of corresponding coarse level (unsigned integers)
	size_t CoarseLimsX[2];
	size_t CoarseLimsY[2];
	size_t CoarseLimsZ[2];

	// 1D arrays
	std::vector<int> XInd; // Vectors of indices
	std::vector<int> YInd;
	std::vector<int> ZInd;
	std::vector<double> XPos; // Vectors of positions of sites
	std::vector<double> YPos;
	std::vector<double> ZPos;

	// Vector nodal properties
	// Flattened 4D arrays (i,j,k,vel)
	std::vector<double> f;
	std::vector<double> feq;
	std::vector<double> u;

	// Scalar nodal properties
	// Flattened 3D arrays (i,j,k)
	std::vector<double> rho;
	std::vector<int> LatTyp;

	// Grid scalars
	double omega;
	double dx;
	double dy;
	double dz;
	int level;
	int region_number;

	// Time step
public :
	double dt;

	/*	
	***************************************************************************************************************
	********************************************* Member Methods **************************************************
	***************************************************************************************************************
	*/

public :

	// Initialisation functions
	void LBM_init_vel();		// Initialise the velocity field
	void LBM_init_rho();		// Initialise the density field
	void LBM_init_grid();		// Initialise top level grid with a velocity and denstiy field
	void LBM_init_subgrid(double offsetX, double offsetY, double offsetZ, double dx0, double omega_coarse);	// Initialise subgrid with all quantities

	// LBM operations
	void LBM_multi();							// Launch the multi-grid kernel
	void LBM_collide(bool core_flag);			// Apply collision + 1 overload
	double LBM_collide(int i, int j, int k, int v);
	void LBM_stream();							// Stream populations
	void LBM_macro();							// Compute macroscopic quantities
	void LBM_boundary(int bc_type_flag);		// Apply boundary conditions

	// Multi-grid operations
	void LBM_explode(int RegionNumber);					// Explode populations from coarse to fine
	void LBM_coalesce(int RegionNumber);				// Coalesce populations from fine to coarse

	// Potentially deprecated methods -- remove in future release
	bool isEdge(size_t i, size_t j, size_t k, int RegionNumber);	// Check whether point is edge of subGrid[RegionNum]
	bool isWithin(size_t i, size_t j, size_t k, int RegionNumber);	// Check whether point lies within particular subGrid[RegionNum]

	// Add subgrid
	void LBM_addSubGrid(int RegionNumber);		// Add and initialise subgrid structure for a given region number

	// IO methods
	void lbm_write3(int t);		// Writes out the contents of the class as well as any subgrids
	// EnsightGold methods
	void genCase(int nsteps, int saveEvery);	// Generate case file
	void genGeo();								// Generate geometry file
	void genVec(int fileNum);					// Generate vectors file
	void genScal(int fileNum);					// Generate scalars file

	
};

