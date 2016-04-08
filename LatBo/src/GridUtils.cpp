#include "../inc/stdafx.h"
#include "../inc/GridObj.h"
#include <sstream>		// String stream package
#include <iostream>
#include "../inc/definitions.h"
#include "../inc/globalvars.h"
#include "../inc/MpiManager.h"

// Default Constructor/Destructor
GridUtils::GridUtils(void)
{
}

GridUtils::~GridUtils(void)
{
}

// ***************************************************************************************************
// Returns a vector which is the cross product of a and b
std::vector<double> GridUtils::crossprod(std::vector<double> a, std::vector<double> b) {

	// Declare resulting vector
	std::vector<double> result;

	result.push_back( a[1]*b[2] - a[2]*b[1] );
	result.push_back( a[2]*b[0] - a[0]*b[2] );
	result.push_back( a[0]*b[1] - a[1]*b[0] );

	return result;

}

// ***************************************************************************************************
// Returns a vector which is a minus b
std::vector<double> GridUtils::subtract(std::vector<double> a, std::vector<double> b) {

	std::vector<double> result;
	for (size_t i = 0; i < a.size(); i++) {
		result.push_back( a[i] - b[i] );
	}
	return result;

}

// ***************************************************************************************************
// Returns a vector which is sum of a and b
std::vector<double> GridUtils::add(std::vector<double> a, std::vector<double> b) {

	std::vector<double> result;
	for (size_t i = 0; i < a.size(); i++) {
		result.push_back( a[i] + b[i] );
	}
	return result;

}

// ***************************************************************************************************
// Returns a vector which is scalar multiplied by vector
std::vector<double> GridUtils::vecmultiply(double scalar, std::vector<double> vec) {

	std::vector<double> result;
	for (size_t i = 0; i < vec.size(); i++) {
		result.push_back( vec[i] * scalar );
	}
	return result;
}
	

// ***************************************************************************************************
// Returns a vector with n uniformly spaced values between min and max
std::vector<double> GridUtils::linspace(double min, double max, int n)
{
	// Declare resulting vector
	std::vector<double> result;

	// Set counter to zero
	int count = 0;

	// Number of values
	int numvals = n - 1; // Cast n to a double to use floor

	// Loop
	for (int i = 0; i <= n-2; i++)
	{
		double temp = min + i * (max - min) / numvals;
		result.insert(result.begin() + count, temp); // Insert element
		count += 1;
	}

	// Add last element
	result.insert(result.begin() + count, max);

	// Return vector
	return result;
}

// ***************************************************************************************************

// Like linspace but spaces elements by 1
std::vector<int> GridUtils::onespace(int min, int max)
{
	// Declare resulting array
	std::vector<int> result;

	// Loop and insert elements
	for (int i = 0; i <= (max-min); i++)
	{
		result.insert( result.begin() + i, min + i ); // Insert element

	}

	// Return vector
	return result;
}

// ***************************************************************************************************

// Functions to compute the magnitude of a vector. Overloads allow different input types
// 2D vector with arguments supplied separately
double GridUtils::vecnorm( double val1, double val2 )
{
	double result;

	result = sqrt( pow(val1,2) + pow(val2,2) );

	return result;
}

// 3D vector with arguments supplied separately
double GridUtils::vecnorm( double val1, double val2, double val3 )
{
	double result;

	result = sqrt( pow(val1,2) + pow(val2,2) + pow(val3,2) );

	return result;
}

// Supplied as a vector
double GridUtils::vecnorm( double vec[] )
{
	double result;

#if (dims == 3)

		result = sqrt( pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2) );

#else

		result = sqrt( pow(vec[0],2) + pow(vec[1],2) );

#endif

	return result;
}

// Supplied as a std::vector
double GridUtils::vecnorm( std::vector<double> vec )
{
	double result = 0.0;

	for (size_t d = 0; d < vec.size(); d++) {

		result += pow(vec[d],2);

	}

	return sqrt(result);
}

// ***************************************************************************************************

// Routine to map the global index of a coarse grid site to a corresponding fine site on the level below
std::vector<int> GridUtils::getFineIndices(int coarse_i, int x_start, int coarse_j, int y_start, int coarse_k, int z_start) {

	// Initialise result
	std::vector<int> fine_ind;

	// Map indices
	fine_ind.insert(fine_ind.begin(), 2*(coarse_i - x_start + 1) - 2 );
	fine_ind.insert(fine_ind.begin() + 1, 2*(coarse_j - y_start + 1) - 2 );
#if (dims == 3)
	fine_ind.insert(fine_ind.begin() + 2, 2*(coarse_k - z_start + 1) - 2 );
#else
	fine_ind.insert(fine_ind.begin() + 2, 0 );
#endif

	return fine_ind;
}
// ***************************************************************************************************

// Routine to map the global index of a fine grid site to parent coarse grid site on the level above
std::vector<int> GridUtils::getCoarseIndices(int fine_i, int x_start, int fine_j, int y_start, int fine_k, int z_start) {

	// Initialise result
	std::vector<int> coarse_ind;

	// Convert to top corner index if necessary
    if ((fine_i % 2) != 0) {
        fine_i = fine_i - 1;
    }
    if ((fine_j % 2) != 0) {
        fine_j = fine_j - 1;
    }
	if ((fine_k % 2) != 0) {
        fine_k = fine_k - 1;
    }

	// Reverse map indices
	coarse_ind.insert( coarse_ind.begin(), (fine_i / 2) + x_start );
	coarse_ind.insert( coarse_ind.begin() + 1, (fine_j / 2) + y_start );
#if (dims == 3)
	coarse_ind.insert( coarse_ind.begin() + 2, (fine_k / 2) + z_start );
#else
	coarse_ind.insert( coarse_ind.begin() + 2, 0 );
#endif

	return coarse_ind;
}
// ***************************************************************************************************

// Dot Product
double GridUtils::dotprod(std::vector<double> vec1, std::vector<double> vec2) {

	// Declare scalar answer
    double answer = 0.0;

	// Do dot product
    for (size_t i = 0; i < vec1.size(); i++) {
        answer += vec1[i] * vec2[i];
    }

	// Return answer
    return answer;
}

// ***************************************************************************************************

// Multiplies matrix A by vector x.
std::vector<double> GridUtils::matrix_multiply(const std::vector< std::vector<double> >& A, const std::vector<double>& x) {

	// Check to makes sure dimensions are correct
	if (A[0].size() != x.size()) {
		std::cout << "Error: See Log File" << std::endl;
		*logfile << "Dimension mismatch -- cannot proceed. Exiting." << std::endl;
		int ignore = system("pause");
		exit(LATBO_FAILED);
	}

	// Initialise answer
	std::vector<double> product (x.size(), 0.0);

	// Do multiplication
    for (size_t row = 0; row < A.size(); row++) {
        for (size_t col = 0; col < x.size(); col++) {
            // Multiply the row of A by the column of B to get the row, column of product.
			product[row] += A[row][col] * x[col];
		}
	}

	return product;
}

// ***************************************************************************************************
// Change index to position
double GridUtils::indexToPosition(int index, double dx) {

	return dx * ( static_cast<double>(index) + 0.5 );

}

// ***************************************************************************************************
// Routine to compute the opposite direction of the one supplied based on D2Q9 or D3Q19 numbering
size_t GridUtils::getOpposite(size_t direction) {

	size_t direction_opposite;

	// If rest particle then opposite is simply itself
	if (direction == nVels-1) {

		direction_opposite = direction;

	} else {

		/*	If direction is even, then opposite is direction+1.
			If direction is odd, then opposite is direction-1.
			e.g. direction 0 (+x direction) has opposite 1 (-x direction) --> +1
			however, direction 1 has opposite 0 --> -1
			Hence we can add (-1 ^ direction) so it alternates between +/-1
		*/

		direction_opposite = direction + (int)pow(-1,direction);

	}

	return direction_opposite;

}

// ***************************************************************************************************
// Function to find whether the recv layer containing local site i,j,k links to an adjacent or periodic neighbour rank.
// Takes in the site indices (local) and the lattice direction in which to check.
bool GridUtils::isOverlapPeriodic(int i, int j, int k, const GridObj& pGrid) {

	// Local declarations
	int exp_MPI_coords[dims], act_MPI_coords[dims], MPI_dims[dims];
	int shift[3] = {0, 0, 0};

	// Initialise local variables
	MPI_dims[0] = Xcores;
	MPI_dims[1] = Ycores;
#if (dims == 3)
	MPI_dims[2] = Zcores;
#endif

	// Define shifts based on which overlap we are on

	// X
	if (GridUtils::isOnRecvLayer(pGrid.XPos[i],0,1)) {
		shift[0] = 1;
	} else if (GridUtils::isOnRecvLayer(pGrid.XPos[i],0,0)) {
		shift[0] = -1;
	}

	// Y
	if (GridUtils::isOnRecvLayer(pGrid.YPos[j],1,1)) {
		shift[1] = 1;
	} else if (GridUtils::isOnRecvLayer(pGrid.YPos[j],1,0)) {
		shift[1] = -1;
	}

#if (dims == 3)
	// Z
	if (GridUtils::isOnRecvLayer(pGrid.ZPos[k],2,1)) {
		shift[2] = 1;
	} else if (GridUtils::isOnRecvLayer(pGrid.ZPos[k],2,0)) {
		shift[2] = -1;
	}
#endif

	// Loop over each Cartesian direction
	for (int d = 0; d < dims; d++) {
		// Define expected (non-periodic) MPI coordinates of neighbour rank
		exp_MPI_coords[d] = MpiManager::MPI_coords[d] + shift[d];

		// Define actual MPI coordinates of neighbour rank (accounting for periodicity)
		act_MPI_coords[d] = (MpiManager::MPI_coords[d] + shift[d] + MPI_dims[d]) % MPI_dims[d];

		// If there is a difference then rank is periodically linked to its neighbour and the overlap
		// site is from a periodic rank so return early
		if (exp_MPI_coords[d] != act_MPI_coords[d]) return true;
	}
	

	// Expected and actual are the same so the neighbour rank is not periodically linked and the overlap site
	// is from an adjacent neighbour not a periodic one.
	return false;

}

// ***************************************************************************************************
// Function to find whether a site with global indices provided is on a given grid or not
bool GridUtils::isOnThisRank(int gi, int gj, int gk, const GridObj& pGrid) {
	
	auto found_x = std::find(pGrid.XInd.begin(), pGrid.XInd.end(), gi);
	auto found_y = std::find(pGrid.YInd.begin(), pGrid.YInd.end(), gj);
	auto found_z = std::find(pGrid.ZInd.begin(), pGrid.ZInd.end(), gk);

	if (	found_x != pGrid.XInd.end() && found_y != pGrid.YInd.end()


#if (dims == 3)
		&& found_z != pGrid.ZInd.end()
#endif
		) {

		return true;

	} else {

		return false;
	}

}

// ***************************************************************************************************
// Overloaded function to find whether a global index gl == (i,j, or k) is on a given grid or not
bool GridUtils::isOnThisRank(int gl, int xyz, const GridObj& pGrid) {

	switch (xyz) {

	case 0:
		{
		// X direction
		auto found_x = std::find(pGrid.XInd.begin(), pGrid.XInd.end(), gl);

		if (found_x != pGrid.XInd.end()) return true;		
		else return false;
		}

	case 1:
		{
		// Y direction
		auto found_y = std::find(pGrid.YInd.begin(), pGrid.YInd.end(), gl);

		if (found_y != pGrid.YInd.end()) return true;		
		else return false;
		}

	case 2:
		{
		// Z direction
		auto found_z = std::find(pGrid.ZInd.begin(), pGrid.ZInd.end(), gl);

		if (found_z != pGrid.ZInd.end()) return true;		
		else return false;
		}

	}

	return false;

}

// ***************************************************************************************************
// Routine to see whether the specified refined region intersects with the span of the provided parent grid
bool GridUtils::hasThisSubGrid(const GridObj& pGrid, int RegNum) {

	// Loop over over X range of subgrid and check for matching index on parent grid
	for (size_t i = RefXstart[pGrid.level][RegNum]; i <= RefXend[pGrid.level][RegNum]; i++) {
		auto found_i = std::find(pGrid.XInd.begin(), pGrid.XInd.end(), i);
		if (found_i != pGrid.XInd.end()) break;		// If a match is found then chance that range intersects parent grid indices
		else if (i == RefXend[pGrid.level][RegNum] && found_i == pGrid.XInd.end()) return false;	// Got to the end and X is not intersecting
	}

	// Loop over over Y range
	for (size_t j = RefYstart[pGrid.level][RegNum]; j <= RefYend[pGrid.level][RegNum]; j++) {
		auto found_j = std::find(pGrid.YInd.begin(), pGrid.YInd.end(), j);
		if (found_j != pGrid.YInd.end()) break;
		else if (j == RefYend[pGrid.level][RegNum] && found_j == pGrid.YInd.end()) return false;
	}

#if (dims == 3)
	// Loop over over Z range
	for (size_t k = RefZstart[pGrid.level][RegNum]; k <= RefZend[pGrid.level][RegNum]; k++) {
		auto found_k = std::find(pGrid.ZInd.begin(), pGrid.ZInd.end(), k);
		if (found_k != pGrid.ZInd.end()) break;
		else if (k == RefZend[pGrid.level][RegNum] && found_k == pGrid.ZInd.end()) return false;
	}
#endif

	return true;
	
}

// ***************************************************************************************************
// Pass null pointer (ptr*) by reference and update it when matching grid is found in hierarchy pointed 
// to by Grids* which we also pass by reference to save copying it
void GridUtils::getGrid(GridObj*& Grids, int level, int region, GridObj*& ptr) {

	// Check supplied grid for a match
	if (Grids->level == level && Grids->region_number == region) {
		ptr = Grids;
		return;

	} else {

		// Loop through array of subgrids on this grid getting each one by reference
		for (GridObj& g : Grids->subGrid) {

			// Create pointer to subgrid
			GridObj* G = &g;

			// Region must match
			if (G->region_number == region) {

				// Look for match on this subgrid
				GridUtils::getGrid(G, level, region, ptr);

				// If a match found then return the non-null pointer
				if (ptr != NULL) {
					return;
				}

			}

		}

	}
	
	// Specified grid has not been found
	return;

}
// ***************************************************************************************************
// Wrapper routine to check whether a site is on the inner MPI overlap based on its position by checking 
// all possible sender layer locations
bool GridUtils::isOnSenderLayer(double pos_x, double pos_y, double pos_z) {

	/* Checks whether X,Y or Z are on a sender layer first.
	 * If true, ensures that neither of the other two are on a receiver layer.
	 * If true then site must be on a sender layer. */

	if (
	(
		// X on sender
		(GridUtils::isOnSenderLayer(pos_x,0,0) || GridUtils::isOnSenderLayer(pos_x,0,1)) &&

		// Y and Z not recv
		(!GridUtils::isOnRecvLayer(pos_y,1,0) && !GridUtils::isOnRecvLayer(pos_y,1,1)
#if (dims == 3)
		&& !GridUtils::isOnRecvLayer(pos_z,2,0) && !GridUtils::isOnRecvLayer(pos_z,2,1)
#endif
		)
	) || (
		// Y on sender
		(GridUtils::isOnSenderLayer(pos_y,1,0) || GridUtils::isOnSenderLayer(pos_y,1,1)) &&

		// X and Z not recv
		(!GridUtils::isOnRecvLayer(pos_x,0,0) && !GridUtils::isOnRecvLayer(pos_x,0,1)
#if (dims == 3)
		&& !GridUtils::isOnRecvLayer(pos_z,2,0) && !GridUtils::isOnRecvLayer(pos_z,2,1)

		)
	) || (
		// Z on sender
		(GridUtils::isOnSenderLayer(pos_z,2,0) || GridUtils::isOnSenderLayer(pos_z,2,1)) &&

		// X and Y not recv
		(!GridUtils::isOnRecvLayer(pos_x,0,0) && !GridUtils::isOnRecvLayer(pos_x,0,1)	&& 
		!GridUtils::isOnRecvLayer(pos_y,1,0) && !GridUtils::isOnRecvLayer(pos_y,1,1)
#endif
		)
	)
	) {
		return true;
	}

	return false;

}

// ***************************************************************************************************
// Routine to check whether a site is in the inner MPI overlap of the coarsest grid based on its position,
// whether it is an "x" = 0, "y" = 1 or "z" = 2 coordinate and which edge of the rank you are checking 
// either "min" = 0 or "max" = 1.
bool GridUtils::isOnSenderLayer(double site_position, int dir, int maxmin) {

	// Get instance of MPI manager
	MpiManager* mpim = MpiManager::getInstance();
	
	// Do checks on rank inner overlap regions dictated by the coarsest rank
	if (dir == 0) {

		if (maxmin == 1) {
			if (site_position > mpim->sender_layer_pos.X[2] && site_position < mpim->sender_layer_pos.X[3] ) return true;

		} else if (maxmin == 0) {
			if (site_position > mpim->sender_layer_pos.X[0] && site_position < mpim->sender_layer_pos.X[1] ) return true;

		}

	} else if (dir == 1) {

		if (maxmin == 1) {
			if (site_position > mpim->sender_layer_pos.Y[2] && site_position < mpim->sender_layer_pos.Y[3] ) return true;

		} else if (maxmin == 0) {
			if (site_position > mpim->sender_layer_pos.Y[0] && site_position < mpim->sender_layer_pos.Y[1] ) return true;

		}

	} else if (dir == 2) {

		if (maxmin == 1) {
			if (site_position > mpim->sender_layer_pos.Z[2] && site_position < mpim->sender_layer_pos.Z[3] ) return true;

		} else if (maxmin == 0) {
			if (site_position > mpim->sender_layer_pos.Z[0] && site_position < mpim->sender_layer_pos.Z[1] ) return true;

		}

	} else {

		// Invalid string indicating direction
		std::cout << "Error: See Log File" << std::endl;
		*logfile << "Invalid direction specified in GridUtils::isOnSenderLayer() / GridUtils::isOnRecvLayer(). Exiting." << std::endl;
		exit(LATBO_FAILED);

	}


	return false;

}

// ***************************************************************************************************
// Wrapper routine to check whether a site is on the outer MPI overlap based on its position by checking 
// all possible recv layer locations
bool GridUtils::isOnRecvLayer(double pos_x, double pos_y, double pos_z) {

	/* If any of the coordinates is in the receiver layer then the site will always 
	 * be a receiver layer site regardless of the other coordinates so the logic is
	 * simple. */

	if (	GridUtils::isOnRecvLayer(pos_x,0,0) || GridUtils::isOnRecvLayer(pos_x,0,1) ||
			GridUtils::isOnRecvLayer(pos_y,1,0) || GridUtils::isOnRecvLayer(pos_y,1,1)
#if (dims == 3)
			|| GridUtils::isOnRecvLayer(pos_z,2,0) || GridUtils::isOnRecvLayer(pos_z,2,1)
#endif
		) {
			return true;
	}

	return false;

}

// ***************************************************************************************************
// Routine to check whether a site is in the outer MPI overlap of the coarsest grid based on its position,
// whether it is an "x" = 0, "y" = 1 or "z" = 2 coordinate and which edge of the rank you are checking 
// either "min" = 0 or "max" = 1.
bool GridUtils::isOnRecvLayer(double site_position, int dir, int maxmin) {

	// Get instance of MPI manager
	MpiManager* mpim = MpiManager::getInstance();
	
	// Do checks on rank outer overlap regions dictated by the coarsest rank
	if (dir == 0) {

		if (maxmin == 1) {
			if (site_position > mpim->recv_layer_pos.X[2] && site_position < mpim->recv_layer_pos.X[3] ) return true;

		} else if (maxmin == 0) {
			if (site_position > mpim->recv_layer_pos.X[0] && site_position < mpim->recv_layer_pos.X[1] ) return true;

		}

	} else if (dir == 1) {

		if (maxmin == 1) {
			if (site_position > mpim->recv_layer_pos.Y[2] && site_position < mpim->recv_layer_pos.Y[3] ) return true;

		} else if (maxmin == 0) {
			if (site_position > mpim->recv_layer_pos.Y[0] && site_position < mpim->recv_layer_pos.Y[1] ) return true;

		}
		
	} else if (dir == 2) {

		if (maxmin == 1) {
			if (site_position > mpim->recv_layer_pos.Z[2] && site_position < mpim->recv_layer_pos.Z[3] ) return true;

		} else if (maxmin == 0) {
			if (site_position > mpim->recv_layer_pos.Z[0] && site_position < mpim->recv_layer_pos.Z[1] ) return true;

		}

	} else {

		// Invalid string indicating direction
		std::cout << "Error: See Log File" << std::endl;
		*logfile << "Invalid direction specified in GridUtils::isOnSenderLayer() / GridUtils::isOnRecvLayer(). Exiting." << std::endl;
		exit(LATBO_FAILED);

	}

	return false;

}

// ***************************************************************************************************
// Creates output directory with filename and path passed as a string and returns int specifying whether
// creation succeeded or failed.
int GridUtils::createOutputDirectory(std::string path_str) {

	int result = 9; // Return code of directory creation

	// Create output directory if it does not already exist
	std::string command = "mkdir -p " + path_str;

	// Only get rank 0 to create output directory
#ifdef BUILD_FOR_MPI

	#ifdef _WIN32   // Running on Windows
		if (MpiManager::my_rank == 0)
			result = CreateDirectoryA((LPCSTR)path_str.c_str(), NULL);
	#else   // Running on Unix system
		if (MpiManager::my_rank == 0)
			result = system(command.c_str());
	#endif // _WIN32

#else // BUILD_FOR_MPI

	#ifdef _WIN32   // Running on Windows
		result = CreateDirectoryA((LPCSTR)path_str.c_str(), NULL);
	#else   // Running on Unix system
		result = system(command.c_str());
	#endif // _WIN32

#endif // BUILD_FOR_MPI

	return result;
}