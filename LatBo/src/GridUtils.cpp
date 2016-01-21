#include "../inc/stdafx.h"
#include "../inc/GridUtils.h"
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
// Returns a vector with n uniformly spaced values between min and max
std::vector<double> GridUtils::linspace(double min, double max, int n)
{
	// Decalre resulting vector
	std::vector<double> result;

	// Set counter to zero
	int count = 0;

	// Loop
	for (int i = 0; i <= n-2; i++)
	{
		double temp = min + i*(max-min)/(floor( (double)n ) - 1); // Cast n to a double to use floor
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

// Supplied as a std::vector (as reference just in case vector is big)
double GridUtils::vecnorm( std::vector<double>& vec )
{
	double result = 0.0;

	for (size_t d = 0; d < vec.size(); d++) {

		result += pow(vec[d],2);

	}

	return sqrt(result);
}

// ***************************************************************************************************

// Routine to map the index of a coarse grid site to a corresponding fine site on the level below
std::vector<int> GridUtils::indmapref(int coarse_i, int x_start, int coarse_j, int y_start, int coarse_k, int z_start) {

	// Initialise result
	std::vector<int> fine_ind;

	// Map indices
	fine_ind.insert(fine_ind.begin(), 2*(coarse_i - x_start + 1) - 2 );
	fine_ind.insert(fine_ind.begin() + 1, 2*(coarse_j - y_start + 1) - 2 );
	fine_ind.insert(fine_ind.begin() + 2, 2*(coarse_k - z_start + 1) - 2 );

	return fine_ind;
}
// ***************************************************************************************************

// Routine to map the index of a coarse grid site to a corresponding fine site on the level below
std::vector<int> GridUtils::revindmapref(int fine_i, int x_start, int fine_j, int y_start, int fine_k, int z_start) {

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
	coarse_ind.insert( coarse_ind.begin() + 2, (fine_k / 2) + z_start );

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
std::vector<double> GridUtils::matrix_multiply(std::vector< std::vector<double> >& A, std::vector<double>& x) {

	// Check to makes sure dimensions are correct
	if (A[0].size() != x.size()) {
		std::cout << "Error: See Log File" << std::endl;
		*logfile << "Dimension mismatch -- cannot proceed. Exiting." << std::endl;
		int ignore = system("pause");
		exit(EXIT_FAILURE);
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
// Function to find whether a site with local indices i,j,k is on the edge of the supplied GridObj
// On L0 this is equivalent to being on the MPI overlap.
bool GridUtils::isOnEdge(unsigned int i, unsigned int j, unsigned int k, GridObj& pGrid) {

	// Get Grid size
	unsigned int N_lim = pGrid.XInd.size();
	unsigned int M_lim = pGrid.YInd.size();
	unsigned int K_lim = pGrid.ZInd.size();


	// Check each coordinate to see if it is on edge
	if	(
			// Case 1: X on edge, Y,Z take any value
			(
			(i == N_lim - 1 || i == 0) && ( (j < M_lim && j >= 0)
#if (dims == 3)
		&&	(k < K_lim && k >= 0)
#endif
		)	)

			// Case 2: Y on edge, X,Z take any value
		||	(
			(j == M_lim - 1 || j == 0) && ( (i < N_lim && i >= 0)
#if (dims == 3)
		&&	(k < K_lim && k >= 0)

		)	)

			// Case 3: Z on edge, Y,X take any value
		||	(
			(k == K_lim - 1 || k == 0) && (	(j < M_lim && j >= 0)
													&&	(i < N_lim && i >= 0)	)

			)
#else
		)	)
#endif

		) {

			return true;

	} else {

		return false;
	}


}

// ***************************************************************************************************
// Function to find whether the overlap belonging to local site i,j,k is from an adjacent or periodic neighbour rank.
// Takes in the site indices (local) and the lattice direction in which to check.
bool GridUtils::isOverlapPeriodic(unsigned int i, unsigned int j, unsigned int k, 
								  GridObj& pGrid, unsigned int lattice_dir) {

	// Local declarations
	int exp_MPI_coords[dims], act_MPI_coords[dims], MPI_dims[dims], Lims[dims], Ind[dims];

	// Initialise local variables
	MPI_dims[0] = Xcores; Lims[0] = pGrid.XInd.size(); Ind[0] = i;
	MPI_dims[1] = Ycores; Lims[1] = pGrid.YInd.size(); Ind[1] = j;
#if (dims == 3)
	MPI_dims[2] = Zcores; Lims[2] = pGrid.ZInd.size(); Ind[2] = k;
#endif

	// Loop over each Cartesian direction
	for (int d = 0; d < dims; d++) {
		// If on X/Y/Z-overlap (should be)
		if (Ind[d] == Lims[d] - 1 || Ind[d] == 0) {

			// Define expected MPI coordinates of neighbour rank
			exp_MPI_coords[d] = MpiManager::MPI_coords[d] + c[d][lattice_dir];

			// Define actual MPI coordinates of neighbour rank (accounting for periodicity)
			act_MPI_coords[d] = (MpiManager::MPI_coords[d] + c[d][lattice_dir] + MPI_dims[d]) % MPI_dims[d];

			// If there is a difference then rank is periodically linked to its neighbour and the overlap
			// site is from a periodic rank so return early
			if (exp_MPI_coords[d] != act_MPI_coords[d]) return true;
		}
	}

	// Expected and actual are the same so the neighbour rank is not periodically linked and the overlap site
	// is from an adjacent neighbour not a periodic one.
	return false;

}

// ***************************************************************************************************
// Function to find whether a site with global indices provided is on a given grid or not (doesn't include the overlap)
bool GridUtils::isOnThisRank(unsigned int gi, unsigned int gj, unsigned int gk, GridObj& pGrid) {

	if (
		// Different conditions when using MPI due to extra overlap cells
#ifdef BUILD_FOR_MPI
		((int)gi <= pGrid.XInd[pGrid.XInd.size()-2] + 1 && (int)gi >= pGrid.XInd[1] - 1 ) &&
		((int)gj <= pGrid.YInd[pGrid.YInd.size()-2] + 1 && (int)gj >= pGrid.YInd[1] - 1 )
#if (dims == 3)
		&& ((int)gk <= pGrid.ZInd[pGrid.ZInd.size()-2] + 1 && (int)gk >= pGrid.ZInd[1] - 1 )
#endif

#else

		((int)gi <= pGrid.XInd[pGrid.XInd.size()-1] && (int)gi >= pGrid.XInd[0] ) &&
		((int)gj <= pGrid.YInd[pGrid.YInd.size()-1] && (int)gj >= pGrid.YInd[0] )
#if (dims == 3)
		&& ((int)gk <= pGrid.ZInd[pGrid.ZInd.size()-1] && (int)gk >= pGrid.ZInd[0] )
#endif


#endif

		) {

			return true;

	} else {

		return false;
	}


}

// ***************************************************************************************************
// Overloaded function to find whether a global index gl == (i,j, or k) is on a given grid or not (doesn't include the overlap)
bool GridUtils::isOnThisRank(unsigned int gl, unsigned int xyz, GridObj& pGrid) {

	switch (xyz) {

	case 0:
		// X direction
#ifdef BUILD_FOR_MPI	// Allow for overlap but ignore periodicity
		if ((int)gl <= pGrid.XInd[pGrid.XInd.size()-2] + 1 && (int)gl >= pGrid.XInd[1] - 1)
#else
		if ((int)gl <= pGrid.XInd[pGrid.XInd.size()-1] && (int)gl >= pGrid.XInd[0])
#endif
		{
			return true;
		
		} else {

			return false;
		}

	case 1:
		// Y direction
#ifdef BUILD_FOR_MPI
		if ((int)gl <= pGrid.YInd[pGrid.YInd.size()-2] + 1 && (int)gl >= pGrid.YInd[1] - 1)
#else
		if ((int)gl <= pGrid.YInd[pGrid.YInd.size()-1] && (int)gl >= pGrid.YInd[0])
#endif
		{
			return true;
		
		} else {

			return false;
		}

	case 2:
		// Z direction
#ifdef BUILD_FOR_MPI
		if ((int)gl <= pGrid.ZInd[pGrid.ZInd.size()-2] + 1 && (int)gl >= pGrid.ZInd[1] - 1)
#else
		if ((int)gl <= pGrid.ZInd[pGrid.ZInd.size()-1] && (int)gl >= pGrid.ZInd[0])
#endif
		{
			return true;
		
		} else {

			return false;
		}

	}

}

// ***************************************************************************************************
// Routine to see whether the specified refined region intersects with the span of the provided parent grid
bool GridUtils::hasThisSubGrid(GridObj& pGrid, int RegNum) {


	// Loop through every global point on the given grid and if one 
	// of them exists within the refined region then return true.
	for (unsigned int i : pGrid.XInd) {
		for (unsigned int j : pGrid.YInd) {
			for (unsigned int k : pGrid.ZInd) {

				if	(
					(i >= RefXstart[pGrid.level][RegNum] && i <= RefXend[pGrid.level][RegNum]) &&
					(j >= RefYstart[pGrid.level][RegNum] && j <= RefYend[pGrid.level][RegNum])
#if (dims == 3)
					&& (k >= RefZstart[pGrid.level][RegNum] && k <= RefZend[pGrid.level][RegNum])
#endif
				) {				
					return true;
				}

			}
		}
	}

	return false;
	
}

// ***************************************************************************************************
// Pass pointer (ptr*) by reference and update it when matching grid is found in hierarchy pointed to by Grids*
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
// Routine to check whether a site is in the inner MPI overlap of the coarsest grid based on its position,
// whether it is an 'x','y' or 'z' coordinate and which edge of the rank you are checking either "min" or "max".
bool GridUtils::isOnSenderLayer(double site_position, char dir, char* maxmin) {

	// Get instance of MPI manager
	MpiManager* mpim = MpiManager::getInstance();
	
	// Do checks on rank inner overlap regions dictated by the coarsest rank
	switch (dir) {

	case 'x':

		if (!strcmp(maxmin,"max") ) {	// Note that strcmp() returns 0 if they are equal
			if (site_position > mpim->sender_layer_pos.X[2] && site_position < mpim->sender_layer_pos.X[3] ) return true;

		} else if (!strcmp(maxmin,"min")) {
			if (site_position > mpim->sender_layer_pos.X[0] && site_position < mpim->sender_layer_pos.X[1] ) return true;

		}
		break;

	case 'y':

		if (!strcmp(maxmin,"max") ) {
			if (site_position > mpim->sender_layer_pos.Y[2] && site_position < mpim->sender_layer_pos.Y[3] ) return true;

		} else if (!strcmp(maxmin,"min")) {
			if (site_position > mpim->sender_layer_pos.Y[0] && site_position < mpim->sender_layer_pos.Y[1] ) return true;

		}
		break;


	case 'z':

		if (!strcmp(maxmin,"max") ) {
			if (site_position > mpim->sender_layer_pos.Z[2] && site_position < mpim->sender_layer_pos.Z[3] ) return true;

		} else if (!strcmp(maxmin,"min")) {
			if (site_position > mpim->sender_layer_pos.Z[0] && site_position < mpim->sender_layer_pos.Z[1] ) return true;

		}
		break;

	}

	return false;

}

// ***************************************************************************************************
// Routine to check whether a site is in the outer MPI overlap of the coarsest grid based on its position,
// whether it is an 'x','y' or 'z' coordinate and which edge of the rank you are checking either "min" or "max".
bool GridUtils::isOnRecvLayer(double site_position, char dir, char* maxmin) {

	// Get instance of MPI manager
	MpiManager* mpim = MpiManager::getInstance();
	
	// Do checks on rank outer overlap regions dictated by the coarsest rank
	switch (dir) {

	case 'x':

		if (!strcmp(maxmin,"max") ) {
			if (site_position > mpim->recv_layer_pos.X[2] && site_position < mpim->recv_layer_pos.X[3] ) return true;

		} else if (!strcmp(maxmin,"min")) {
			if (site_position > mpim->recv_layer_pos.X[0] && site_position < mpim->recv_layer_pos.X[1] ) return true;

		}
		break;

	case 'y':

		if (!strcmp(maxmin,"max") ) {
			if (site_position > mpim->recv_layer_pos.Y[2] && site_position < mpim->recv_layer_pos.Y[3] ) return true;

		} else if (!strcmp(maxmin,"min")) {
			if (site_position > mpim->recv_layer_pos.Y[0] && site_position < mpim->recv_layer_pos.Y[1] ) return true;

		}
		break;


	case 'z':

		if (!strcmp(maxmin,"max") ) {
			if (site_position > mpim->recv_layer_pos.Z[2] && site_position < mpim->recv_layer_pos.Z[3] ) return true;

		} else if (!strcmp(maxmin,"min")) {
			if (site_position > mpim->recv_layer_pos.Z[0] && site_position < mpim->recv_layer_pos.Z[1] ) return true;

		}
		break;

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

