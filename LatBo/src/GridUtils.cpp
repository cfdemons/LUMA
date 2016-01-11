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
// Function to find whether a site is on the overlap region of an MPI grid
bool GridUtils::isOnOverlap(unsigned int i, unsigned int j, unsigned int k, unsigned int N_lim, unsigned int M_lim, unsigned int K_lim) {

	// Source in overlap
	if	(
			// Case 1: X source in overlap, Y,Z take any value
			(
			(i == N_lim - 1 || i == 0) && ( (j < M_lim && j >= 0)
#if (dims == 3)
		&&	(k < K_lim && k >= 0)
#endif
		)	)

			// Case 2: Y source in overlap, X,Z take any value
		||	(
			(j == M_lim - 1 || j == 0) && ( (i < N_lim && i >= 0)
#if (dims == 3)
		&&	(k < K_lim && k >= 0)

		)	)

			// Case 3: Z source in overlap, Y,X take any value
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
// Function to find whether the specified site is on overlap from an adjacent or periodic neighbour rank.
// Takes in the site indices (MPI local) and the lattice direction in which to check.
bool GridUtils::isOverlapPeriodic(unsigned int i, unsigned int j, unsigned int k,
								  unsigned int N_lim, unsigned int M_lim, unsigned int K_lim,
								  unsigned int lattice_dir) {

	// Local declarations
	int exp_MPI_coords[dims], act_MPI_coords[dims], MPI_dims[dims], Lims[dims], Ind[dims];

	// Initialise local variables
	MPI_dims[0] = Xcores; Lims[0] = N_lim; Ind[0] = i;
	MPI_dims[1] = Ycores; Lims[1] = M_lim; Ind[1] = j;
#if (dims == 3)
	MPI_dims[2] = Zcores; Lims[2] = K_lim; Ind[2] = k;
#endif

	// Loop over each Cartesian direction
	for (int d = 0; d < dims; d++) {
		// If on X/Y/Z-overlap
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
// Function to find whether a site wiht global indices provided is on a given grid or not
bool GridUtils::isOnThisRank(unsigned int i, unsigned int j, unsigned int k, GridObj& pGrid) {

	if (
		// Different conditions when using MPI
#ifdef BUILD_FOR_MPI
		((int)i <= pGrid.XInd[pGrid.XInd.size()-2] && (int)i >= pGrid.XInd[1] ) &&
		((int)j <= pGrid.YInd[pGrid.YInd.size()-2] && (int)j >= pGrid.YInd[1] )
#if (dims == 3)
		&& ((int)k <= pGrid.ZInd[pGrid.ZInd.size()-2] && (int)k >= pGrid.ZInd[1] )
#endif

#else

		((int)i <= pGrid.XInd[pGrid.XInd.size()-1] && (int)i >= pGrid.XInd[1] ) &&
		((int)j <= pGrid.YInd[pGrid.YInd.size()-1] && (int)j >= pGrid.YInd[0] )
#if (dims == 3)
		&& ((int)k <= pGrid.ZInd[pGrid.ZInd.size()-1] && (int)k >= pGrid.ZInd[0] )
#endif


#endif

		) {

			return true;

	} else {

		return false;
	}


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

