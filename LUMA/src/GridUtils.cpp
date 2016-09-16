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

#include "../inc/stdafx.h"
#include "../inc/GridUtils.h"
#include "../inc/MpiManager.h"
#include "../inc/GridObj.h"

// Mappings of directions for specular reflection: col == normal direction, row == velocity
// Constants so initialise outside the class but in scope.
#if (L_dims == 3)

// THIS IS WRONG AS WE HAVE CHANGED TO D3Q27 NOW!!!!
const int GridUtils::dir_reflect[L_dims * 2][L_nVels] = 
	{
		{1, 0, 2, 3, 4, 5, 9, 8, 7, 6, 10, 11, 12, 13, 16, 17, 14, 15, 18}, 
		{1, 0, 2, 3, 4, 5, 9, 8, 7, 6, 10, 11, 12, 13, 16, 17, 14, 15, 18},
		{0, 1, 3, 2, 4, 5, 8, 9, 6, 7, 13, 12, 11, 10, 14, 15, 16, 17, 18},
		{0, 1, 3, 2, 4, 5, 8, 9, 6, 7, 13, 12, 11, 10, 14, 15, 16, 17, 18},
		{0, 1, 2, 3, 5, 4, 6, 7, 8, 9, 12, 13, 10, 11, 17, 16, 15, 14, 18},
		{0, 1, 2, 3, 5, 4, 6, 7, 8, 9, 12, 13, 10, 11, 17, 16, 15, 14, 18}
	};
#else
const int GridUtils::dir_reflect[L_dims * 2][L_nVels] = 
	{
		{1, 0, 2, 3, 7, 6, 5, 4, 8}, 
		{1, 0, 2, 3, 4, 6, 5, 4, 8},
		{0, 1, 3, 2, 6, 7, 4, 5, 8},
		{0, 1, 3, 2, 6, 7, 4, 5, 8}
	};
#endif

/// Default constructor
GridUtils::GridUtils(void)
{
}

/// Default destructor
GridUtils::~GridUtils(void)
{
}

// *****************************************************************************
/// \brief	Computes vector product.
/// \param	a	a vector.
/// \param	b	a second vector.
/// \return	a vector which is the cross product of a and b.
std::vector<double> GridUtils::crossprod(std::vector<double> a, std::vector<double> b) {

	// Declare resulting vector
	std::vector<double> result;

	result.push_back( a[1]*b[2] - a[2]*b[1] );
	result.push_back( a[2]*b[0] - a[0]*b[2] );
	result.push_back( a[0]*b[1] - a[1]*b[0] );

	return result;

}

// *****************************************************************************
/// \brief	Subtracts two vectors.
/// \param	a	a vector.
/// \param	b	a second vector.
/// \return	a vector which is a - b.
std::vector<double> GridUtils::subtract(std::vector<double> a, std::vector<double> b) {

	std::vector<double> result;
	for (size_t i = 0; i < a.size(); i++) {
		result.push_back( a[i] - b[i] );
	}
	return result;

}

// *****************************************************************************
/// \brief	Adds two vectors.
/// \param	a	a vector.
/// \param	b	a second vector.
/// \return vector which is a + b.
std::vector<double> GridUtils::add(std::vector<double> a, std::vector<double> b) {

	std::vector<double> result;
	for (size_t i = 0; i < a.size(); i++) {
		result.push_back( a[i] + b[i] );
	}
	return result;

}

// *****************************************************************************
/// \brief	Multiplies a scalar by a vector.
/// \param	scalar	a scalar double.
/// \param	vec		a vector double.
/// \return a vector which is a scalar multiplied by a vector.
std::vector<double> GridUtils::vecmultiply(double scalar, std::vector<double> vec) {

	std::vector<double> result;
	for (size_t i = 0; i < vec.size(); i++) {
		result.push_back( vec[i] * scalar );
	}
	return result;
}
	

// *****************************************************************************
/// \brief	Creates a linearly-spaced vector of values.
/// \param	min	starting value of output vector.
/// \param	max	ending point of output vector.
/// \param	n	number of values in output vector.
/// \return	a vector with n uniformly spaced values between min and max.
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

// *****************************************************************************
/// \brief	Creates a linearly-spaced vector of integers.
/// \param	min	starting value of output vector.
/// \param	max	ending point of output vector.
/// \return	a vector with uniformly spaced integer values between min and max.
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

// *****************************************************************************
/// \brief	Computes the L2 norm using the vector components supplied.
/// \param	val1	first vector component.
/// \param	val2	second vector component.
/// \return	the L2 norm.
double GridUtils::vecnorm( double val1, double val2 )
{
	double result;

	result = sqrt( pow(val1,2) + pow(val2,2) );

	return result;
}

// *****************************************************************************
/// \brief	Computes the L2 norm using the vector components supplied.
/// \param	val1	first vector component.
/// \param	val2	second vector component.
/// \param	val3	third vector component.
/// \return	the L2 norm.
double GridUtils::vecnorm( double val1, double val2, double val3 )
{
	double result;

	result = sqrt( pow(val1,2) + pow(val2,2) + pow(val3,2) );

	return result;
}

// *****************************************************************************
/// \brief	Computes the L2 norm using the vector supplied.
/// \param	vec	old-style C array representing a vector with the same number of
///				number of components as the problem dimension.
/// \return	the L2 norm.
double GridUtils::vecnorm( double vec[L_dims] )
{
	double result;

#if (L_dims == 3)

		result = sqrt( pow(vec[0],2) + pow(vec[1],2) + pow(vec[2],2) );

#else

		result = sqrt( pow(vec[0],2) + pow(vec[1],2) );

#endif

	return result;
}

// *****************************************************************************
/// \brief	Computes the L2 norm using the vector supplied.
/// \param	vec	C++ std::vector.
/// \return	the L2 norm.
double GridUtils::vecnorm( std::vector<double> vec )
{
	double result = 0.0;

	for (size_t d = 0; d < vec.size(); d++) {

		result += pow(vec[d],2);

	}

	return sqrt(result);
}


// *****************************************************************************
/// \brief	Computes the scalar product of two vectors.
/// \param	vec1	a vector.
/// \param	vec2	a second vector.
/// \return	the dot product of the two vectors.
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

// *****************************************************************************
/// \brief	Multiplies matrix A by vector x.
/// \param	A	a matrix represented as a vector or vectors.
/// \param	x	a vector.
/// \return	a vector which is A * x.
std::vector<double> GridUtils::matrix_multiply(const std::vector< std::vector<double> >& A, const std::vector<double>& x) {

	// Check to makes sure dimensions are correct
	if (A[0].size() != x.size()) {
		std::cout << "Error: See Log File" << std::endl;
		*logfile << "Dimension mismatch -- cannot proceed. Exiting." << std::endl;
		int ignore = system("pause");
		exit(LUMA_FAILED);
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

// *****************************************************************************
/// \brief	Gets the indices of the fine site given the coarse site.
///
///			Maps the indices of a coarse grid site to a corresponding fine 
///			site on the level below.
///
/// \param	coarse_i	local i-index of coarse site to be mapped.
/// \param	x_start		local x-index of start of refined region.
/// \param	coarse_j	local j-index of coarse site to be mapped.
/// \param	y_start		local y-index of start of refined region.
/// \param	coarse_k	local k-index of coarse site to be mapped.
/// \param	z_start		local z-index of start of refined region.
/// \return	local indices of the fine grid site.
std::vector<int> GridUtils::getFineIndices(int coarse_i, int x_start, int coarse_j, int y_start, int coarse_k, int z_start) {

	// Initialise result
	std::vector<int> fine_ind;

	// Map indices
	fine_ind.insert(fine_ind.begin(), 2 * (coarse_i - x_start + 1) - 2);
	fine_ind.insert(fine_ind.begin() + 1, 2 * (coarse_j - y_start + 1) - 2);
#if (L_dims == 3)
	fine_ind.insert(fine_ind.begin() + 2, 2 * (coarse_k - z_start + 1) - 2);
#else
	fine_ind.insert(fine_ind.begin() + 2, 0);
#endif

	return fine_ind;
}

// *****************************************************************************
/// \brief	Gets the indices of the coarse site given the fine site.
///
///			Maps the indices of a fine grid site to a corresponding coarse 
///			site on the level above.
///
/// \param	fine_i	local i-index of fine site to be mapped.
/// \param	x_start	local x-index of start of refined region on the grid above.
/// \param	fine_j	local j-index of fine site to be mapped.
/// \param	y_start	local y-index of start of refined region on the grid above.
/// \param	fine_k	local k-index of fine site to be mapped.
/// \param	z_start	local z-index of start of refined region on the grid above.
/// \return	local indices of the coarse grid site.
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
	coarse_ind.insert(coarse_ind.begin(), (fine_i / 2) + x_start);
	coarse_ind.insert(coarse_ind.begin() + 1, (fine_j / 2) + y_start);
#if (L_dims == 3)
	coarse_ind.insert(coarse_ind.begin() + 2, (fine_k / 2) + z_start);
#else
	coarse_ind.insert(coarse_ind.begin() + 2, 0);
#endif

	return coarse_ind;
}

// *****************************************************************************
/// \brief	Gets the opposite lattice direction to the one supplied.
///
///			This is model independent as long as the model directions are 
///			specified such that the oppoiste direction is either one vector on
///			or one vector back in the listing depending on whether the direction
///			supplied is even or odd.
///
/// \param	direction	direction to be reversed.
/// \return	opposite direction in lattice model.
int GridUtils::getOpposite(int direction) {

	int direction_opposite;

	// If rest particle then opposite is simply itself
	if (direction == L_nVels-1) {

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

// *****************************************************************************
/// \brief	Finds out whether halo containng i,j,k links to neighbour rank periodically.
///
///			Checks the receiver layer containing local site i,j,k and determines 
///			from the MPI topology information whether this layer couples to an 
///			adjacent or periodic neighbour rank. I.e. if the neighbour is physically 
///			next to the rank or whether it is actaully at the other side of the domain.
///
/// \param	i	local i-index of recv layer site being queried.
/// \param	j	local j-index of recv layer site being queried.
/// \param	k	local k-index of recv layer site being queried.
/// \param	g	grid on which point being queried resides.
/// \return	boolean answer.
bool GridUtils::isOverlapPeriodic(int i, int j, int k, const GridObj& g) {

	// Local declarations
	int exp_MPI_coords[L_dims], act_MPI_coords[L_dims], MPI_dims[L_dims];
	int shift[3] = {0, 0, 0};

	// Initialise local variables
	MPI_dims[0] = L_Xcores;
	MPI_dims[1] = L_Ycores;
#if (L_dims == 3)
	MPI_dims[2] = L_Zcores;
#endif

	// Define shifts based on which overlap we are on

	// X
	if (GridUtils::isOnRecvLayer(g.XPos[i],eXDirection,eMaximum)) {
		shift[0] = 1;
	} else if (GridUtils::isOnRecvLayer(g.XPos[i],eXDirection,eMinimum)) {
		shift[0] = -1;
	}

	// Y
	if (GridUtils::isOnRecvLayer(g.YPos[j],eYDirection,eMaximum)) {
		shift[1] = 1;
	} else if (GridUtils::isOnRecvLayer(g.YPos[j],eYDirection,eMinimum)) {
		shift[1] = -1;
	}

#if (L_dims == 3)
	// Z
	if (GridUtils::isOnRecvLayer(g.ZPos[k],eZDirection,eMaximum)) {
		shift[2] = 1;
	} else if (GridUtils::isOnRecvLayer(g.ZPos[k],eZDirection,eMinimum)) {
		shift[2] = -1;
	}
#endif

	// Loop over each Cartesian direction
	for (int d = 0; d < L_dims; d++) {
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

// *****************************************************************************
/// \brief	Finds out whether site with supplied index in on the current rank.
/// \param	gi		global i-index of site.
/// \param	gj		global j-index of site.
/// \param	gk		global k-index of site.
/// \param	grid	grid being queried.
/// \return	boolean answer.
bool GridUtils::isOnThisRank(int gi, int gj, int gk, const GridObj& grid) {
	
	auto found_x = std::find(grid.XInd.begin(), grid.XInd.end(), gi);
	auto found_y = std::find(grid.YInd.begin(), grid.YInd.end(), gj);
	auto found_z = std::find(grid.ZInd.begin(), grid.ZInd.end(), gk);

	if (	found_x != grid.XInd.end() && found_y != grid.YInd.end()


#if (L_dims == 3)
		&& found_z != grid.ZInd.end()
#endif
		) {

		return true;

	} else {

		return false;
	}

}

// *****************************************************************************
/// \brief	Finds out whether global index can be found on the current rank.
/// \param	gl		global index (i,j or k).
/// \param	xyz		cartesian direction of interest.
/// \param	grid	grid being queried.
/// \return	boolean answer.
bool GridUtils::isOnThisRank(int gl, enum eCartesianDirection xyz, const GridObj& grid) {

	switch (xyz) {

	case eXDirection:
		{
		// X direction
		auto found_x = std::find(grid.XInd.begin(), grid.XInd.end(), gl);

		if (found_x != grid.XInd.end()) return true;		
		else return false;
		}

	case eYDirection:
		{
		// Y direction
		auto found_y = std::find(grid.YInd.begin(), grid.YInd.end(), gl);

		if (found_y != grid.YInd.end()) return true;		
		else return false;
		}

	case eZDirection:
		{
		// Z direction
		auto found_z = std::find(grid.ZInd.begin(), grid.ZInd.end(), gl);

		if (found_z != grid.ZInd.end()) return true;		
		else return false;
		}

	}

	return false;

}

// *****************************************************************************
/// \brief	Finds out whether specified refined region is on the grid provided.
/// \param	pGrid	parent grid at appropriate level.
/// \param	RegNum	region number desired.
/// \return	boolean answer.
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

#if (L_dims == 3)
	// Loop over over Z range
	for (size_t k = RefZstart[pGrid.level][RegNum]; k <= RefZend[pGrid.level][RegNum]; k++) {
		auto found_k = std::find(pGrid.ZInd.begin(), pGrid.ZInd.end(), k);
		if (found_k != pGrid.ZInd.end()) break;
		else if (k == RefZend[pGrid.level][RegNum] && found_k == pGrid.ZInd.end()) return false;
	}
#endif

	return true;
	
}

// ****************************************************************************
/// \brief	Get a pointer to a given grid in the hierarchy.
///
///			Takes a NULL pointer by reference and updates it when matching grid 
///			is found in hierarchy on this rank. If grid not found, pointer is 
///			returned without change and stays NULL. Can be used to test for the 
///			existence of a grid on a rank by passing in a NULL pointer and 
///			checking if a NULL pointer is returned.
///
/// \param		Grids	x-position of site.
/// \param		level	y-position of site.
/// \param		region	z-position of site.
/// \param[out] ptr		pointer containing address of grid in hierarchy.
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

// ****************************************************************************
/// \brief	Check whether site is on an inner (sender) halo.
///
///			Wrapper which checks every halo region of the rank for intersection
///			with supplied site position.
///
/// \param	pos_x	x-position of site.
/// \param	pos_y	y-position of site.
/// \param	pos_z	z-position of site.
/// \return	boolean answer.
bool GridUtils::isOnSenderLayer(double pos_x, double pos_y, double pos_z) {

	/* Checks whether X,Y or Z are on a sender layer first.
	 * If true, ensures that neither of the other two are on a receiver layer.
	 * If true then site must be on a sender layer. */

	if (
	(
		// X on sender
		(GridUtils::isOnSenderLayer(pos_x,eXDirection,eMinimum) || GridUtils::isOnSenderLayer(pos_x,eXDirection,eMaximum)) &&

		// Y and Z not recv
		(!GridUtils::isOnRecvLayer(pos_y,eYDirection,eMinimum) && !GridUtils::isOnRecvLayer(pos_y,eYDirection,eMaximum)
#if (L_dims == 3)
		&& !GridUtils::isOnRecvLayer(pos_z,eZDirection,eMinimum) && !GridUtils::isOnRecvLayer(pos_z,eZDirection,eMaximum)
#endif
		)
	) || (
		// Y on sender
		(GridUtils::isOnSenderLayer(pos_y,eYDirection,eMinimum) || GridUtils::isOnSenderLayer(pos_y,eYDirection,eMaximum)) &&

		// X and Z not recv
		(!GridUtils::isOnRecvLayer(pos_x,eXDirection,eMinimum) && !GridUtils::isOnRecvLayer(pos_x,eXDirection,eMaximum)
#if (L_dims == 3)
		&& !GridUtils::isOnRecvLayer(pos_z,eZDirection,eMinimum) && !GridUtils::isOnRecvLayer(pos_z,eZDirection,eMaximum)

		)
	) || (
		// Z on sender
		(GridUtils::isOnSenderLayer(pos_z,eZDirection,eMinimum) || GridUtils::isOnSenderLayer(pos_z,eZDirection,eMaximum)) &&

		// X and Y not recv
		(!GridUtils::isOnRecvLayer(pos_x,eXDirection,eMinimum) && !GridUtils::isOnRecvLayer(pos_x,eXDirection,eMaximum)	&& 
		!GridUtils::isOnRecvLayer(pos_y,eYDirection,eMinimum) && !GridUtils::isOnRecvLayer(pos_y,eYDirection,eMaximum)
#endif
		)
	)
	) {
		return true;
	}

	return false;

}

// ****************************************************************************
/// \brief	Check whether site is on an inner (sender) halo.
///
///			Wrapper available which checks every halo. This method only checks 
///			the halo specified by the Cartesian direction and whether it is the 
///			left/bottom/front (minimum) or right/top/back (maximum) edge of the 
///			block.
///
/// \param	site_position	position of site.
/// \param	dir				cartesian direction.
/// \param maxmin			choice of edge in given direction.
/// \return	boolean answer.
bool GridUtils::isOnSenderLayer(double site_position, enum eCartesianDirection dir, enum eMinMax maxmin) {

	// Get instance of MPI manager
	MpiManager* mpim = MpiManager::getInstance();
	
	// Do checks on rank inner overlap regions dictated by the coarsest rank
	if (dir == eXDirection) {

		if (maxmin == eMaximum) {
			if (site_position > mpim->sender_layer_pos.X[2] && site_position < mpim->sender_layer_pos.X[3] ) return true;

		} else if (maxmin == eMinimum) {
			if (site_position > mpim->sender_layer_pos.X[0] && site_position < mpim->sender_layer_pos.X[1] ) return true;

		}

	} else if (dir == eYDirection) {

		if (maxmin == eMaximum) {
			if (site_position > mpim->sender_layer_pos.Y[2] && site_position < mpim->sender_layer_pos.Y[3] ) return true;

		} else if (maxmin == eMinimum) {
			if (site_position > mpim->sender_layer_pos.Y[0] && site_position < mpim->sender_layer_pos.Y[1] ) return true;

		}

	} else if (dir == eZDirection) {

		if (maxmin == eMaximum) {
			if (site_position > mpim->sender_layer_pos.Z[2] && site_position < mpim->sender_layer_pos.Z[3] ) return true;

		} else if (maxmin == eMinimum) {
			if (site_position > mpim->sender_layer_pos.Z[0] && site_position < mpim->sender_layer_pos.Z[1] ) return true;

		}

	} else {

		// Invalid string indicating direction
		std::cout << "Error: See Log File" << std::endl;
		*logfile << "Invalid direction specified in GridUtils::isOnSenderLayer() / GridUtils::isOnRecvLayer(). Exiting." << std::endl;
		exit(LUMA_FAILED);

	}


	return false;

}

// ****************************************************************************
/// \brief	Check whether site is on an outer (receiver) halo.
///
///			Wrapper which checks every halo region of the rank for intersection
///			with supplied site position.
///
/// \param	pos_x	x-position of site.
/// \param	pos_y	y-position of site.
/// \param	pos_z	z-position of site.
/// \return	boolean answer.
bool GridUtils::isOnRecvLayer(double pos_x, double pos_y, double pos_z) {

	/* If any of the coordinates is in the receiver layer then the site will always 
	 * be a receiver layer site regardless of the other coordinates so the logic is
	 * simple. */

	if (	GridUtils::isOnRecvLayer(pos_x,eXDirection,eMinimum) || GridUtils::isOnRecvLayer(pos_x,eXDirection,eMaximum) ||
			GridUtils::isOnRecvLayer(pos_y,eYDirection,eMinimum) || GridUtils::isOnRecvLayer(pos_y,eYDirection,eMaximum)
#if (L_dims == 3)
			|| GridUtils::isOnRecvLayer(pos_z,eZDirection,eMinimum) || GridUtils::isOnRecvLayer(pos_z,eZDirection,eMaximum)
#endif
		) {
			return true;
	}

	return false;

}

// ****************************************************************************
/// \brief	Check whether site is on an outer (receiver) halo.
///
///			Wrapper available which checks every halo. This method only checks 
///			the halo specified by the Cartesian direction and whether it is the 
///			left/bottom/front (minimum) or right/top/back (maximum) edge of the 
///			block.
///
/// \param	site_position	position of site.
/// \param	dir				cartesian direction.
/// \param maxmin			choice of edge in given direction.
/// \return	boolean answer.
bool GridUtils::isOnRecvLayer(double site_position, enum eCartesianDirection dir, enum eMinMax maxmin) {

	// Get instance of MPI manager
	MpiManager* mpim = MpiManager::getInstance();
	
	// Do checks on rank outer overlap regions dictated by the coarsest rank
	if (dir == eXDirection) {

		if (maxmin == eMaximum) {
			if (site_position > mpim->recv_layer_pos.X[2] && site_position < mpim->recv_layer_pos.X[3] ) return true;

		} else if (maxmin == eMinimum) {
			if (site_position > mpim->recv_layer_pos.X[0] && site_position < mpim->recv_layer_pos.X[1] ) return true;

		}

	} else if (dir == eYDirection) {

		if (maxmin == eMaximum) {
			if (site_position > mpim->recv_layer_pos.Y[2] && site_position < mpim->recv_layer_pos.Y[3] ) return true;

		} else if (maxmin == eMinimum) {
			if (site_position > mpim->recv_layer_pos.Y[0] && site_position < mpim->recv_layer_pos.Y[1] ) return true;

		}
		
	} else if (dir == eZDirection) {

		if (maxmin == eMaximum) {
			if (site_position > mpim->recv_layer_pos.Z[2] && site_position < mpim->recv_layer_pos.Z[3] ) return true;

		} else if (maxmin == eMinimum) {
			if (site_position > mpim->recv_layer_pos.Z[0] && site_position < mpim->recv_layer_pos.Z[1] ) return true;

		}

	} else {

		// Invalid string indicating direction
		std::cout << "Error: See Log File" << std::endl;
		*logfile << "Invalid direction specified in GridUtils::isOnSenderLayer() / GridUtils::isOnRecvLayer(). Exiting." << std::endl;
		exit(LUMA_FAILED);

	}

	return false;

}

// ****************************************************************************
/// \brief	Create output directory.
///
///			Compatible with both Windows and Linux. Filename and path passed as 
///			a single string. Returns 9 if the directory creation was not 
///			attempted due to not being rank 0. Returns platform specific codes
///			for everything else.
///
/// \param	path_str	full path and filename as string.
/// \return indicator of status of action.
int GridUtils::createOutputDirectory(std::string path_str) {

	int result = 9; // Return code of directory creation

	// Create output directory if it does not already exist
	std::string command = "mkdir -p " + path_str;

	// Only get rank 0 to create output directory
#ifdef L_BUILD_FOR_MPI

	#ifdef _WIN32   // Running on Windows
		if (MpiManager::my_rank == 0)
			result = CreateDirectoryA((LPCSTR)path_str.c_str(), NULL);
	#else   // Running on Unix system
		if (MpiManager::my_rank == 0)
			result = system(command.c_str());
	#endif // _WIN32

#else // L_BUILD_FOR_MPI

	#ifdef _WIN32   // Running on Windows
		result = CreateDirectoryA((LPCSTR)path_str.c_str(), NULL);
	#else   // Running on Unix system
		result = system(command.c_str());
	#endif // _WIN32

#endif // L_BUILD_FOR_MPI

	return result;
}

// ****************************************************************************
/// \brief	Tests whether a site is on a given grid.
/// \param	i	local i-index.
/// \param	j	local j-index.
/// \param	k	local k-index.
/// \param	g	grid on which to check.
/// \return boolean answer.
bool GridUtils::isOffGrid(int i, int j, int k, GridObj& g) {

	if (	(i >= g.N_lim || i < 0) ||
			(j >= g.M_lim || j < 0) ||
			(k >= g.K_lim || k < 0)
			) {
				return true;

#ifdef L_BUILD_FOR_MPI
		// When using MPI, equivalent to off-grid is when destination is in 
		// periodic recv layer with periodic boundaries disabled.
		} else if ( GridUtils::isOnRecvLayer(g.XPos[i],g.YPos[j],g.ZPos[k]) 
			&& GridUtils::isOverlapPeriodic(i,j,k,g) ) {

			return true;		

#endif	// L_BUILD_FOR_MPI

		}

		return false;
}

// ****************************************************************************
/// \brief	Get global voxel indices
///
///			Will return the voxel indices of the nearest voxel on the lattice 
///			for a given point in global space. Assumes that the physical lattice 
///			size is 1 unit big which relies on the position being scaled. If not,
///			method acts as a natural voxel grid filter.
///
/// \param	x	global x-position.
/// \param	y	global y-position.
/// \param	z	global z-position.
/// \return vector of indices of the nearest voxel.
std::vector<int> GridUtils::getVoxInd(double x, double y, double z) {

	std::vector<int> vox;


	// Check whether point is closer to the floor or the ceiling value and add correct one
	if (x - (int)std::floor(x) > 0.5) vox.push_back((int)std::ceil(x));
	else vox.push_back((int)std::floor(x));

	if (y - (int)std::floor(y) > 0.5) vox.push_back((int)std::ceil(y));
	else vox.push_back((int)std::floor(y));

	if (z - (int)std::floor(z) > 0.5) vox.push_back((int)std::ceil(z));
	else vox.push_back((int)std::floor(z));

	return vox;

}

// ****************************************************************************
/// \brief	Get global voxel index
///
////		Will return the voxel indices of the nearest voxel on the lattice
///			for a given point in global space. Assumes that the physical lattice 
///			size is 1 unit big which relies on the position being scaled. If not,
///			method acts as a natural voxel grid filter.
///
/// \param	p	global position.
/// \return corresponding global index.
int GridUtils::getVoxInd(double p) {

	if (p - (int)std::floor(p) > 0.5) return (int)std::ceil(p);
	else return (int)std::floor(p);

}
