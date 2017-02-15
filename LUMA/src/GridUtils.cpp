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
#include "../inc/GridObj.h"

// Mappings of directions for specular reflection: col == normal direction, row == velocity
// Constants so initialise outside the class but in scope.
#if (L_DIMS == 3)

// THIS IS WRONG AS WE HAVE CHANGED TO D3Q27 NOW!!!!
const int GridUtils::dir_reflect[L_DIMS * 2][L_NUM_VELS] = 
	{
		{1, 0, 2, 3, 4, 5, 9, 8, 7, 6, 10, 11, 12, 13, 16, 17, 14, 15, 18}, 
		{1, 0, 2, 3, 4, 5, 9, 8, 7, 6, 10, 11, 12, 13, 16, 17, 14, 15, 18},
		{0, 1, 3, 2, 4, 5, 8, 9, 6, 7, 13, 12, 11, 10, 14, 15, 16, 17, 18},
		{0, 1, 3, 2, 4, 5, 8, 9, 6, 7, 13, 12, 11, 10, 14, 15, 16, 17, 18},
		{0, 1, 2, 3, 5, 4, 6, 7, 8, 9, 12, 13, 10, 11, 17, 16, 15, 14, 18},
		{0, 1, 2, 3, 5, 4, 6, 7, 8, 9, 12, 13, 10, 11, 17, 16, 15, 14, 18}
	};
#else
const int GridUtils::dir_reflect[L_DIMS * 2][L_NUM_VELS] = 
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
double GridUtils::vecnorm( double vec[L_DIMS] )
{
	double result;

#if (L_DIMS == 3)

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
		L_ERROR("Dimension mismatch -- cannot proceed. Exiting.", logfile);
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
#if (L_DIMS == 3)
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
#if (L_DIMS == 3)
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
	if (direction == L_NUM_VELS-1) {

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
bool GridUtils::isOverlapPeriodic(int i, int j, int k, GridObj const & g) {

	// Local declarations
	int exp_rank_coords[L_DIMS], act_rank_coords[L_DIMS], dimensions[L_DIMS];
	int shift[3] = {0, 0, 0};

	// Get MpiManager instance
	MpiManager *mpim = MpiManager::getInstance();

	// Initialise local variables
	dimensions[0] = L_MPI_XCORES;
	dimensions[1] = L_MPI_YCORES;
#if (L_DIMS == 3)
	dimensions[2] = L_MPI_ZCORES;
#endif

	// Define shifts based on which overlap we are on

	// X
	if (GridUtils::isOnRecvLayer(g.XPos[i],eXMax)) {
		shift[0] = 1;
	} else if (GridUtils::isOnRecvLayer(g.XPos[i],eXMin)) {
		shift[0] = -1;
	}

	// Y
	if (GridUtils::isOnRecvLayer(g.YPos[j],eYMax)) {
		shift[1] = 1;
	} else if (GridUtils::isOnRecvLayer(g.YPos[j],eYMin)) {
		shift[1] = -1;
	}

#if (L_DIMS == 3)
	// Z
	if (GridUtils::isOnRecvLayer(g.ZPos[k],eZMax)) {
		shift[2] = 1;
	} else if (GridUtils::isOnRecvLayer(g.ZPos[k],eZMin)) {
		shift[2] = -1;
	}
#endif

	// Loop over each Cartesian direction
	for (int d = 0; d < L_DIMS; d++) {
		// Define expected (non-periodic) MPI coordinates of neighbour rank
		exp_rank_coords[d] = mpim->rank_coords[d] + shift[d];

		// Define actual MPI coordinates of neighbour rank (accounting for periodicity)
		act_rank_coords[d] = (mpim->rank_coords[d] + shift[d] + dimensions[d]) % dimensions[d];

		// If there is a difference then rank is periodically linked to its neighbour and the overlap
		// site is from a periodic rank so return early
		if (exp_rank_coords[d] != act_rank_coords[d]) return true;
	}
	

	// Expected and actual are the same so the neighbour rank is not periodically linked and the overlap site
	// is from an adjacent neighbour not a periodic one.
	return false;

}

// *****************************************************************************
/// \brief	Finds out whether site with supplied position is on the current rank.
///
///			Will return true if the site is in the halo as well (send or recv).
///			Location information provided to indicate where point is. Returns
///			eNone enumeration if not request or if query is false. If a grid is 
///			supplied, will only return true if site is on the grid supplied.
///			If you want to exclude the sites that belong to the halo you can
///			call isOnRecvLayer() or isOnSenderLayer() on the same site.
///
/// \param		x		x-position of site.
/// \param		y		y-position of site.
/// \param		z		z-position of site.
/// \param[out]	pos		pointer to the start of a vector in which local indices are returned.
/// \param		grid	grid being queried.
/// \param[out]	loc		description of the location of the point.
/// \return	boolean answer.
bool GridUtils::isOnThisRank(double x, double y, double z, eLocationOnRank *loc, 
	GridObj const * const grid, std::vector<int> *pos) {
	
	// Initialise result
	bool result = false;
	bool new_vector = false;

#ifdef L_BUILD_FOR_MPI
	MpiManager *mpim = MpiManager::getInstance();
	int rank = GridUtils::safeGetRank();
	
	// Check whether point within the edges of grid core
	if (
		mpim->rank_core_edge[eXMin][rank] <= x && x < mpim->rank_core_edge[eXMax][rank] &&
		mpim->rank_core_edge[eYMin][rank] <= y && y < mpim->rank_core_edge[eYMax][rank]

#if (L_DIMS == 3)
		&&
		mpim->rank_core_edge[eZMin][rank] <= z && z < mpim->rank_core_edge[eZMax][rank]
#endif		
		)
	{
		if (loc != nullptr) *loc = eCore;
		result = true;
	}

	// Check whether point within receiver layers
	else if (GridUtils::isOnRecvLayer(x, y, z))
	{
		if (loc != nullptr) *loc = eHalo;
		result = true;
	}

#else

	// In serial always on Core or not on grid at all
	GridManager *gm = GridManager::getInstance();

	// Check with coarsest grid limits
	if (
		gm->global_edges[eXMin][0] <= x && x < gm->global_edges[eXMax][0] &&
		gm->global_edges[eYMin][0] <= y && y < gm->global_edges[eYMax][0]
#if (L_DIMS == 3)
		&&
		gm->global_edges[eZMin][0] <= z && z < gm->global_edges[eZMax][0]
#endif
		)
	{
		if (loc != nullptr) *loc = eCore;
		result = true;
	}


#endif	// L_BUILD_FOR_MPI

	// Not on either core or halo so not on grid
	else
	{
		if (loc != nullptr) *loc = eNone;
		result = false;
		return result;
	}

	// If a grid is supplied then get its index on the grid
	if (grid != nullptr)
	{
		// If location vector is null then create one and remember to delete it
		if (pos == nullptr) {
			new_vector = true;
			pos = new std::vector<int>();
		}

		// Get voxel index
		getEnclosingVoxel(x, y, z, grid, pos);

		// If voxel indices are not within the grid range then off-grid
		if (GridUtils::isOffGrid((*pos)[0], (*pos)[1], (*pos)[2], grid)) result = false;
	}
	
	if (new_vector) delete pos;
	return result;
}

// *****************************************************************************
/// \brief	Finds out whether the supplied position can be found on the current rank.
///
///			Direction-specific version of the overload.
///
/// \param	xyz		position (x, y or z)
/// \param	dir		cartesian direction of interest (x, y or z).
/// \param[out]	loc	description of the location of the point.
/// \param	grid	grid being queried.
/// \param[out]	pos	the local index of the found site.
/// \return	boolean answer.
bool GridUtils::isOnThisRank(double xyz, eCartesianDirection dir, 
	eLocationOnRank *loc, GridObj const * const grid, int *pos) {

	// Initialise result
	bool result = false;
	int lims[2];

	// Specific activity per direction
	switch (dir) {

	case eXDirection:
		lims[0] = eXMin;
		lims[1] = eXMax;
		break;

	case eYDirection:
		lims[0] = eYMin;
		lims[1] = eYMax;
		break;

	case eZDirection:
		lims[0] = eZMin;
		lims[1] = eZMax;
		break;
	}


#ifdef L_BUILD_FOR_MPI

	MpiManager *mpim = MpiManager::getInstance();
	int rank = GridUtils::safeGetRank();

	// Is on core?
	if (mpim->rank_core_edge[lims[0]][rank] <= xyz && xyz < mpim->rank_core_edge[lims[1]][rank])
	{
		if (loc != nullptr) *loc = eCore;
		result = true;
	}

	// Is on halo?
	else if (
		GridUtils::isOnRecvLayer(xyz, static_cast<eCartMinMax>(lims[0])) || 
		GridUtils::isOnRecvLayer(xyz, static_cast<eCartMinMax>(lims[1]))
		)
	{
		if (loc != nullptr) *loc = eHalo;
		result = true;
	}

#else

	// In serial always on Core or not on grid at all
	GridManager *gm = GridManager::getInstance();

	// Check with coarsest grid limits
	if (gm->global_edges[lims[0]][0] <= xyz && xyz < gm->global_edges[lims[1]][0])
	{
		if (loc != nullptr) *loc = eCore;
		result = true;
	}

#endif	// L_BUILD_FOR_MPI

	// Not on grid
	else
	{
		if (loc != nullptr) *loc = eNone;
		return false;
	}

	// If a grid is supplied then get its index on the grid
	if (grid != nullptr)
	{
		GridUtils::getEnclosingVoxel(xyz, grid, dir, pos);
	}

	return result;

}

// *****************************************************************************
/// \brief	Finds out whether all or part of specified refined region intersects
///			with the space occupied by the grid provided.
///
///			Prinicpal use is for sub-grid initialisation to determine
///			whether a sub-grid needs adding or not. This decision is made based
///			on whether any part of the grid is covered by the discrete voxels
///			of existing grids on the rank.
///
/// \param	pGrid	parent grid at appropriate level.
/// \param	RegNum	region number desired.
/// \return	boolean answer.
bool GridUtils::intersectsRefinedRegion(GridObj const & pGrid, int RegNum) {


	/* To check this condition we need at least one of the voxel centres on the 
	 * grid to have a position within the edges of the refined region. The GM
	 * has the grid edge information against which we can check the voxels. */
	GridManager *gm = GridManager::getInstance();
	bool result = false;
	size_t i;
		
	// X //
	// Check to see whether the x-edges on rank	
	for (i = 0; i < pGrid.N_lim; ++i)
	{
		// If voxel centre intersects range of refined region
		if (pGrid.XPos[i] > gm->global_edges[eXMin][(pGrid.level + 1) + RegNum * L_NUM_LEVELS] &&
			pGrid.XPos[i] < gm->global_edges[eXMax][(pGrid.level + 1) + RegNum * L_NUM_LEVELS])
		{
			result = true;	// Possibility of intersection
			break; 
		}
	}
	// If result has not be changed to true then no X intersect so return false result
	if (!result) return result;
	else result = false;	// If not reset search result flag

	// Y //
	for (i = 0; i < pGrid.M_lim; ++i)
	{
		if (pGrid.YPos[i] > gm->global_edges[eYMin][(pGrid.level + 1) + RegNum * L_NUM_LEVELS] &&
			pGrid.YPos[i] < gm->global_edges[eYMax][(pGrid.level + 1) + RegNum * L_NUM_LEVELS])
		{
			result = true;
			break;
		}
	}
	if (!result) return result;
	
#if (L_DIMS == 3)
	else result = false;

	// Z //
	for (i = 0; i < pGrid.K_lim; ++i)
	{
		if (pGrid.ZPos[i] > gm->global_edges[eZMin][(pGrid.level + 1) + RegNum * L_NUM_LEVELS] &&
			pGrid.ZPos[i] < gm->global_edges[eZMax][(pGrid.level + 1) + RegNum * L_NUM_LEVELS])
		{
			result = true;
			break;
		}
	}
	if (!result) return result;
#endif

	return result;
	
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
/// \param		Grids	constant pointer to the grid at which to start searching.
/// \param		level	level desired.
/// \param		region	region desried.
/// \param[out] ptr		reference to pointer where address of grid matching in hierarchy will be assigned.
void GridUtils::getGrid(GridObj* const Grids, int level, int region, GridObj*& ptr) {

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
		(GridUtils::isOnSenderLayer(pos_x,eXMin) || GridUtils::isOnSenderLayer(pos_x,eXMax)) &&

		// Y and Z not recv
		(
		!GridUtils::isOnRecvLayer(pos_y,eYMin) && !GridUtils::isOnRecvLayer(pos_y,eYMax)
#if (L_DIMS == 3)
		&& !GridUtils::isOnRecvLayer(pos_z,eZMin) && !GridUtils::isOnRecvLayer(pos_z,eZMax)
#endif
		)
	) || (
		// Y on sender
		(GridUtils::isOnSenderLayer(pos_y,eYMin) || GridUtils::isOnSenderLayer(pos_y,eYMax)) &&

		// X and Z not recv
		(
		!GridUtils::isOnRecvLayer(pos_x,eXMin) && !GridUtils::isOnRecvLayer(pos_x,eXMax)
#if (L_DIMS == 3)
		&& !GridUtils::isOnRecvLayer(pos_z,eZMin) && !GridUtils::isOnRecvLayer(pos_z,eZMax)

		)
	) || (
		// Z on sender
		(GridUtils::isOnSenderLayer(pos_z,eZMin) || GridUtils::isOnSenderLayer(pos_z,eZMax)) &&

		// X and Y not recv
		(
		!GridUtils::isOnRecvLayer(pos_x,eXMin) && !GridUtils::isOnRecvLayer(pos_x,eXMax)	&& 
		!GridUtils::isOnRecvLayer(pos_y,eYMin) && !GridUtils::isOnRecvLayer(pos_y,eYMax)
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
/// \param	edge			combination of cartesian direction and choice of edge.
/// \return	boolean answer.
bool GridUtils::isOnSenderLayer(double site_position, enum eCartMinMax edge) {

	// Get instance of MPI manager
	MpiManager* mpim = MpiManager::getInstance();
	
	// Do checks on rank inner overlap regions dictated by the coarsest rank
	if (edge == eXMax)
	{
		if (site_position >= mpim->sender_layer_pos.X[eRightMin] && site_position < mpim->sender_layer_pos.X[eRightMax]) return true;
	}
	else if (edge == eXMin)
	{
		if (site_position >= mpim->sender_layer_pos.X[eLeftMin] && site_position < mpim->sender_layer_pos.X[eLeftMax]) return true;
	}
	else if (edge == eYMax)
	{
		if (site_position >= mpim->sender_layer_pos.Y[eRightMin] && site_position < mpim->sender_layer_pos.Y[eRightMax]) return true;
	}
	else if (edge == eYMin)
	{
		if (site_position >= mpim->sender_layer_pos.Y[eLeftMin] && site_position < mpim->sender_layer_pos.Y[eLeftMax]) return true;
	}
	else if (edge == eZMax)
	{
		if (site_position >= mpim->sender_layer_pos.Z[eRightMin] && site_position < mpim->sender_layer_pos.Z[eRightMax]) return true;
	}
	else if (edge == eZMin)
	{
		if (site_position >= mpim->sender_layer_pos.Z[eLeftMin] && site_position < mpim->sender_layer_pos.Z[eLeftMax]) return true;
	}
	else
	{
		// Invalid string indicating direction
		L_ERROR("Invalid direction specified in GridUtils::isOnSenderLayer(). Exiting.", logfile);
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

	if (
		GridUtils::isOnRecvLayer(pos_x,eXMin) || GridUtils::isOnRecvLayer(pos_x,eXMax) ||
		GridUtils::isOnRecvLayer(pos_y,eYMin) || GridUtils::isOnRecvLayer(pos_y,eYMax)
#if (L_DIMS == 3)
		||
		GridUtils::isOnRecvLayer(pos_z,eZMin) || GridUtils::isOnRecvLayer(pos_z,eZMax)
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
/// \param	edge			combination of cartesian direction and choice of edge.
/// \return	boolean answer.
bool GridUtils::isOnRecvLayer(double site_position, enum eCartMinMax edge) {

	// Get instance of MPI manager
	MpiManager* mpim = MpiManager::getInstance();
	
	// Do checks on rank outer overlap regions dictated by the coarsest rank
	if (edge == eXMax)
	{
		if (site_position >= mpim->recv_layer_pos.X[eRightMin] && site_position < mpim->recv_layer_pos.X[eRightMax]) return true;
	}
	else if (edge == eXMin)
	{
		if (site_position >= mpim->recv_layer_pos.X[eLeftMin] && site_position < mpim->recv_layer_pos.X[eLeftMax]) return true;
	}
	else if (edge == eYMax)
	{
		if (site_position >= mpim->recv_layer_pos.Y[eRightMin] && site_position < mpim->recv_layer_pos.Y[eRightMax]) return true;
	}
	else if (edge == eYMin)
	{
		if (site_position >= mpim->recv_layer_pos.Y[eLeftMin] && site_position < mpim->recv_layer_pos.Y[eLeftMax]) return true;
	}
	else if (edge == eZMax)
	{
		if (site_position >= mpim->recv_layer_pos.Z[eRightMin] && site_position < mpim->recv_layer_pos.Z[eRightMax]) return true;
	}
	else if (edge == eZMin)
	{
		if (site_position >= mpim->recv_layer_pos.Z[eLeftMin] && site_position < mpim->recv_layer_pos.Z[eLeftMax]) return true;
	}
	else
	{
		// Invalid string indicating direction
		L_ERROR("Invalid direction specified in GridUtils::isOnRecvLayer(). Exiting.", logfile);
	}

	return false;

}

// ****************************************************************************
/// \brief	Check whether site is on a TL
///
///			Wrapper which checks every possible TL location on the grid supplied.
///
/// \param	pos_x	x-position of site.
/// \param	pos_y	y-position of site.
/// \param	pos_z	z-position of site.
/// \param	grid	given grid on which to check.
/// \return	boolean answer.
bool GridUtils::isOnTransitionLayer(double pos_x, double pos_y, double pos_z, GridObj const  * const grid)
{

	if (
		GridUtils::isOnTransitionLayer(pos_x, eXMin, grid) ||
		GridUtils::isOnTransitionLayer(pos_x, eXMax, grid) ||
		GridUtils::isOnTransitionLayer(pos_y, eYMin, grid) ||
		GridUtils::isOnTransitionLayer(pos_y, eYMax, grid)
#if (L_DIMS == 3)
		||

		GridUtils::isOnTransitionLayer(pos_z, eZMin, grid) ||
		GridUtils::isOnTransitionLayer(pos_z, eZMax, grid)
#endif
		) return true;

	return false;

}

// ****************************************************************************
/// \brief	Check whether site is on a specific TL (to upper).
///
///			Wrapper available which checks every TL. This method only checks 
///			the TL specified by the Cartesian direction and whether it is the 
///			left/bottom/front (minimum) or right/top/back (maximum) edge of the 
///			supplied grid.
///
/// \param	position	position of point.
/// \param	edge			combination of cartesian direction and choice of edge.
/// \param	grid		given grid on which to check.
/// \return	boolean answer.
bool GridUtils::isOnTransitionLayer(double position, enum eCartMinMax edge, GridObj const * grid) {

	// L0 has no TL so always return false
	if (grid->level == 0) return false;

	// Get GM instance
	GridManager *gm = GridManager::getInstance();

	/* TL defined by absolute position of inner edge relative to the refined limits (outer edge).
	 * The actual location of this region is obtainable through by adding / subtracting the width
	 * of the TL from the grid edges which the MPIM holds from initialisation. */
	double left_edge, right_edge;

	// Get indexing
	int idx = grid->level + grid->region_number * L_NUM_LEVELS;

	// Switch on direction
	if (edge == eXMax)
	{
		// If no TL as indicated at grid initialisation time then return
		if (!gm->subgrid_tlayer_key[eXMax][idx - 1]) return false;

		// Define TL
		right_edge = gm->global_edges[eXMax][idx];
		left_edge = right_edge - 2.0 * grid->dh;
	}
	else if (edge == eXMin)
	{
		if (!gm->subgrid_tlayer_key[eXMin][idx - 1]) return false;
		left_edge = gm->global_edges[eXMin][idx];
		right_edge = left_edge + 2.0 * grid->dh;
	}
	else if (edge == eYMax)
	{
		if (!gm->subgrid_tlayer_key[eYMax][idx - 1]) return false;
		right_edge = gm->global_edges[eYMax][idx];
		left_edge = right_edge - 2.0 * grid->dh;
	}
	else if (edge == eYMin)
	{
		if (!gm->subgrid_tlayer_key[eYMin][idx - 1]) return false;
		left_edge = gm->global_edges[eYMin][idx];
		right_edge = left_edge + 2.0 * grid->dh;
	}
	else if (edge == eZMax)
	{
		if (!gm->subgrid_tlayer_key[eZMax][idx - 1]) return false;
		right_edge = gm->global_edges[eZMax][idx];
		left_edge = right_edge - 2.0 * grid->dh;

	}
	else if (edge == eZMin)
	{
		if (!gm->subgrid_tlayer_key[eZMin][idx - 1]) return false;
		left_edge = gm->global_edges[eZMin][idx];
		right_edge = left_edge + 2.0 * grid->dh;

	}
	else
	{
		// Invalid string indicating direction
		L_ERROR("Invalid direction specified in GridUtils::isOnTL(). Exiting.", logfile);
	}

	if (position >= left_edge && position < right_edge) return true;
	else return false;

}

// ****************************************************************************
/// \brief	Create output directory.
///
///			Compatible with both Windows and Linux. Filename and path passed as 
///			a single string. Returns nothing at the moment.
///
/// \param	path_str	full path and filename as string.
/// \return indicator of status of action.
void GridUtils::createOutputDirectory(std::string path_str) {

	// Create output directory if it does not already exist
	std::string command = "mkdir -p " + path_str;

#ifdef _WIN32   // Running on Windows
			CreateDirectoryA((LPCSTR)path_str.c_str(), NULL);
#else
			system(command.c_str());
#endif // _WIN32

	return;	// TODO: Handle directory creation errors
}

// ****************************************************************************
/// \brief	Reads coordinates and velocity data from a file. 
///
///			The file format is x_coord  y_coord  z_coord ux uy uz 
///			It expects the file to have a column for z_coord and a column for uz even with L_DIMS = 2
///
/// \param	path_str	full path and filename as string.
/// \param  x_coord     vector where the x coordinate of each point will be stored
/// \param  y_coord     vector where the y coordinate of each point will be stored
/// \param  z_coord     vector where the z coordinate of each point will be stored
/// \param  ux          vector where the x component of the velocity will be stored
/// \param  uy          vector where the y component of the velocity will be stored
/// \param  uz          vector where the z component of the velocity will be stored
void GridUtils::readVelocityFromFile(std::string path_str, IVector<double>& x_coord, IVector<double>& y_coord, IVector<double>& z_coord, IVector<double>& ux, IVector<double>& uy, IVector<double>& uz)
{
	// Indicate to log
	*GridUtils::logfile << "Reading file..." << std::endl;

	double tmp;
	
	// Buffer information from file
	std::ifstream datafile;
	datafile.open(path_str, std::ios::in);
	if (!datafile.is_open()) {
		// Error opening file
		L_ERROR("Cannot open velocity data file named " + path_str + ". Exiting.", GridUtils::logfile);

	}
	else {

		std::string line_in;	// String to store line
		std::istringstream iss;	// Buffer stream

		while (!datafile.eof()) {

			// Get line and put in buffer
			std::getline(datafile, line_in, '\n');
			iss.str(line_in);
			iss.seekg(0); // Reset buffer position to start of buffer

			// Get x position
			iss >> tmp;
			x_coord.push_back(tmp);

			// Get y position
			iss >> tmp;
			y_coord.push_back(tmp);

			// Get z position
			iss >> tmp;
			z_coord.push_back(tmp);

			// Get x velocity
			iss >> tmp;
			ux.push_back(tmp);

			// Get y velocity
			iss >> tmp;
			uy.push_back(tmp);

			// Get z velocity
			iss >> tmp;
			uz.push_back(tmp);

		}

	}

}

// ****************************************************************************
/// \brief	Tests whether a site is on a given grid.
/// \param	i	local i-index.
/// \param	j	local j-index.
/// \param	k	local k-index.
/// \param	g	grid on which to check.
/// \return boolean answer.
bool GridUtils::isOffGrid(int i, int j, int k, GridObj const * const g) {

	if (
		(i >= g->N_lim || i < 0) ||
		(j >= g->M_lim || j < 0) ||
		(k >= g->K_lim || k < 0)
		)
	{
		return true;
	}

	return false;
}

// ****************************************************************************
/// \brief	Get direction in MPI topology from unit vector.
///
/// \param	offset_vector	unit vector pointing away from current rank.
/// \return	MPI direction.
int GridUtils::getMpiDirection(int offset_vector[])
{

	// Get Mpi Manager
	MpiManager *mpim = MpiManager::getInstance();

	// Loop over the directions
	for (int d = 0; d < L_MPI_DIRS; ++d)
	{
		if (offset_vector[0] == mpim->neighbour_vectors[0][d] &&
			offset_vector[1] == mpim->neighbour_vectors[1][d]

#if (L_DIMS == 3)
			&& offset_vector[2] == mpim->neighbour_vectors[2][d]
#endif
			)
		{
			return d;
		}
	}

	return -1;

}

// ****************************************************************************
/// \brief	Get local voxel indices on grid in which provided position lies.
///
///			Wrapper for the overload which concentates all check into a vector.
///
/// \param	x	x-position.
/// \param	y	y-position.
/// \param	z	z-position.
/// \param	g	lattice on which to look for enclosing voxel.
/// \param	ijk	pointer to vector where indices are to be placed.
void GridUtils::getEnclosingVoxel(double x, double y, double z, GridObj const * const g, std::vector<int> *ijk) {

	// Declarations
	ijk->clear();
	int temp;

	// X //
	getEnclosingVoxel(x, g, eXDirection, &temp);
	ijk->push_back(temp);

	// Y //
	getEnclosingVoxel(y, g, eYDirection, &temp);
	ijk->push_back(temp);

#if (L_DIMS == 3)
	// Z //
	getEnclosingVoxel(z, g, eZDirection, &temp);
	ijk->push_back(temp);
#else
	ijk->push_back(0);
#endif

	return;

}

// ****************************************************************************
/// \brief	Get local voxel indices on grid in which provided position lies.
///
///			Will return the 1D voxel index of the voxel on the lattice 
///			provided within which point with position (xyz) lies.
///			This is done by rounding the position to obtain how many voxels
///			in from the grid core edge it is, then accounting for whether the
///			grid starts on another rank, in the halo, or further into the grid
///			by offsetting the original index by this amount. This approach saves
///			expensive seraches of the position vectors on each grid.
///			This method can be used as a position -> voxel converter.
///			The index may be off grid so it is advisable to call isOnThisRank instead.
///
/// \param	xyz	x, y or z-position.
/// \param	g	lattice on which to look for enclosing voxel.
///	\param	dir	1D direction.
/// \param	ijk	pointer to local index storage location.
void GridUtils::getEnclosingVoxel(double xyz, GridObj const * const g, eCartesianDirection dir, int *ijk) {

	// Declarations
	int offset, idxLower, idxUpper;
	double offset_baseline, local_edge;

	// Get GM instance
	GridManager *gm = GridManager::getInstance();
	
	// Find how far point is from edge of grid core and round to nearest voxel.

	if (dir == eXDirection)
	{
		idxLower = eXMin;
		idxUpper = eXMax;
	}
	else if (dir == eYDirection)
	{
		idxLower = eYMin;
		idxUpper = eYMax;
	}
	else if (dir == eZDirection)
	{
		idxLower = eZMin;
		idxUpper = eZMax;
	}

	// Set offset baseline to grid start edge of grid for serial builds
	offset_baseline = gm->global_edges[idxLower][g->level + g->region_number * L_NUM_LEVELS];
	local_edge = offset_baseline;

#ifdef L_BUILD_FOR_MPI
	MpiManager *mpim = MpiManager::getInstance();
	int rank = GridUtils::safeGetRank();

	/* If less than the lower grid core edge but is on the upper halo
	 * then must be a TL site in a periodically wrapped upper halo. */
	if (xyz < mpim->rank_core_edge[idxLower][rank] && 
		GridUtils::isOnRecvLayer(xyz, static_cast<eCartMinMax>(idxUpper)))
	{
		// Set to projected position off top of rank core to get correct ijk
		xyz = mpim->rank_core_edge[idxUpper][rank] + 
			abs(xyz - gm->global_edges[idxLower][g->level]);

		// Set offset baseline to the upper edge of rank core
		offset_baseline = mpim->rank_core_edge[idxUpper][rank];
	}

	/* If greater than the upper grid core edge but is on the lower halo
	* then must be a TL site in a periodically wrapped lower halo. */
	else if (xyz > mpim->rank_core_edge[idxUpper][rank] && 
		GridUtils::isOnRecvLayer(xyz, static_cast<eCartMinMax>(idxLower)))
	{
		// Set to projected position off bottom of rank core to get correct ijk
		xyz = mpim->rank_core_edge[idxLower][rank] - 
			abs(xyz - gm->global_edges[idxUpper][g->level]);

		// Set offset baseline to the upper edge of rank core
		offset_baseline = mpim->rank_core_edge[idxLower][rank];
	}

	// Get local lower edge
	local_edge = mpim->rank_core_edge[idxLower][rank];

#endif // L_BUILD_FOR_MPI

	// Compute number of complete cells between point and edge of rank core
	*ijk = static_cast<int>(
		std::floor((xyz - local_edge) / g->dh)
		);

	/* Compute offset from global edge i.e. number of voxels between grid edge and 
	 * rank egde. We round it as although they are exact (as they are MPIM quantities
	 * there may be round-off errors which might give an non-zero answer when zero 
	 * was expected) */
	offset = static_cast<int>(
		std::round((local_edge - offset_baseline) / g->dh)
		);

	/* If the origin is not on the halo but outside this rank to the left then offset
	 * is at most the width of the halo in voxels. This offset can be negative if the 
	 * grid starts somewhere on this rank. Does not check though to see whether it is 
	 * on the grid though and may return indices that are negative if the point is off
	 * the grid supplied. */
	if (offset > static_cast<int>(1.0 / g->refinement_ratio))
		offset = static_cast<int>(1.0 / g->refinement_ratio);

	// If on L0, always add a halo offset of 1 for MPI builds
#ifdef L_BUILD_FOR_MPI
	if (g->level == 0)
		offset = 1;
#endif

	// Correct the ijk position
	*ijk += offset;

	return;

}


// ****************************************************************************
/// \brief	Safe method to get the rank number.
///
///			This is a serial/parallel agnostic method to get the rank number.
///			This is necessary as often we just want to access the rank number for
///			logging purposes and don't want to have to wrap every call to avoid
///			attempts to use the MPI manager in non-MPI code.
///
///	\returns	integer specifying the rank number. Zero if using serial code.
int GridUtils::safeGetRank()
{
	
#ifdef L_BUILD_FOR_MPI

	// Get MPI Manager Instance
	MpiManager *mpim = MpiManager::getInstance();
	return mpim->my_rank;

#else

	return 0;

#endif

}