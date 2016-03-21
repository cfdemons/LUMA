#include "../inc/stdafx.h"
#include "../inc/BFLBody.h"
#include "GridObj.h"
#include "../inc/globalvars.h"
#include "../inc/MpiManager.h"


// Implementation of BFL body class //

/*********************************************/
// Default constructor / destructor
BFLBody::BFLBody(void)
{
}


BFLBody::~BFLBody(void)
{
}


/*********************************************/
// Custom constructor to populate body from array of points
BFLBody::BFLBody(PCpts* _PCpts, GridObj* g) {

	// Assign pointer to owning grid
	GridUtils::getGrid(g,bfl_on_grid_lev,bfl_on_grid_reg,this->_Owner);

	// Voxel grid filter //

	*GridUtils::logfile << "ObjectManagerBFL: Applying voxel grid filter..." << std::endl;

	// Place first BFL marker
	addMarker(_PCpts->x[0], _PCpts->y[0], _PCpts->z[0]);

	// Increment counters
	int curr_marker = 0;
	std::vector<int> counter;
	counter.push_back(1);

	// Loop over array of points
	for (size_t a = 1; a < _PCpts->x.size(); a++) {
    
		// Pass to point builder
		bflMarkerAdder(_PCpts->x[a], _PCpts->y[a], _PCpts->z[a], curr_marker, counter);

	}

	*GridUtils::logfile << "ObjectManagerBFL: Object represented by " << std::to_string(markers.size()) << 
		" markers using 1 marker / voxel voxelisation." << std::endl;

#ifdef BFL_DEBUG
	std::ofstream file;
	file.open(GridUtils::path_str + "/marker_data_rank" + std::to_string(MpiManager::my_rank) + ".out",std::ios::out);
	for (size_t n = 0; n < markers.size(); n++) {
		file << std::to_string(n) << ", " << 
			markers[n].position[0] << ", " << markers[n].position[1] << ", " << markers[n].position[2] << ", " <<
			markers[n].supp_i[0] << ", " << markers[n].supp_j[0] << ", " << markers[n].supp_k[0] << std::endl;
	}
	file.close();
#endif



	// Labelling //
	*GridUtils::logfile << "ObjectManagerBFL: Labelling lattice voxels..." << std::endl;

	int N_lim = _Owner->XInd.size();
	int M_lim = _Owner->YInd.size();
	int K_lim = _Owner->ZInd.size();

	// Get with times and start using for_each instead of old-fashioned C syntax...
	for (Marker& m : markers) {

		// When using MPI need to convert the global indices of the support sites to local indices for array access
#ifdef BUILD_FOR_MPI
		std::vector<int> locals; GridUtils::global_to_local(m.supp_i[0],m.supp_j[0],m.supp_k[0],_Owner,locals);
		_Owner->LatTyp(locals[0],locals[1],locals[2],M_lim,K_lim) = 10;
#else
		_Owner->LatTyp(m.supp_i[0],m.supp_j[0],m.supp_k[0],M_lim,K_lim) = 10;
#endif

	}


	// Compute Q //
	*GridUtils::logfile << "ObjectManagerBFL: Computing Q..." << std::endl;

	// Initialise Q stores as max double value
	Q = std::vector< std::vector<double> > (nVels * 2, std::vector<double> ( 
		markers.size(), std::numeric_limits<double>::max() ) );

	// Loop over local grid and inspect the streaming operations
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {

				// If site is a BFL voxel
				if (_Owner->LatTyp(i,j,k,M_lim,K_lim) == 10) {

					// Compute Q for all stream vectors storing on source voxel BFL marker
#if (dims == 3)
					computeQ(i,j,k,N_lim,M_lim,K_lim,_Owner);
#else
					computeQ(i,j,N_lim,M_lim,_Owner);
#endif

				}


			}
		}
	}

	// Computation of Q complete
	*GridUtils::logfile << "ObjectManagerBFL: Q computation complete." << std::endl;

#ifdef BFL_DEBUG
	file.open(GridUtils::path_str + "/marker_Qs_rank" + std::to_string(MpiManager::my_rank) + ".out",std::ios::out);
	for (std::vector<double> i : Q) {
		for (size_t n = 0; n < markers.size(); n++) {
			file << i[n] << '\t';
		}
		file << std::endl;
	}
	file.close();
#endif

}

/*********************************************/
// Routine to add markers to the current BFL body
void BFLBody::bflMarkerAdder(double x, double y, double z, int& curr_mark, std::vector<int>& counter) {

	// If point in current BFL voxel
	if (isInVoxel(x,y,z,curr_mark,this)) {

		// Increment point counter
		counter[curr_mark]++;

		// Update position of marker in current BFL voxel
		markers[curr_mark].position[0] = 
			( (markers[curr_mark].position[0] * (counter[curr_mark] - 1)) + x) / counter[curr_mark];
		markers[curr_mark].position[1] = 
			( (markers[curr_mark].position[1] * (counter[curr_mark] - 1)) + y) / counter[curr_mark];
		markers[curr_mark].position[2] = 
			( (markers[curr_mark].position[2] * (counter[curr_mark] - 1)) + z) / counter[curr_mark];


	// If point is in an existing BFL voxel
	} else if (isVoxelBflVoxel(x,y,z,this)) {

		// Recover voxel number
		MarkerData m_data = getMarkerData(x,y,z,this);
		curr_mark = m_data.ID;

		// Increment point counter
		counter[curr_mark]++;

		// Update position of marker in current BFL voxel
		markers[curr_mark].position[0] = 
			( (markers[curr_mark].position[0] * (counter[curr_mark] - 1)) + x) / counter[curr_mark];
		markers[curr_mark].position[1] = 
			( (markers[curr_mark].position[1] * (counter[curr_mark] - 1)) + y) / counter[curr_mark];
		markers[curr_mark].position[2] = 
			( (markers[curr_mark].position[2] * (counter[curr_mark] - 1)) + z) / counter[curr_mark];


	// Must be in a new BFL voxel
	} else {

		// Reset counter and increment voxel index
		curr_mark = counter.size();
		counter.push_back(1);        

		// Create new marker as this is a new BFL voxel
		addMarker(x,y,z);

	}


}

/*********************************************/
// Routine to compute Q values at a given local voxel for each application of the BFL BC.
void BFLBody::computeQ(int i, int j, int k, int N_lim, int M_lim, int K_lim, GridObj* g) {

	// Declarations
	int dest_i, dest_j, dest_k, ib, jb, kb, storeID;
	MarkerData m_data;

	/* Get voxel IDs of self and stencil required to specify planes
	 *
	 * 7 -- 16 -- 25       8 -- 17 -- 26		9 -- 18 -- 27    
	 * |    |     |        |     |     |		|     |     |
	 * 4 -- 13 -- 22       5 -- 14 -- 23		6 -- 15 -- 24
	 * |    |     |        |     |     |		|     |     |
	 * 1 -- 10 -- 19       2 -- 11 -- 20		3 -- 12 -- 21
	 *
	 *   Front                 Middle               Back
	 */


	// Set stencil centre data
	ib = i;
	jb = j;
	kb = k;

	// Get marker data associated with this local site
#ifdef BUILD_FOR_MPI
	// Convert to global indices for marker access
	std::vector<int> globals; GridUtils::local_to_global(i,j,k,g,globals);
	m_data = getMarkerData(globals[0],globals[1],globals[2],this);
#else
	m_data = getMarkerData(i,j,k,this);
#endif

	storeID = m_data.ID;

	// Get list of IDs of neighbour vertices for plane construction
	std::vector<int> V;
	for (int ii = ib - 1; ii <= ib + 1; ii++) {
		for (int jj = jb - 1; jj <= jb + 1; jj++) {
			for (int kk = kb - 1; kk <= kb + 1; kk++) {

				// If indices are valid
				if (ib >= 0 && ib < N_lim && jb >= 0 && jb < M_lim && kb >= 0 && kb < K_lim) {

					// Fetch data if available //
#ifdef BUILD_FOR_MPI
					// Convert to global indices for marker access
					globals.clear(); GridUtils::local_to_global(ii,jj,kk,g,globals);
					m_data = getMarkerData(globals[0],globals[1],globals[2],this);
#else
					m_data = getMarkerData(ii,jj,kk,this);
#endif

					// If data valid, then store ID
					if (!_isnan(m_data.x)) V.push_back(m_data.ID);
				}

			}
		}
	}


	// Build a triangular plane for each combination of vertices //

	// Cannot compute Q if not enough neighbour markers to make a triangle
	if (V.size() < 3) return;

	// Get unique combinations of vertex IDs
	std::vector< std::vector<int> > combo;

	// Sort IDs in ascending order
	std::sort(V.begin(), V.end());

	// Loop over sorted array and store ID combinations
	std::vector<int> tmp;
	std::vector<bool> flags(V.size());
	std::fill(flags.begin(), flags.end() - V.size() + 3, true);

	// Loop over permutations
	do {

		// Reset tmp array
		tmp.clear();

		// Loop over permutation and pull corresponding 
		// values from V into combo array
		for (size_t i = 0; i < V.size(); ++i) {
			if (flags[i]) {
				tmp.push_back(V[i]);
			}
		}
		if (tmp.size() == 3) combo.push_back(tmp);

	} while (std::prev_permutation(flags.begin(), flags.end()));

	// Loop over each triangle
	for (std::vector<int>& tri : combo) {

		// Consider only triangles where indices are single lattice site apart
		if (
			( GridUtils::vecnorm( std::abs(markers[tri[1]].supp_i[0] - markers[tri[0]].supp_i[0]),
			std::abs(markers[tri[1]].supp_j[0] - markers[tri[0]].supp_j[0]), 
			std::abs(markers[tri[1]].supp_k[0] - markers[tri[0]].supp_k[0]) ) < sqrt(2) )
			&&
			( GridUtils::vecnorm( std::abs(markers[tri[2]].supp_i[0] - markers[tri[1]].supp_i[0]),
			std::abs(markers[tri[2]].supp_j[0] - markers[tri[1]].supp_j[0]), 
			std::abs(markers[tri[2]].supp_k[0] - markers[tri[1]].supp_k[0]) ) < sqrt(2) )
			&&
			( GridUtils::vecnorm( std::abs(markers[tri[2]].supp_i[0] - markers[tri[0]].supp_i[0]),
			std::abs(markers[tri[2]].supp_j[0] - markers[tri[0]].supp_j[0]), 
			std::abs(markers[tri[2]].supp_k[0] - markers[tri[0]].supp_k[0]) ) < sqrt(2) )
			) {

			// Perform 3D line-triangle intersection test to get Q //

			// Loop over even velocities and ignore rest distribution to save computing Q twice
			for (int vel = 0; vel < nVels - 1; vel+=2) {

				// Compute destination coordinates
				dest_i = (i + c[0][vel] + N_lim) % N_lim;
				dest_j = (j + c[1][vel] + M_lim) % M_lim;
				dest_k = (k + c[2][vel] + K_lim) % K_lim;

				// Define vectors
				std::vector<double> u, v, local_origin;
				u.push_back(markers[tri[1]].position[0] - markers[tri[0]].position[0]);
				u.push_back(markers[tri[1]].position[1] - markers[tri[0]].position[1]);
				u.push_back(markers[tri[1]].position[2] - markers[tri[0]].position[2]);
				v.push_back(markers[tri[2]].position[0] - markers[tri[0]].position[0]);
				v.push_back(markers[tri[2]].position[1] - markers[tri[0]].position[1]);
				v.push_back(markers[tri[2]].position[2] - markers[tri[0]].position[2]);
				local_origin.push_back( markers[tri[0]].position[0] );
				local_origin.push_back( markers[tri[0]].position[1] );
				local_origin.push_back( markers[tri[0]].position[2] );

				// Cross product gives normal vector to plane
				std::vector<double> n = GridUtils::crossprod(u,v);

				if (GridUtils::vecnorm(n) == 0) continue; // Triangle degenerate

				// Create global vectors to source and destination end of streaming vector
				std::vector<double> src; 
				GridUtils::local_to_global(i,j,k,g,src);
				std::vector<double> dest; 
				GridUtils::local_to_global(dest_i,dest_j,dest_k,g,dest);

				std::vector<double> dir = GridUtils::subtract(dest, src);
				std::vector<double> w0 = GridUtils::subtract(src,local_origin);
				double a = -GridUtils::dotprod(n,w0);
				double b = GridUtils::dotprod(n,dir);

				if (abs(b) < 1e-9) {
					if (a == 0) continue;	// Triangle and line are in the same plane
					else continue;			// Triangle and line are disjoint
				}


				// Get intersect point
				double r = a / b;

				if (r < 0 || r > 1) continue; // No intersect

				std::vector<double> intersect = GridUtils::add(src, GridUtils::vecmultiply(r,dir) );    // Intersect point

				double uu = GridUtils::dotprod(u,u);
				double uv = GridUtils::dotprod(u,v);
				double vv = GridUtils::dotprod(v,v);
				std::vector<double> w = GridUtils::subtract(intersect, local_origin);
				double wu = GridUtils::dotprod(w,u);
				double wv = GridUtils::dotprod(w,v);
				double D = uv * uv - uu * vv;

				double s = (uv * wv - vv * wu) / D;
				double t = (uv * wu - uu * wv) / D;

				if (s < 0.0 || s > 1.0)	continue;				// I is outside T
				else if (t < 0.0 || (s + t) > 1.0) continue;	// I is outside T
				else {
					// Inside T so compute Q
					double q = GridUtils::vecnorm( GridUtils::subtract(intersect,src) ) / GridUtils::vecnorm(dir);
					if (q < Q[vel][storeID]) {

						// Set outgoing Q value
						Q[vel][storeID] = q;

						// Incoming Q value (in destination store at opposite direction) is 1 minus the incoming
						Q[GridUtils::getOpposite(vel) + nVels][storeID] = 1 - q;
					}
				}

			}

		}
	}

}

// Overloaded compute Q routine for 2D
void BFLBody::computeQ(int i, int j, int N_lim, int M_lim, GridObj* g) {

	/* For each possible line extending from the marker voxel to an 
	 * immediate neighbour we can check for an intersection between the 
	 * stream vector and the constructed line. The distance of intersection
	 * is the value of q. As with the 3D case, 1-q is assigned to the 
	 * destination store after we have added q to the source store. */

	// Declarations
	int dest_i, dest_j;
	MarkerData m_data;

	// Set stencil centre data
	int ib = i;
	int jb = j;
	
	// Get marker data associated with this local site
#ifdef BUILD_FOR_MPI
	// Convert to global indices for marker access
	std::vector<int> globals; GridUtils::local_to_global(i,j,0,g,globals);
	m_data = getMarkerData(globals[0],globals[1],0,this);
#else
	m_data = getMarkerData(i,j,0,this);
#endif

	int storeID = m_data.ID;


	// Get list of IDs of neighbour vertices for line construction
	std::vector<int> V;
	for (int ii = ib - 1; ii <= ib + 1; ii++) {
		for (int jj = jb - 1; jj <= jb + 1; jj++) {

			// If indices are valid
			if (ib >= 0 && ib < N_lim && jb >= 0 && jb < M_lim) {
			

				// Fetch data if available //
#ifdef BUILD_FOR_MPI
				// Convert to global indices for marker access
				globals.clear(); GridUtils::local_to_global(ii,jj,0,g,globals);
				m_data = getMarkerData(globals[0],globals[1],0,this);
#else
				m_data = getMarkerData(ii,jj,0,this);
#endif

				// If data valid, then store ID
				if (!_isnan(m_data.x)) V.push_back(m_data.ID);

			}

		}
	}

	// Can only continue if we have at least 2 points
	if (V.size() < 2) return;

	// Get unique combinations of vertex IDs
	std::vector< std::vector<int> > combo;

	// Sort IDs in ascending order
	std::sort(V.begin(), V.end());

	// Loop over sorted array and store ID combinations
	std::vector<int> tmp;
	std::vector<bool> flags(V.size());
	std::fill(flags.begin(), flags.end() - V.size() + 2, true);

	// Loop over permutations
	do {

		// Reset tmp array
		tmp.clear();

		// Loop over permutation and pull corresponding 
		// values from V into combo array
		for (size_t i = 0; i < V.size(); ++i) {
			if (flags[i]) {
				tmp.push_back(V[i]);
			}
		}
		if (tmp.size() == 2) combo.push_back(tmp);

	} while (std::prev_permutation(flags.begin(), flags.end()));



	// Loop through valid marker combinations
	for (std::vector<int>& line : combo) {

		// Only check intersections when the surface line is between immediate neighbours
		if ( GridUtils::vecnorm(
				std::abs(markers[line[1]].supp_i[0] - markers[line[0]].supp_i[0]),
				std::abs(markers[line[1]].supp_j[0] - markers[line[0]].supp_j[0]), 
				std::abs(markers[line[1]].supp_k[0] - markers[line[0]].supp_k[0]) 
				) < sqrt(2)			
			) {


			// Perform line intersection test //

			// Loop over velocities (ignore rest distribution)
			for (int vel = 0; vel < nVels - 1; vel++) {

				// Compute destination coordinates
				dest_i = (i + c[0][vel] + N_lim) % N_lim;
				dest_j = (j + c[1][vel] + M_lim) % M_lim;

				// Compute vectors
				std::vector<double> p;
				GridUtils::local_to_global(i,j,0,g,p);	// Position of source site
				std::vector<double> ppr;
				GridUtils::local_to_global(dest_i,dest_j,0,g,ppr);	// Position of destination site

				std::vector<double> q;	// Position of first marker
				q.push_back(markers[line[0]].position[0]);
				q.push_back(markers[line[0]].position[1]);
				q.push_back(markers[line[0]].position[2]);

				std::vector<double> r = GridUtils::subtract(ppr,p);	// Length of streaming vector

				std::vector<double> s;	// Length of marker vector
				s.push_back(markers[line[1]].position[0] - q[0]);
				s.push_back(markers[line[1]].position[1] - q[1]);
				s.push_back(markers[line[1]].position[2] - q[2]);
    
				// Compute parameters
				double rXs = GridUtils::vecnorm(GridUtils::crossprod( r, s ));
				double qpXr = GridUtils::vecnorm(GridUtils::crossprod( GridUtils::subtract(q, p), r ));
				double qpXs = GridUtils::vecnorm(GridUtils::crossprod( GridUtils::subtract(q, p), s ));
				double t = qpXs / rXs;
				double u = qpXr / rXs;    
    
				// Check for the intersecting case
				if (rXs != 0 && t >= 0 && t <= 1 && u >= 0 && u <= 1) {
        
					// Lines intersect at p + tr = q + us
					std::vector<double> intersect;
					intersect.push_back(p[0] + t * r[0]);
					intersect.push_back(p[1] + t * r[1]);
					intersect.push_back(0);
        
					// Compute Q
					double q = GridUtils::vecnorm( GridUtils::subtract(intersect, p) ) / GridUtils::vecnorm(r);
					if (q < Q[vel][storeID]) {

						// Set outgoing Q value
						Q[vel][storeID] = q;

						// Incoming Q value (in destination store at opposite direction) is 1 minus the incoming
						Q[GridUtils::getOpposite(vel) + nVels][storeID] = 1 - q;
					}
        
				}

				// Maybe include the collinear case if it cuts through a vertex?

			}

		}


	}



}


/*********************************************/
// Utilities//

// Return global voxel indices for a given point in global space
std::vector<int> BFLBody::getVoxInd(double x, double y, double z) {

	std::vector<int> vox;

	if (x - (int)std::floor(x) > 0.5) vox.push_back( (int)std::ceil(x) );
	else vox.push_back( (int)std::floor(x) );

	if (y - (int)std::floor(y) > 0.5) vox.push_back( (int)std::ceil(y) );
	else vox.push_back( (int)std::floor(y) );

	if (z - (int)std::floor(z) > 0.5) vox.push_back( (int)std::ceil(z) );
	else vox.push_back( (int)std::floor(z) );

	return vox;

}

// Overload of above for a single index
int BFLBody::getVoxInd(double p) {

	if (p - (int)std::floor(p) > 0.5) return (int)std::ceil(p);
	else return (int)std::floor(p);

}

// Return marker and voxel data associated with global position supplied by querying given body*
MarkerData BFLBody::getMarkerData(double x, double y, double z, BFLBody* body) {

	// Get indices of voxel associated with the supplied position
	std::vector<int> vox = getVoxInd(x,y,z);

    // Find all marker IDs that match these indices
	std::vector<int> found_i, found_j, found_k;
	for (size_t i = 0; i < body->markers.size(); ++i) {
		if (body->markers[i].supp_i[0] == vox[0]) found_i.push_back(i);
		if (body->markers[i].supp_j[0] == vox[1]) found_j.push_back(i);
		if (body->markers[i].supp_k[0] == vox[2]) found_k.push_back(i);
	}

    // Loop over marker IDs and check for matches in the other two found vectors
    for (size_t ind_x = 0; ind_x < found_i.size(); ++ind_x) {

		int matching_index = found_i[ind_x];

		// If a match exists in both other vectors (lambda expressions for C++)
        if (	std::any_of(found_j.begin(), found_j.end(), [matching_index](int index) {
					return index == matching_index; }) && 
				std::any_of(found_k.begin(), found_k.end(), [matching_index](int index) { 
					return index == matching_index; }) ) {

			// "matching_index" is a valid ID so create new MarkerData store
			MarkerData* m_MarkerData = new MarkerData(	body->markers[matching_index].supp_i[0],
														body->markers[matching_index].supp_j[0],
														body->markers[matching_index].supp_k[0],
														body->markers[matching_index].position[0],
														body->markers[matching_index].position[1],
														body->markers[matching_index].position[2],
														matching_index
														);
			return *m_MarkerData;	// Return the store information through a copy
		}
	
	}

	// Create empty MarkerData store using default constructor
	MarkerData* m_MarkerData = new MarkerData();
	return *m_MarkerData;

}

// Returns boolean as to whether a given point is in the voxel associated with <curr_mark>
bool BFLBody::isInVoxel(double x, double y, double z, int curr_mark, BFLBody* body) {

	try {

		// Try to retrieve the position of the marker identified by <curr_mark>
		int i = body->markers[curr_mark].supp_i[0];
		int j = body->markers[curr_mark].supp_j[0];
		int k = body->markers[curr_mark].supp_k[0];

		// Test within
		if (	(x > i - 0.5 && x <= i + 0.5) && 
				(y > j - 0.5 && y <= j + 0.5) && 
				(z > k - 0.5 && z <= k + 0.5)
			) return true;

	// Catch all
	} catch (...) {

		// If failed, marker probably doesn't exist so return false
		return false;

	}

	return false;

}

// Returns boolean as to whether a given point is in an existing Bfl voxel
bool BFLBody::isVoxelBflVoxel(double x, double y, double z, BFLBody* body) {

	// Try get the MarkerData store
	MarkerData m_data = getMarkerData(x,y,z,body);

	// True if the data store is not empty
	if (!_isnan(m_data.x)) {

		return true;

	}

	return false;

}
