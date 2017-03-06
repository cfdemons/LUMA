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
#include "../inc/BFLBody.h"
#include "../inc/PCpts.h"
#include "../inc/GridObj.h"


// Implementation of BFL body class //

/*********************************************/
/// Default constructor
BFLBody::BFLBody(void)
{
}

/// Default destructor
BFLBody::~BFLBody(void)
{
}


/*********************************************/
/// \brief Custom constructor to populate body from array of points.
/// \param g		hierarchy pointer to grid hierarchy
/// \param bodyID	ID of body in array of bodies.
/// \param _PCpts	pointer to point cloud data
BFLBody::BFLBody(GridObj* g, int bodyID, PCpts* _PCpts) 
	: Body(g, bodyID, _PCpts)
{

	// Initialiser list ensures the correct super class constructor is called first

#ifdef L_BFL_DEBUG
	std::ofstream file;
	file.open(GridUtils::path_str + "/marker_data_rank" + std::to_string(rank) + ".out",std::ios::out);
	file.precision(L_OUTPUT_PRECISION);
	for (size_t n = 0; n < markers.size(); n++) {
		file << std::to_string(n) << ", " << 
			markers[n].position[0] << ", " << markers[n].position[1] << ", " << markers[n].position[2] << ", " <<
			markers[n].supp_i[0] << ", " << markers[n].supp_j[0] << ", " << markers[n].supp_k[0] << std::endl;
}
	file.close();
#endif

	// Labelling //
	*GridUtils::logfile << "ObjectManagerBFL: Labelling lattice voxels..." << std::endl;

	int N_lim = _Owner->N_lim;
	int M_lim = _Owner->M_lim;
	int K_lim = _Owner->K_lim;

	// Get each marker in turn
	for (Marker& m : markers) {

		// Label as BFL site
		_Owner->LatTyp(m.supp_i[0],m.supp_j[0],m.supp_k[0],M_lim,K_lim) = eBFL;
	}


	// Compute Q //
	*GridUtils::logfile << "ObjectManagerBFL: Computing Q..." << std::endl;

	// Initialise Q stores to the "invalid" value
	Q = std::vector< std::vector<double> > (L_NUM_VELS * 2, std::vector<double> ( 
		markers.size(), -1.0 ) );

	// Loop over local grid and inspect the streaming operations
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {

				// If site is a BFL voxel
				if (_Owner->LatTyp(i,j,k,M_lim,K_lim) == eBFL) {

					// Compute Q for all stream vectors storing on source voxel BFL marker
#if (L_DIMS == 3)
					computeQ(i,j,k,_Owner);
#else
					computeQ(i,j,_Owner);
#endif

				}
			}
		}
	}

	// Computation of Q complete
	*GridUtils::logfile << "ObjectManagerBFL: Q computation complete." << std::endl;

#ifdef L_BFL_DEBUG
	file.open(GridUtils::path_str + "/marker_Qs_rank" + std::to_string(rank) + ".out",std::ios::out);
	file.precision(L_OUTPUT_PRECISION);
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
/// \brief 	Custom constructor for building prefab filament
/// \param g				hierarchy pointer to grid hierarchy
/// \param bodyID			ID of body in array of bodies.
/// \param start_position	start position of base of filament
/// \param length			length of filament
/// \param angles			angle of filament
BFLBody::BFLBody(GridObj* g, int bodyID, std::vector<double> &start_position,
		double length, std::vector<double> &angles) : Body(g, bodyID, start_position, length, angles)
{

#ifdef L_BFL_DEBUG
	std::ofstream file;
	file.open(GridUtils::path_str + "/marker_data_rank" + std::to_string(rank) + ".out",std::ios::out);
	file.precision(L_OUTPUT_PRECISION);
	for (size_t n = 0; n < markers.size(); n++) {
		file << std::to_string(n) << ", " <<
			markers[n].position[0] << ", " << markers[n].position[1] << ", " << markers[n].position[2] << ", " <<
			markers[n].supp_i[0] << ", " << markers[n].supp_j[0] << ", " << markers[n].supp_k[0] << std::endl;
}
	file.close();
#endif

	// Labelling //
	*GridUtils::logfile << "ObjectManagerBFL: Labelling lattice voxels..." << std::endl;

	int N_lim = _Owner->N_lim;
	int M_lim = _Owner->M_lim;
	int K_lim = _Owner->K_lim;

	// Get each marker in turn
	for (Marker& m : markers) {

		// Label as BFL site
		_Owner->LatTyp(m.supp_i[0],m.supp_j[0],m.supp_k[0],M_lim,K_lim) = eBFL;
	}


	// Compute Q //
	*GridUtils::logfile << "ObjectManagerBFL: Computing Q..." << std::endl;

	// Initialise Q stores to the "invalid" value
	Q = std::vector< std::vector<double> > (L_NUM_VELS * 2, std::vector<double> (
		markers.size(), -1.0 ) );

	// Loop over local grid and inspect the streaming operations
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {

				// If site is a BFL voxel
				if (_Owner->LatTyp(i,j,k,M_lim,K_lim) == eBFL) {

					// Compute Q for all stream vectors storing on source voxel BFL marker
#if (L_DIMS == 3)
					computeQ(i,j,k,_Owner);
#else
					computeQ(i,j,_Owner);
#endif

				}
			}
		}
	}

	// Computation of Q complete
	*GridUtils::logfile << "ObjectManagerBFL: Q computation complete." << std::endl;

#ifdef L_BFL_DEBUG
	file.open(GridUtils::path_str + "/marker_Qs_rank" + std::to_string(rank) + ".out",std::ios::out);
	file.precision(L_OUTPUT_PRECISION);
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
/// \brief 	Custom constructor for building prefab circle/sphere
/// \param g				hierarchy pointer to grid hierarchy
/// \param bodyID			ID of body in array of bodies.
/// \param centre_point		centre point of circle
/// \param radius			radius of circle
BFLBody::BFLBody(GridObj* g, int bodyID, std::vector<double> &centre_point,
		double radius) : Body(g, bodyID, centre_point, radius)
{

#ifdef L_BFL_DEBUG
	std::ofstream file;
	file.open(GridUtils::path_str + "/marker_data_rank" + std::to_string(rank) + ".out",std::ios::out);
	file.precision(L_OUTPUT_PRECISION);
	for (size_t n = 0; n < markers.size(); n++) {
		file << std::to_string(n) << ", " <<
			markers[n].position[0] << ", " << markers[n].position[1] << ", " << markers[n].position[2] << ", " <<
			markers[n].supp_i[0] << ", " << markers[n].supp_j[0] << ", " << markers[n].supp_k[0] << std::endl;
}
	file.close();
#endif

	// Labelling //
	*GridUtils::logfile << "ObjectManagerBFL: Labelling lattice voxels..." << std::endl;

	int N_lim = _Owner->N_lim;
	int M_lim = _Owner->M_lim;
	int K_lim = _Owner->K_lim;

	// Get each marker in turn
	for (Marker& m : markers) {

		// Label as BFL site
		_Owner->LatTyp(m.supp_i[0],m.supp_j[0],m.supp_k[0],M_lim,K_lim) = eBFL;
	}


	// Compute Q //
	*GridUtils::logfile << "ObjectManagerBFL: Computing Q..." << std::endl;

	// Initialise Q stores to the "invalid" value
	Q = std::vector< std::vector<double> > (L_NUM_VELS * 2, std::vector<double> (
		markers.size(), -1.0 ) );

	// Loop over local grid and inspect the streaming operations
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {

				// If site is a BFL voxel
				if (_Owner->LatTyp(i,j,k,M_lim,K_lim) == eBFL) {

					// Compute Q for all stream vectors storing on source voxel BFL marker
#if (L_DIMS == 3)
					computeQ(i,j,k,_Owner);
#else
					computeQ(i,j,_Owner);
#endif

				}
			}
		}
	}

	// Computation of Q complete
	*GridUtils::logfile << "ObjectManagerBFL: Q computation complete." << std::endl;

#ifdef L_BFL_DEBUG
	file.open(GridUtils::path_str + "/marker_Qs_rank" + std::to_string(rank) + ".out",std::ios::out);
	file.precision(L_OUTPUT_PRECISION);
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
/// \brief 	Custom constructor for building square/cuboid
/// \param g				hierarchy pointer to grid hierarchy
/// \param bodyID			ID of body in array of bodies.
/// \param centre_point		centre point of square
/// \param dimensions		dimensions of square
/// \param angles			angle of square
BFLBody::BFLBody(GridObj* g, int bodyID, std::vector<double> &centre_point,
		std::vector<double> &dimensions, std::vector<double> &angles) : Body(g, bodyID, centre_point, dimensions, angles)
{

#ifdef L_BFL_DEBUG
	std::ofstream file;
	file.open(GridUtils::path_str + "/marker_data_rank" + std::to_string(rank) + ".out",std::ios::out);
	file.precision(L_OUTPUT_PRECISION);
	for (size_t n = 0; n < markers.size(); n++) {
		file << std::to_string(n) << ", " <<
			markers[n].position[0] << ", " << markers[n].position[1] << ", " << markers[n].position[2] << ", " <<
			markers[n].supp_i[0] << ", " << markers[n].supp_j[0] << ", " << markers[n].supp_k[0] << std::endl;
}
	file.close();
#endif

	// Labelling //
	*GridUtils::logfile << "ObjectManagerBFL: Labelling lattice voxels..." << std::endl;

	int N_lim = _Owner->N_lim;
	int M_lim = _Owner->M_lim;
	int K_lim = _Owner->K_lim;

	// Get each marker in turn
	for (Marker& m : markers) {

		// Label as BFL site
		_Owner->LatTyp(m.supp_i[0],m.supp_j[0],m.supp_k[0],M_lim,K_lim) = eBFL;
	}


	// Compute Q //
	*GridUtils::logfile << "ObjectManagerBFL: Computing Q..." << std::endl;

	// Initialise Q stores to the "invalid" value
	Q = std::vector< std::vector<double> > (L_NUM_VELS * 2, std::vector<double> (
		markers.size(), -1.0 ) );

	// Loop over local grid and inspect the streaming operations
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {

				// If site is a BFL voxel
				if (_Owner->LatTyp(i,j,k,M_lim,K_lim) == eBFL) {

					// Compute Q for all stream vectors storing on source voxel BFL marker
#if (L_DIMS == 3)
					computeQ(i,j,k,_Owner);
#else
					computeQ(i,j,_Owner);
#endif

				}
			}
		}
	}

	// Computation of Q complete
	*GridUtils::logfile << "ObjectManagerBFL: Q computation complete." << std::endl;

#ifdef L_BFL_DEBUG
	file.open(GridUtils::path_str + "/marker_Qs_rank" + std::to_string(rank) + ".out",std::ios::out);
	file.precision(L_OUTPUT_PRECISION);
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
/// \brief	Routine to compute wall distance Q.
///
///			Computes Q values in 3D at a given local voxel for each application of 
///			the BFL BC. Performs a line-plane intersection algorithm for every 
///			possible triangular plane constructed out of the marker in the voxel
///			and its nearest neighbours.
///
/// \param i local i-index of BFL voxel
/// \param j local j-index of BFL voxel
/// \param k local k-index of BFL voxel
/// \param g pointer to owner grid
void BFLBody::computeQ(int i, int j, int k, GridObj* g) {

	// Declarations
	int dest_i, dest_j, dest_k, ib, jb, kb, storeID;
	MarkerData* m_data;

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
	m_data = getMarkerData(g->XPos[i], g->YPos[j], g->ZPos[k]);

	storeID = m_data->ID;
	delete m_data;

	// Get list of IDs of neighbour vertices for plane construction
	std::vector<int> V;
	for (int ii = ib - 1; ii <= ib + 1; ii++) {
		for (int jj = jb - 1; jj <= jb + 1; jj++) {
			for (int kk = kb - 1; kk <= kb + 1; kk++) {

				// If indices are valid and not equal to centre voxel

				if (	ib >= 0 && ib < g->N_lim
					&&	jb >= 0 && jb < g->M_lim
					&&	kb >= 0 && kb < g->K_lim
					) {

					// Fetch data if available
					m_data = getMarkerData(g->XPos[ii], g->YPos[jj], g->ZPos[kk]);

					// If data valid, then store ID
					if (m_data->ID != -1) V.push_back(m_data->ID);
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

		// Perform 3D line-triangle intersection test to get Q //

		// Define vectors for triangle vertices
		std::vector<double> u, v, local_origin;
		u.push_back(markers[tri[1]].position[0] - markers[tri[0]].position[0]);
		u.push_back(markers[tri[1]].position[1] - markers[tri[0]].position[1]);
		u.push_back(markers[tri[1]].position[2] - markers[tri[0]].position[2]);
		v.push_back(markers[tri[2]].position[0] - markers[tri[0]].position[0]);
		v.push_back(markers[tri[2]].position[1] - markers[tri[0]].position[1]);
		v.push_back(markers[tri[2]].position[2] - markers[tri[0]].position[2]);
		local_origin.push_back(markers[tri[0]].position[0]);
		local_origin.push_back(markers[tri[0]].position[1]);
		local_origin.push_back(markers[tri[0]].position[2]);

		// Create global position of start of streaming vector
		std::vector<double> src;
		src.push_back(_Owner->XPos[i]);
		src.push_back(_Owner->YPos[j]);
		src.push_back(_Owner->ZPos[k]);

		// Loop over even velocities and ignore rest distribution to save computing Q twice
		for (int vel = 0; vel < L_NUM_VELS - 1; vel+=2) {

			// Compute destination coordinates
			dest_i = (i + c[0][vel] + g->N_lim) % g->N_lim;
			dest_j = (j + c[1][vel] + g->M_lim) % g->M_lim;
			dest_k = (k + c[2][vel] + g->K_lim) % g->K_lim;

			// Cross product gives normal vector to plane
			std::vector<double> n = GridUtils::crossprod(u,v);

			if (GridUtils::vecnorm(n) == 0) continue; // Triangle degenerate

			// Global position of end of streaming vector
			std::vector<double> dest;
			dest.push_back(_Owner->XPos[dest_i]);
			dest.push_back(_Owner->YPos[dest_j]);
			dest.push_back(_Owner->ZPos[dest_k]);

			std::vector<double> dir = GridUtils::subtract(dest, src);
			std::vector<double> w0 = GridUtils::subtract(src,local_origin);
			double a = -GridUtils::dotprod(n,w0);
			double b = GridUtils::dotprod(n,dir);

			if (abs(b) < L_SMALL_NUMBER) {
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

			if (s < 0.0 || s > 1.0)	continue;			
			else if (t < 0.0 || (s + t) > 1.0) continue;
			else {
				// Inside so compute Q
				double q = GridUtils::vecnorm( GridUtils::subtract(intersect,src) ) / GridUtils::vecnorm(dir);

				// On first pass, set to valid value
				if (Q[vel][storeID] == -1) Q[vel][storeID] = std::numeric_limits<double>::max();

				if (q < Q[vel][storeID]) {

					// Set outgoing Q value
					Q[vel][storeID] = q;

					// Incoming Q value (in destination store at opposite direction) is 1 minus the incoming
					Q[GridUtils::getOpposite(vel) + L_NUM_VELS][storeID] = 1 - q;
				}
			}
		}
	}

}

/// \brief	Routine to compute wall distance Q.
///
///			Computes Q values in 2D at a given local voxel for each application of 
///			the BFL BC. Performs a line-line intersection algorithm for each line 
///			segment either side of the voxel marker.
///
/// \param i local i-index of BFL voxel
/// \param j local j-index of BFL voxel
/// \param g pointer to owner grid
void BFLBody::computeQ(int i, int j, GridObj* g) {

	/* For each possible line extending from the marker voxel to an 
	 * immediate neighbour we can check for an intersection between the 
	 * stream vector and the constructed line. The distance of intersection
	 * is the value of q. As with the 3D case, 1-q is assigned to the 
	 * destination store after we have added q to the source store. */

	// Declarations
	int dest_i, dest_j;
	MarkerData* m_data;

	// Set stencil centre data
	int ib = i;
	int jb = j;
	
	// Get marker data associated with this local site
	m_data = getMarkerData(g->XPos[i], g->YPos[j], g->ZPos[0]);

	int storeID = m_data->ID;
	delete m_data;


	// Get list of IDs of neighbour vertices for line construction
	std::vector<int> V;
	for (int ii = ib - 1; ii <= ib + 1; ii++) {
		for (int jj = jb - 1; jj <= jb + 1; jj++) {

			// Check local indices on the grid
			if (	ib >= 0 && ib < g->N_lim
				&&	jb >= 0 && jb < g->M_lim
				) {			

				// Fetch data if available
				m_data = getMarkerData(g->XPos[ii], g->YPos[jj], g->ZPos[0]);

				// If data valid, then store ID
				if (m_data->ID != -1) V.push_back(m_data->ID);

				// Clean-up after getMarkerData call
				delete m_data;

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

		// Perform line intersection test //

		std::vector<double> q;	// Position of adjacent marker
		q.push_back(markers[line[0]].position[0]);
		q.push_back(markers[line[0]].position[1]);
		q.push_back(markers[line[0]].position[2]);

		std::vector<double> s;	// Length of marker-marker vector
		s.push_back(markers[line[1]].position[0] - q[0]);
		s.push_back(markers[line[1]].position[1] - q[1]);
		s.push_back(markers[line[1]].position[2] - q[2]);

		std::vector<double> p;	// Position of source site
		p.push_back(_Owner->XPos[i]);
		p.push_back(_Owner->XPos[j]);
		p.push_back(0.0);

		// Loop over velocities (ignore rest distribution)
		for (int vel = 0; vel < L_NUM_VELS - 1; vel++) {

			// Compute destination coordinates
			dest_i = (i + c[0][vel] + g->N_lim) % g->N_lim;
			dest_j = (j + c[1][vel] + g->M_lim) % g->M_lim;

			// Position of destination site
			std::vector<double> ppr;
			ppr.push_back(_Owner->XPos[dest_i]);
			ppr.push_back(_Owner->XPos[dest_j]);
			ppr.push_back(0.0);			

			std::vector<double> r = GridUtils::subtract(ppr,p);	// Length of streaming vector			
    
			// Compute parameters of intersection
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

				// On first pass, set to valid value
				if (Q[vel][storeID] == -1) Q[vel][storeID] = std::numeric_limits<double>::max();
				
				if (q < Q[vel][storeID]) {

					// Set outgoing Q value
					Q[vel][storeID] = q;

					// Incoming Q value (in destination store at opposite direction) is 1 minus the incoming
					Q[GridUtils::getOpposite(vel) + L_NUM_VELS][storeID] = 1 - q;
				}
        
			}

			// Maybe include the collinear case if it cuts through a vertex?
		}
	}
}
