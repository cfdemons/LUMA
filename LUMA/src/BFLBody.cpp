/*
* --------------------------------------------------------------
*
* ------ Lattice Boltzmann @ The University of Manchester ------
*
* -------------------------- L-U-M-A ---------------------------
*
* Copyright 2019 The University of Manchester
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
* http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.*
*/

#include "../inc/stdafx.h"
#include "../inc/BFLBody.h"
#include "../inc/PCpts.h"
#include "../inc/GridObj.h"


// Implementation of BFL body class //

/******************************************************************************/
/// Default constructor
BFLBody::BFLBody(void)
{
}

/******************************************************************************/
/// Default destructor
BFLBody::~BFLBody(void)
{
}

/******************************************************************************/
/// \brief	Initialiser for a general BFL body.
///
///			This initialiser performs debugging IO and wraps the calling of the
///			labelling and Q computation. It is called immediately after the 
///			appropriate constructor by all derived constructors.
void BFLBody::initialise()
{

	// Write out marker data now body has been built
#ifdef L_BFL_DEBUG
	std::ofstream file;
	file.open(GridUtils::path_str + "/BflMarkerData_Rank" + std::to_string(GridUtils::safeGetRank()) + ".out", std::ios::out);
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
		_Owner->LatTyp(m.supp_i[0], m.supp_j[0], m.supp_k[0], M_lim, K_lim) = eBFL;
	}

	// Close Body //
	*GridUtils::logfile << "ObjectManagerBFL: Checking surface integrity..." << std::endl;
	enforceSurfaceClosure();

	// Write out closed body marker data
#ifdef L_BFL_DEBUG
	file;
	file.open(GridUtils::path_str + "/BflMarkerDataClosed_Rank" + std::to_string(GridUtils::safeGetRank()) + ".out", std::ios::out);
	file.precision(L_OUTPUT_PRECISION);
	for (size_t n = 0; n < markers.size(); n++) {
		file << std::to_string(n) << ", " <<
			markers[n].position[0] << ", " << markers[n].position[1] << ", " << markers[n].position[2] << ", " <<
			markers[n].supp_i[0] << ", " << markers[n].supp_j[0] << ", " << markers[n].supp_k[0] << std::endl;
	}
	file.close();
#endif

	// Compute Q //
	*GridUtils::logfile << "ObjectManagerBFL: Computing Q..." << std::endl;

	// Initialise Q stores to the "invalid" value
	Q.resize(L_NUM_VELS * markers.size(), -1.0);

	// Loop over local grid and inspect the streaming operations
	for (int i = 0; i < N_lim; i++) {
		for (int j = 0; j < M_lim; j++) {
			for (int k = 0; k < K_lim; k++) {

				// If site is a BFL voxel
				if (_Owner->LatTyp(i, j, k, M_lim, K_lim) == eBFL) {

					// Compute Q for all stream vectors storing on source voxel BFL marker
#if (L_DIMS == 3)
					computeQ(i, j, k, _Owner);
#else
					computeQ(i, j, _Owner);
#endif

				}
			}
		}
	}

	// Computation of Q complete
	*GridUtils::logfile << "ObjectManagerBFL: Q computation complete." << std::endl;

	// Set valid markers
	validMarkers = GridUtils::onespace(0, static_cast<int>(markers.size()) - 1);

	// Write out Q values for each marker
#ifdef L_BFL_DEBUG
	file.open(GridUtils::path_str + "/BflMarkerQs_Rank" + std::to_string(GridUtils::safeGetRank()) + ".out", std::ios::out);
	file.precision(L_OUTPUT_PRECISION);
	for (size_t n = 0; n < markers.size(); ++n)
	{
		for (size_t v = 0; v < L_NUM_VELS; ++v)
		{		
			file << Q[v + L_NUM_VELS * n] << '\t';
		}
		file << std::endl;
	}
	file.close();
#endif

}

/******************************************************************************/
/// \brief Custom constructor to populate body from array of points.
/// \param g		hierarchy pointer to grid hierarchy
/// \param bodyID	ID of body in array of bodies.
/// \param _PCpts	pointer to point cloud data
BFLBody::BFLBody(GridObj* g, int bodyID, PCpts* _PCpts) 
	: Body(g, bodyID, _PCpts)
{
	initialise();
}


/******************************************************************************/
/// \brief 	Custom constructor for building prefab filament
/// \param g				hierarchy pointer to grid hierarchy
/// \param bodyID			ID of body in array of bodies.
/// \param start_position	start position of base of filament
/// \param length			length of filament
/// \param angles			angle of filament
BFLBody::BFLBody(GridObj* g, int bodyID, std::vector<double> &start_position,
	double length, std::vector<double> &angles)
	: Body(g, bodyID, start_position, length, angles)
{
	initialise();
}


/******************************************************************************/
/// \brief 	Custom constructor for building prefab circle/sphere
/// \param g				hierarchy pointer to grid hierarchy
/// \param bodyID			ID of body in array of bodies.
/// \param centre_point		centre point of circle
/// \param radius			radius of circle
BFLBody::BFLBody(GridObj* g, int bodyID, std::vector<double> &centre_point,
	double radius)
	: Body(g, bodyID, centre_point, radius)
{
	initialise();
}


/******************************************************************************/
/// \brief 	Custom constructor for building square/cuboid
/// \param g				hierarchy pointer to grid hierarchy
/// \param bodyID			ID of body in array of bodies.
/// \param centre_point		centre point of square
/// \param dimensions		dimensions of square
/// \param angles			angle of square
BFLBody::BFLBody(GridObj* g, int bodyID, std::vector<double> &centre_point,
	std::vector<double> &dimensions, std::vector<double> &angles)
	: Body(g, bodyID, centre_point, dimensions, angles)
{
	initialise();
}


/******************************************************************************/
/// \brief 	Custom constructor for building plate
/// \param g				hierarchy pointer to grid hierarchy
/// \param bodyID			ID of body in array of bodies.
/// \param centre_point		centre point of plate
/// \param length			length of plate
/// \param width			width of plate
/// \param angles			angle of plate
BFLBody::BFLBody(GridObj* g, int bodyID, std::vector<double> &centre_point,
	double length, double width, std::vector<double> &angles)
	: Body(g, bodyID, centre_point, length, width, angles)
{
	initialise();
}

/******************************************************************************/
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
void BFLBody::computeQ(int i, int j, int k, GridObj* g)
{

	// Declarations
	int dest_i, dest_j, dest_k, storeID;
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


	// TODO: Update under the restrictions we have done for 2D //

	// Get marker data associated with this local site
	m_data = getMarkerData(g->XPos[i], g->YPos[j], g->ZPos[k]);

	storeID = m_data->ID;
	delete m_data;

	// Get list of IDs of neighbour vertices for plane construction
	std::vector<int> V;
	for (int ii = i - 1; ii <= i + 1; ii++) {
		for (int jj = j - 1; jj <= j + 1; jj++) {
			for (int kk = k - 1; kk <= k + 1; kk++) {

				// If indices are valid (include itself in combinations)
				if (	ii >= 0 && ii < g->N_lim
					&&	jj >= 0 && jj < g->M_lim
					&&	kk >= 0 && kk < g->K_lim
					)
				{

					// Fetch data if available
					m_data = getMarkerData(g->XPos[ii], g->YPos[jj], g->ZPos[kk]);

					// If data valid, then store ID
					if (m_data->isValid()) V.push_back(m_data->ID);
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
				if (Q[vel + L_NUM_VELS * storeID] == -1) Q[vel + L_NUM_VELS * storeID] = std::numeric_limits<double>::max();

				if (q < Q[vel + L_NUM_VELS * storeID]) {

					// Set outgoing Q value
					Q[vel + L_NUM_VELS * storeID] = q;

				}
			}
		}
	}

}

/******************************************************************************/
/// \brief	Routine to compute wall distance Q.
///
///			Computes Q values in 2D at a given local voxel for each application of 
///			the BFL BC. Performs a line-line intersection algorithm for each line 
///			segment either side of the voxel marker.
///
/// \param i local i-index of BFL voxel
/// \param j local j-index of BFL voxel
/// \param g pointer to owner grid
void BFLBody::computeQ(int i, int j, GridObj* g)
{

	/* For each possible line extending from the marker voxel to an 
	 * immediate neighbour we can check for an intersection between the 
	 * stream vector and the constructed line. The distance of intersection
	 * is the value of q. Since the inclusion of the closure model, we can limit
	 * our selection to just vertical and horizontal adjacent markers. */

	// Declarations
	int dest_i, dest_j;
	double s, t, s1_x, s1_y, s2_x, s2_y;
	MarkerData* m_data;
	
	// Get marker data associated with this local site
	m_data = getMarkerData(g->XPos[i], g->YPos[j], g->ZPos[0]);

	int storeID = m_data->ID;
	delete m_data;


	// Get IDs of vertical and horizontal neighbour vertices for line construction
	std::vector< std::pair<int, int> > combo;
	for (int ii = i - 1; ii <= i + 1; ii++) {
		for (int jj = j - 1; jj <= j + 1; jj++) {

			// Must be on grid and only consider v/h neighbours (exclude self)
			if	(	ii >= 0 && ii < g->N_lim &&	
					jj >= 0 && jj < g->M_lim &&
					(ii == i || jj == j) &&
					!(ii == i && jj == j)
				)
			{			

				// Fetch data if available
				m_data = getMarkerData(g->XPos[ii], g->YPos[jj], g->ZPos[0]);

				// If data valid, then store ID
				if (m_data->isValid()) combo.push_back(std::pair<int,int>(storeID, m_data->ID));

				// Clean-up after getMarkerData call
				delete m_data;

			}

		}
	}

	// Can only continue at least 1 pair
	if (combo.size() == 0) return;

	// Get position of marker in this cell
	std::vector<double> q;
	q.push_back(markers[storeID].position[0]);
	q.push_back(markers[storeID].position[1]);
	q.push_back(markers[storeID].position[2]);

	// Loop through valid marker combinations
	for (std::pair<int,int>& line : combo) {

		/* Perform line intersection test according to 2nd answer on:
		 * http://stackoverflow.com/questions/563198/how-do-you-detect-where-two-line-segments-intersect
		 */

		std::vector<double> qps;	// Position of next marker
		qps.push_back(markers[line.second].position[0]);
		qps.push_back(markers[line.second].position[1]);
		qps.push_back(markers[line.second].position[2]);

		std::vector<double> p;	// Position of source site
		p.push_back(_Owner->XPos[i]);
		p.push_back(_Owner->YPos[j]);
		p.push_back(markers[line.first].position[2]);	// In 2D must be in the same plane

		// Loop over velocities (ignore rest distribution)
		for (int vel = 0; vel < L_NUM_VELS - 1; vel++) {

			// Compute destination coordinates
			dest_i = (i + c[0][vel] + g->N_lim) % g->N_lim;
			dest_j = (j + c[1][vel] + g->M_lim) % g->M_lim;

			// Position of destination site
			std::vector<double> ppr;
			ppr.push_back(_Owner->XPos[dest_i]);
			ppr.push_back(_Owner->YPos[dest_j]);
			ppr.push_back(markers[line.first].position[2]);

			// Compute lengths of lines
			s1_x = ppr[0] - p[0];
			s1_y = ppr[1] - p[1];
			s2_x = qps[0] - q[0];
			s2_y = qps[1] - q[1];

			// Cross products
			s = (-s1_y * (p[0] - q[0]) + s1_x * (p[1] - q[1])) / (-s2_x * s1_y + s1_x * s2_y);
			t = (s2_x * (p[1] - q[1]) - s2_y * (p[0] - q[0])) / (-s2_x * s1_y + s1_x * s2_y);

			// Test for intersection
			if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
			{
				// Lines intersect at p + ts_1
				std::vector<double> intersect;
				intersect.push_back(p[0] + t * s1_x);
				intersect.push_back(p[1] + t * s1_y);
				intersect.push_back(p[2]);	// Again same plane in 2D
        
				// Compute Q (normalised)
				double wall_distance = GridUtils::vecnorm(GridUtils::subtract(intersect, p)) / GridUtils::vecnorm(GridUtils::subtract(ppr, p));

				// On first pass, set to valid value
				if (Q[vel + L_NUM_VELS * storeID] == -1)
					Q[vel + L_NUM_VELS * storeID] = std::numeric_limits<double>::max();
				
				// Store if lower than previous
				if (wall_distance < Q[vel + L_NUM_VELS * storeID]) {

					// Set outgoing Q value
					Q[vel + L_NUM_VELS * storeID] = wall_distance;

				}
        
			}

			// Maybe include the collinear case if it cuts through a vertex?
		}
	}
}

/******************************************************************************/
///	\brief	Ensures that BFL body is represented by markers such that there are
///			no holes in the surface caused by the voxel grid filter.
///
///			There is a perculiar problem with BFL bodies where if the boundary 
///			represented by the point cloud passes diagonally between two cells
///			such that there is no point in any vertical or horizontal cell
///			connecting them then the wall is not visible in the direction
///			perpendicular to wall as there is no BFL site to trigger the boundary
///			condition and artefacts are generated. Increasing the input cloud 
///			resolution would eventually fix the problem but could require a lot
///			of computational power to initialise.
///			We can work round this by checking the immediate neighbourhood of 
///			every BFL cell to ensure it has a vertical or horizontal connection to
///			another BFL site. If not, then we can project the line between two 
///			diagonal neigbours to add a new marker.
///
void BFLBody::enforceSurfaceClosure()
{
	// Commonly accessed
	int K_lim = _Owner->K_lim;
	int M_lim = _Owner->M_lim;
	int i_neigh, j_neigh, k_neigh, is, js, ks;
	int av_count;
	bool bAdjacentConnectionFound = false, bNewMarkerRequired = false;
	double start_vec[3];
	double len_vec[3];
	double px, py, pz, av_pos_x, av_pos_y, av_pos_z;
	std::vector<int> ijk_point;
	int num_points_in_projection = 20;
	double projection_spacing = 1.0 / num_points_in_projection;

	/* Loop over every marker in the body (don't use a ranged-for as iterators are 
	 * invalidated if reallocation is performed when new marker is added. */
	for (size_t m = 0; m < markers.size(); ++m)
	{
		// Check diagonals
		for (int i =  -1; i <= 1; i += 2)
		{
			for (int j = -1; j <= 1; j += 2)
			{
#if (L_DIMS == 3)
				for (int k = -1; k <= 1; k += 2)
#else
				int k = 0;
#endif
				{

					// Get indices of the diagonal site
					i_neigh = markers[m].supp_i[0] + i;
					j_neigh = markers[m].supp_j[0] + j;
					k_neigh = markers[m].supp_k[0] + k;

					// If diagonal BFL site found
					if (!GridUtils::isOffGrid(i_neigh, j_neigh, k_neigh, _Owner) &&
						_Owner->LatTyp(i_neigh, j_neigh, k_neigh, M_lim, K_lim) == eBFL)
					{
						// Reset flag
						bAdjacentConnectionFound = false;

						// Check the surrounding values for a connecting site
						for (int i_count = 0; i_count < 2; i_count++)
						{
							for (int j_count = 0; j_count < 2; j_count++)
							{
								for (int k_count = 0; k_count < 2; k_count++)
								{
									// Compute neighbour indices
									is = i - i_count * i;
									js = j - j_count * j;
									ks = k - k_count * k;

									// Ignore the diagonal itself and centre site
									if	(	(is == i && js == j && ks == k) ||
											(is == 0 && js == 0 && ks == 0)
										) continue;

									// Check site label of adjacent site and early break if found
									if (!GridUtils::isOffGrid(markers[m].supp_i[0] + is, markers[m].supp_j[0] + js, markers[m].supp_k[0] + ks, _Owner) &&
										_Owner->LatTyp(markers[m].supp_i[0] + is, markers[m].supp_j[0] + js, markers[m].supp_k[0] + ks, M_lim, K_lim) == eBFL)
									{
										bAdjacentConnectionFound = true;
										break;
									}
								}
								if (bAdjacentConnectionFound) break;
							}
							if (bAdjacentConnectionFound) break;
						}

						// If no site found then need to add new marker
						if (!bAdjacentConnectionFound)
						{
							// Reset flag
							bNewMarkerRequired = false;
							av_count = 0;

							// Create vector of points between the current marker and the illegal diagonal
							start_vec[eXDirection] = markers[m].position[eXDirection];
							start_vec[eYDirection] = markers[m].position[eYDirection];
							start_vec[eZDirection] = markers[m].position[eZDirection];

							MarkerData *m_data = getMarkerData(_Owner->XPos[i_neigh], _Owner->YPos[j_neigh], _Owner->ZPos[k_neigh]);
							len_vec[eXDirection] = markers[m_data->ID].position[eXDirection] - start_vec[eXDirection];
							len_vec[eYDirection] = markers[m_data->ID].position[eYDirection] - start_vec[eYDirection];
							len_vec[eZDirection] = markers[m_data->ID].position[eZDirection] - start_vec[eZDirection];
							delete m_data;

							// Start an iterative projection procedure
							while (!bNewMarkerRequired)
							{
								num_points_in_projection *= 2;
								projection_spacing = 1.0 / num_points_in_projection;

								// For each point see if it ends up in a non-BFL voxel and flag
								for (int p = 0; p < num_points_in_projection + 1; p++)
								{
									px = start_vec[eXDirection] + p * projection_spacing * len_vec[eXDirection];
									py = start_vec[eYDirection] + p * projection_spacing * len_vec[eYDirection];
									pz = start_vec[eZDirection] + p * projection_spacing * len_vec[eZDirection];

									GridUtils::getEnclosingVoxel(px, py, pz, _Owner, &ijk_point);

									// If site not in a BFL voxel then new marker required so average position
									if (_Owner->LatTyp(ijk_point[eXDirection], ijk_point[eYDirection], ijk_point[eZDirection], M_lim, K_lim) != eBFL)
									{
										if (!bNewMarkerRequired)
										{
											bNewMarkerRequired = true;
											av_pos_x = px;
											av_pos_y = py;
											av_pos_z = pz;
										}
										else
										{
											av_pos_x *= av_count;
											av_pos_x += px;
											av_pos_x /= av_count + 1;

											av_pos_y *= av_count;
											av_pos_y += py;
											av_pos_y /= av_count + 1;

											av_pos_z *= av_count;
											av_pos_z += pz;
											av_pos_z /= av_count + 1;
										}
										av_count++;
									}
								}

								// Add marker to body and label voxel
								if (bNewMarkerRequired)
								{
									// Add new marker to the end of the array
									addMarker(av_pos_x, av_pos_y, av_pos_z, static_cast<int>(markers.size()));
									_Owner->LatTyp(markers.back().supp_i[0], markers.back().supp_j[0], markers.back().supp_k[0], M_lim, K_lim) = eBFL;
								}


							}	// Iterative projection

						}	// No adjacent connection found

					}	// Diagonally connected to another BFL site

				}	// k diagonal loop
			}	// j diagonal loop
		}	// i diagonal loop

	}	// Loop over markers in the BFL body

} // Method end
/******************************************************************************/