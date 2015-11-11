#include "../inc/stdafx.h"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include "../inc/definitions.h"
#include "../inc/MPI_manager.h"
#include "../inc/GridObj.h"
#include "../inc/globalvars.h"


// Constructor
MPI_manager::MPI_manager(void)
{
}

// Destructor
MPI_manager::~MPI_manager(void)
{
}

// *************************************************************************************************** //
// Const data member initialised outside class definition
// Define for 3D where first 8 mimic the 2D ones. Opposites are simply the next or previous column in the array.
const int MPI_manager::MPI_cartlab[3][26] =
	{
		{1, -1,  1, -1,	 0,  0, -1,  1,		0,  0,		1, -1,  1, -1,  0,  0, -1,  1, -1,  1, -1,  1,  0,  0,  1, -1},
		{0,  0, -1,  1, -1,  1, -1,  1,		0,  0,		0,  0, -1,  1, -1,  1, -1,  1,  0,  0,  1, -1,  1, -1,  1, -1},
		{0,  0,  0,  0,  0,  0,  0,  0,		1, -1,		1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1}
	};

// *************************************************************************************************** //

// Initialisation routines
void MPI_manager::mpi_init( ) {

	// Create communicator and topology
	int MPI_periodic[dims], MPI_reorder;
	MPI_reorder = true;
	MPI_dims[0] = Xcores; MPI_dims[1] = Ycores;
	MPI_periodic[0] = true;	MPI_periodic[1] = true;
#if (dims == 3)
	MPI_dims[2] = Zcores; MPI_periodic[2] = true;
#endif

	MPI_Cart_create(MPI_COMM_WORLD, dims, MPI_dims, MPI_periodic, MPI_reorder, &my_comm);

	// Get Cartesian topology info
	MPI_Comm_rank( my_comm, &my_rank );
	MPI_Comm_size( my_comm, &num_ranks );

	// Store coordinates in the new topology
	MPI_Cart_coords(my_comm, my_rank, dims, MPI_coords);

#ifdef MPI_VERBOSE
	std::ofstream logout;
	logout.open( timeout_str + "/mpiLog_Rank_" + std::to_string(my_rank) + ".out", std::ios::out );

	// Write out coordinates
	logout << "Coordinates on rank " << my_rank << " are (";
	for (size_t d = 0; d < dims; d++) {
		logout << "\t" << MPI_coords[d];
	}
	logout << "\t)" << std::endl;
	logout.close();
#endif

	// Store global grid size
	global_dims[0] = N;
	global_dims[1] = M;
#if (dims == 3)
	global_dims[2] = K;
#else
	global_dims[2] = 1;
#endif

	MPI_Barrier(my_comm);

#ifdef MPI_VERBOSE
	// State my rank
	logout.open( timeout_str + "/mpiLog_Rank_" + std::to_string(my_rank) + ".out", std::ios::out | std::ios::app );
	logout << "My rank is " << my_rank << ". There are " << num_ranks << " ranks." << std::endl;
	logout.close();
#endif

	// Get neighbour ID //
	// Cyclical data transfer i.e. from 0 to 1, 1 to 2, 3 to 4 etc... so use MPI_Sendrecv

	// Loop over grid direction
	for (int dir = 0; dir < MPI_dir; dir++) {

		MPI_Barrier(my_comm);
#ifdef MPI_VERBOSE
		logout.open( timeout_str + "/mpiLog_Rank_" + std::to_string(my_rank) + ".out", std::ios::out | std::ios::app );

		if (my_rank == 0) {
			logout<< "\nDirection = " << dir << std::endl;
		}
#endif

		// Get coordinates of neighbour (taking into account periodic structure)
		int coord_tmp[dims];
		for (size_t d = 0; d < dims; d++) {
			neighbour_coords[d][dir] = (MPI_coords[d] + MPI_cartlab[d][dir] + MPI_dims[d]) % MPI_dims[d];
			coord_tmp[d] = neighbour_coords[d][dir];	// Store single vector for getting neighbour rank
		}

		// Get rank of neighbour and build vector of neighbour ranks
		int tmp;
		MPI_Cart_rank(my_comm, coord_tmp, &tmp);
		neighbour_rank[dir] = tmp;

		MPI_Barrier(my_comm);

#ifdef MPI_VERBOSE
		// Print out current neighbour coordinates and rank
		logout << "Neighbour of rank " << my_rank << " is (";
		for (size_t d = 0; d < dims; d++) {
			logout << "\t" << neighbour_coords[d][dir];
		}
		logout << "\t): Rank " << neighbour_rank[dir] << std::endl;

#ifdef USE_CUSTOM_MPI_SIZES
	// If using custom sizes, user must set the Zcores to 1
	if (dims == 2 && Zcores != 1) {
		std::cout << "Error: See Log File" << std::endl;
		logout << "Error: Zcores must be set to 1 when using custom MPI sizes in 2D. Exiting." << std::endl;
		logout.close();
		MPI_Finalize();
		exit(EXIT_FAILURE);
	}
#endif

		logout.close();
#endif

	}

	// End Initialisation //

	return;
}

// *************************************************************************************************** //
void MPI_manager::mpi_gridbuild( ) {

	// Global physical dimensions
	double Lx = b_x - a_x;
	double Ly = b_y - a_y;
	double Lz = b_z - a_z;
	double dx = 2*(Lx/(2*N));
	double dy = 2*(Ly/(2*M));
	double dz = 2*(Lz/(2*K));

	// Compute required local grid size
	// Loop over dimensions
	for (size_t d = 0; d < dims; d++) {

		if (MPI_dims[d] == 1) {
			// If only 1 rank in this direction local grid is same size a global grid
			local_size.push_back( global_dims[d] );

		} else if ( global_dims[d] % MPI_dims[d] != 0 ) {
			// If number of cores doesn't allow exact division of grid sites, exit.
			std::ofstream logout;
			logout.open( timeout_str + "/mpiLog_Rank_" + std::to_string(my_rank) + ".out", std::ios::out | std::ios::app );
			std::cout << "Error: See Log File" << std::endl;
			logout << "Grid cannot be divided evenly among the cores. Exiting." << std::endl;
			logout.close();
			MPI_Finalize();
			exit(EXIT_FAILURE);

		} else {

			// Else, find local grid size
#ifdef USE_CUSTOM_MPI_SIZES

			// Get grids sizes from the definitions file
			switch (d)
			{
			case 0:
				local_size.push_back( xRankSize[my_rank] + 2 );
				break;
			case 1:
				local_size.push_back( yRankSize[my_rank] + 2 );
				break;
			case 2:
				local_size.push_back( zRankSize[my_rank] + 2 );
				break;

			default:
				break;
			}



#else
			local_size.push_back( (global_dims[d]/MPI_dims[d]) + 2 ); // Simple uniform decomposition + overlap



#endif
		}
	}

	// Find indices and positions of edges of each grid in the global coordinate space excluding the overlap
	global_edge_ind.resize( 6, std::vector<unsigned int>(num_ranks) );
	global_edge_pos.resize( 6, std::vector<double>(num_ranks) );


#ifdef USE_CUSTOM_MPI_SIZES

	// If using custom sizing need to cumulatively establish how far from origin

	int adj_rank;

	// X edges
	global_edge_ind[1][my_rank] = 0;
	for (int i = 0; i < MPI_coords[0] + 1; i++) {
		// Find i-th rank
#if (dims == 3)
		int adj_coords[dims] = {i, MPI_coords[1], MPI_coords[2]};
#else
		int adj_coords[dims] = {i, MPI_coords[1]};
#endif
		MPI_Cart_rank(my_comm, adj_coords, &adj_rank);
		// Add the number of sites on this rank to total
		global_edge_ind[1][my_rank] += xRankSize[adj_rank];
	}
	global_edge_ind[0][my_rank] = global_edge_ind[1][my_rank] - xRankSize[my_rank];

	// Y edges
	global_edge_ind[3][my_rank] = 0;
	for (int i = 0; i < MPI_coords[1] + 1; i++) {
		// Find i-th rank
#if (dims == 3)
		int adj_coords[dims] = {MPI_coords[0], i, MPI_coords[2]};
#else
		int adj_coords[dims] = {MPI_coords[0], i};
#endif
		MPI_Cart_rank(my_comm, adj_coords, &adj_rank);
		// Add the number of sites on this rank to total
		global_edge_ind[3][my_rank] += yRankSize[adj_rank];
	}
	global_edge_ind[2][my_rank] = global_edge_ind[3][my_rank] - yRankSize[my_rank];

#if (dims == 3)
	// Z edges
	global_edge_ind[5][my_rank] = 0;
	for (int i = 0; i < MPI_coords[2] + 1; i++) {
		// Find i-th rank
		int adj_coords[dims] = {MPI_coords[0], MPI_coords[1], i};
		MPI_Cart_rank(my_comm, adj_coords, &adj_rank);
		// Add the number of sites on this rank to total
		global_edge_ind[5][my_rank] += zRankSize[adj_rank];
	}
	global_edge_ind[4][my_rank] = global_edge_ind[5][my_rank] - zRankSize[my_rank];
#else
	global_edge_ind[5][my_rank] = 1;
	global_edge_ind[4][my_rank] = 0;
#endif


	// Using uniform decomposition
#else

	// Find global indices of edges of grid excluding the overlap
	global_edge_ind[1][my_rank] = (N / Xcores) * (MPI_coords[0] + 1);
	global_edge_ind[0][my_rank] = global_edge_ind[1][my_rank] - (N / Xcores);
	global_edge_ind[3][my_rank] = (M / Ycores) * (MPI_coords[1] + 1);
	global_edge_ind[2][my_rank] = global_edge_ind[3][my_rank] - (M / Ycores);
#if (dims == 3)
	global_edge_ind[5][my_rank] = (K / Zcores) * (MPI_coords[2] + 1);
	global_edge_ind[4][my_rank] = global_edge_ind[5][my_rank] - (K / Zcores);
#else
	global_edge_ind[5][my_rank] = 1;
	global_edge_ind[4][my_rank] = 0;
#endif

#endif



	// Find global positions of edges of grid excluding the overlap from the global indices
	for (int d = 0; d < 6; d++) {
		global_edge_pos[d][my_rank] = global_edge_ind[d][my_rank] * dx;
	}
#if (dims != 3)
	global_edge_pos[5][my_rank] = b_z;
	global_edge_pos[4][my_rank] = a_z;
#endif




	MPI_Barrier(my_comm);

#ifdef MPI_VERBOSE
	// Write out the Grid size vector
	std::ofstream logout;
	logout.open( timeout_str + "/mpiLog_Rank_" + std::to_string(my_rank) + ".out", std::ios::out | std::ios::app );
	logout << "Grid size on rank " << my_rank << " is (";
	for (size_t d = 0; d < dims; d++) {
		logout << "\t" << local_size[d];
	}
	logout << "\t)" << std::endl;

	logout << "Limits of the grid (indices) and (position) are (" <<
		global_edge_ind[0][my_rank] << "-" << global_edge_ind[1][my_rank] <<
		", " << global_edge_ind[2][my_rank] << "-" << global_edge_ind[3][my_rank] <<
		", " << global_edge_ind[4][my_rank] << "-" << global_edge_ind[5][my_rank] <<
		"), (" << global_edge_pos[0][my_rank] << "-" << global_edge_pos[1][my_rank] <<
		", " << global_edge_pos[2][my_rank] << "-" << global_edge_pos[3][my_rank] <<
		", " << global_edge_pos[4][my_rank] << "-" << global_edge_pos[5][my_rank] <<
		")" << std:: endl;

	logout.close();
#endif


	// Check my grid size dimensions with the neighbours to make sure it all lines up
#ifdef USE_CUSTOM_MPI_SIZES

	// 3D check
#if (dims == 3)
	if ( (	zRankSize[neighbour_rank[0]] != local_size[2]-2 || zRankSize[neighbour_rank[1]] != local_size[2]-2 ||
			yRankSize[neighbour_rank[0]] != local_size[1]-2 || yRankSize[neighbour_rank[1]] != local_size[1]-2
		 ) || (
			zRankSize[neighbour_rank[4]] != local_size[2]-2 || zRankSize[neighbour_rank[5]] != local_size[2]-2 ||
			xRankSize[neighbour_rank[4]] != local_size[0]-2 || xRankSize[neighbour_rank[5]] != local_size[0]-2
		 ) || (
			xRankSize[neighbour_rank[8]] != local_size[0]-2 || xRankSize[neighbour_rank[9]] != local_size[0]-2 ||
			yRankSize[neighbour_rank[8]] != local_size[1]-2 || yRankSize[neighbour_rank[9]] != local_size[1]-2
		 )
		) {

			std::ofstream logout;
			logout.open(timeout_str + "/mpiLog_Rank_" + std::to_string(my_rank) + ".out", std::ios::out | std::ios::app );
			std::cout << "Error: See Log File" << std::endl;
			logout << "Error: Block sizes have been specified in the wrong order, faces do not line up. Exiting." << std::endl;

			// Tell user size it should be
			logout <<
				" Z (left/right): " <<
				zRankSize[neighbour_rank[0]] << " needed " << local_size[2]-2 << ", " <<
				zRankSize[neighbour_rank[1]] << " needed " << local_size[2]-2 << ", " <<
				" Z (up/down): " <<
				zRankSize[neighbour_rank[4]] << " needed " << local_size[2]-2 << ", " <<
				zRankSize[neighbour_rank[5]] << " needed " << local_size[2]-2 << ", " <<
				" Y (left/right): " <<
				yRankSize[neighbour_rank[0]] << " needed " << local_size[1]-2 << ", " <<
				yRankSize[neighbour_rank[1]] << " needed " << local_size[1]-2 << ", " <<
				" Y (front/back): " <<
				yRankSize[neighbour_rank[8]] << " needed " << local_size[1]-2 << ", " <<
				yRankSize[neighbour_rank[9]] << " needed " << local_size[1]-2 << ", " <<
				" X (up/down): " <<
				xRankSize[neighbour_rank[4]] << " needed " << local_size[0]-2 << ", " <<
				xRankSize[neighbour_rank[5]] << " needed " << local_size[0]-2 << ", " <<
				" X (front/back): " <<
				xRankSize[neighbour_rank[8]] << " needed " << local_size[0]-2 << ", " <<
				xRankSize[neighbour_rank[9]] << " needed " << local_size[0]-2;

			logout.close();
			MPI_Finalize();
			exit(EXIT_FAILURE);

		 }

#else

	// 2D check
	if ( (	yRankSize[neighbour_rank[0]] != local_size[1]-2 || yRankSize[neighbour_rank[1]] != local_size[1]-2
		 ) || (
			xRankSize[neighbour_rank[4]] != local_size[0]-2 || xRankSize[neighbour_rank[5]] != local_size[0]-2
		 )
		) {

			std::ofstream logout;
			logout.open( timeout_str + "/mpiLog_Rank_" + std::to_string(my_rank) + ".out", std::ios::out | std::ios::app );
			std::cout << "Error: See Log File" << std::endl;
			logout << "Error: Block sizes have been specified in the wrong order, faces do not line up. Exiting." << std::endl;

			// Tell user size it should be
			logout <<
				" Y (left/right): " <<
				yRankSize[neighbour_rank[0]] << " needed " << local_size[1]-2 << ", " <<
				yRankSize[neighbour_rank[1]] << " needed " << local_size[1]-2 << ", " <<
				" X (up/down): " <<
				xRankSize[neighbour_rank[4]] << " needed " << local_size[0]-2 << ", " <<
				xRankSize[neighbour_rank[5]] << " needed " << local_size[0]-2;

			logout.close();
			MPI_Finalize();
			exit(EXIT_FAILURE);

		 }

#endif
#endif








}

// *************************************************************************************************** //

void MPI_manager::writeout_buf( std::string filename ) {

	std::ofstream rankout;
	rankout.open(filename.c_str(), std::ios::out);
	rankout << "f_buffer is of size " << f_buffer.size() << std::endl;

		rankout << "\n\nf_buffer Values" << std::endl;
		for (size_t v = 0; v < f_buffer.size(); v++) {

			rankout << f_buffer[v] << std::endl;

		}

		rankout.close();

	return;
}

// *************************************************************************************************** //
