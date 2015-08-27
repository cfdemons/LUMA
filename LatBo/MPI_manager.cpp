#include "stdafx.h"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include "definitions.h"
#include "MPI_manager.h"
#include "GridObj.h"


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
		{0,  0,  1, -1,  1, -1,  1, -1,		0,  0,		0,  0,  1, -1,  1, -1,  1, -1,  0,  0, -1,  1, -1,  1, -1,  1},
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
	logout.open( "./Output/mpiLog_Rank_" + std::to_string(my_rank) + ".out", std::ios::out );

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
	logout.open( "./Output/mpiLog_Rank_" + std::to_string(my_rank) + ".out", std::ios::out | std::ios::app );
	logout << "My rank is " << my_rank << ". There are " << num_ranks << " ranks." << std::endl;
	logout.close();
#endif

	// Get neighbour ID //
	// Cyclical data transfer i.e. from 0 to 1, 1 to 2, 3 to 4 etc... so use MPI_Sendrecv

	// Declare array for storing coordinates of neighbour
	int neighbour_coords[dims];

	// Loop over grid direction
	for (int dir = 0; dir < MPI_dir; dir++) {

		MPI_Barrier(my_comm);
#ifdef MPI_VERBOSE
		logout.open( "./Output/mpiLog_Rank_" + std::to_string(my_rank) + ".out", std::ios::out | std::ios::app );

		if (my_rank == 0) {
			logout<< "\nDirection = " << dir << std::endl;
		}
#endif

		// Get coordinates of neighbour (taking into account periodic structure)
		for (size_t d = 0; d < dims; d++) {
			neighbour_coords[d] = (MPI_coords[d] + MPI_cartlab[d][dir] + MPI_dims[d]) % MPI_dims[d];
		}		
		
		// Get rank of neighbour and build vector of ranks
		int tmp;
		MPI_Cart_rank(my_comm, neighbour_coords, &tmp);
		neighbour_rank[dir] = tmp;

		MPI_Barrier(my_comm);

#ifdef MPI_VERBOSE
		// Print out current neighbour coordinates and rank
		logout << "Neighbour of rank " << my_rank << " is (";
		for (size_t d = 0; d < dims; d++) {
			logout << "\t" << neighbour_coords[d];
		}
		logout << "\t): Rank " << neighbour_rank[dir] << std::endl;

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
			if (my_rank == 0) {
				std::cout << "Grid cannot be divided evenly among the cores. Exiting." << std::endl;
			}
			MPI_Finalize();
			exit(EXIT_FAILURE);
		} else {
			// Else, find local grid size
			local_size.push_back( (global_dims[d]/MPI_dims[d]) + 2 ); // Simple uniform decomposition + overlap
		}
	}
	
	// Find global indices of edges of grid excluding the overlap
	global_edge_ind.resize( 6, std::vector<unsigned int>(num_ranks) );
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
	
	// Find global positions of edges of grid excluding the overlap
	global_edge_pos.resize( 6, std::vector<double>(num_ranks) );
	global_edge_pos[1][my_rank] = (Lx / Xcores) * (MPI_coords[0] + 1);
	global_edge_pos[0][my_rank] = global_edge_pos[1][my_rank] - (Lx / Xcores);
	global_edge_pos[3][my_rank] = (Ly / Ycores) * (MPI_coords[1] + 1);
	global_edge_pos[2][my_rank] = global_edge_pos[3][my_rank] - (Ly / Ycores);
#if (dims == 3)
	global_edge_pos[5][my_rank] = (Lz / Zcores) * (MPI_coords[2] + 1);
	global_edge_pos[4][my_rank] = global_edge_pos[5][my_rank] - (Lz / Zcores);
#else
	global_edge_pos[5][my_rank] = b_z;
	global_edge_pos[4][my_rank] = a_z;
#endif


	MPI_Barrier(my_comm);

#ifdef MPI_VERBOSE
	// Write out the Grid size vector
	std::ofstream logout;
	logout.open( "./Output/mpiLog_Rank_" + std::to_string(my_rank) + ".out", std::ios::out | std::ios::app );
	logout << "Grid size on rank " << my_rank << " is (";
	for (size_t d = 0; d < dims; d++) {
		logout << "\t" << local_size[d];
	}
	logout << "\t)" << std::endl;
	logout.close();
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