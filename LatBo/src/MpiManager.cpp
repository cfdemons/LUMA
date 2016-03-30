#include "../inc/stdafx.h"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include "../inc/definitions.h"
#include "../inc/MpiManager.h"
#include "../inc/GridObj.h"
#include "../inc/globalvars.h"

// Static declarations
MpiManager* MpiManager::me;
std::ofstream* MpiManager::logout;
GridObj* MpiManager::Grids;

// *************************************************************************************************** //
// Constructor (private)
MpiManager::MpiManager(void)
{
	// Resize buffer arrays based on number of MPI directions
	f_buffer_send.resize(MPI_dir, std::vector<double>(0));
	f_buffer_recv.resize(MPI_dir, std::vector<double>(0));
}

// Destructor
MpiManager::~MpiManager(void)
{

	// Close the logfile
	MpiManager::logout->close();
}

// *************************************************************************************************** //
// Instance creator
MpiManager* MpiManager::getInstance() {

	if (!me) me = new MpiManager;	// Private construction
	return me;						// Return pointer to new object

}

// *************************************************************************************************** //
// Const data member initialised outside class definition
// Define for 3D where first 8 mimic the 2D ones. Opposites are simply the next or previous column in the array.
const int MpiManager::MPI_cartlab[3][26] =
	{
		{1, -1,  1, -1,	 0,  0, -1,  1,		0,  0,		1, -1,  1, -1,  0,  0, -1,  1, -1,  1, -1,  1,  0,  0,  1, -1},
		{0,  0,  1, -1,  1, -1,  1, -1,		0,  0,		0,  0,  1, -1,  1, -1,  1, -1,  0,  0, -1,  1, -1,  1, -1,  1},
		{0,  0,  0,  0,  0,  0,  0,  0,		1, -1,		1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1}
	};

// *************************************************************************************************** //

// Initialisation routines
void MpiManager::mpi_init( ) {
	
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

	// Open logfile now my_rank has been assigned
	logout->open( GridUtils::path_str + "/mpi_log_rank" + std::to_string(my_rank) + ".out", std::ios::out );

#ifdef MPI_VERBOSE
	// Write out coordinates
	*MpiManager::logout << "Coordinates on rank " << my_rank << " are (";
	for (size_t d = 0; d < dims; d++) {
		*MpiManager::logout << "\t" << MPI_coords[d];
	}
	*MpiManager::logout << "\t)" << std::endl;
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
	*MpiManager::logout << "My rank is " << my_rank << ". There are " << num_ranks << " ranks." << std::endl;
#endif

	// Get neighbour ID //
	// Cyclical data transfer i.e. from 0 to 1, 1 to 2, 3 to 4 etc... so use MPI_Sendrecv

	// Loop over grid direction
	for (int dir = 0; dir < MPI_dir; dir++) {

		MPI_Barrier(my_comm);

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
		*MpiManager::logout << "Neighbour in direction " << dir << " of rank " << my_rank << " is (";
		for (size_t d = 0; d < dims; d++) {
			*MpiManager::logout << "\t" << neighbour_coords[d][dir];
		}
		*MpiManager::logout << "\t): Rank " << neighbour_rank[dir] << std::endl;

#ifdef USE_CUSTOM_MPI_SIZES
	// If using custom sizes, user must set the Zcores to 1
	if (dims == 2 && Zcores != 1) {
		std::cout << "Error: See Log File" << std::endl;
		*MpiManager::logout << "Error: Zcores must be set to 1 when using custom MPI sizes in 2D. Exiting." << std::endl;
		MpiManager::logout->close();
		MPI_Finalize();
		exit(10000);
	}
#endif
#endif

	}

	// End Initialisation //

	return;
}

// *************************************************************************************************** //
// Routine to do the domain decomposition and generate parameters for GridObj constructors
void MpiManager::mpi_gridbuild( ) {

	// Global physical dimensions
	double Lx = b_x - a_x;
	double Ly = b_y - a_y;
	double Lz = b_z - a_z;
	double dx = 2 * (Lx / (2 * static_cast<double>(N)));
	double dy = 2 * (Ly / (2 * static_cast<double>(M)));
	double dz = 2 * (Lz / (2 * static_cast<double>(K)));

	// Compute required local grid size
	// Loop over dimensions
	for (size_t d = 0; d < dims; d++) {

		if (MPI_dims[d] == 1) {
			// If only 1 rank in this direction local grid is same size a global grid
			local_size.push_back( global_dims[d] );

		} else if ( fmod(static_cast<double>(global_dims[d]) , static_cast<double>(MPI_dims[d])) ) {
			// If number of cores doesn't allow exact division of grid sites, exit.
			std::cout << "Error: See Log File" << std::endl;
			*MpiManager::logout << "Grid cannot be divided evenly among the cores. Exiting." << std::endl;
			MpiManager::logout->close();
			MPI_Finalize();
			exit(10000);

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
	global_edge_ind.resize( 6, std::vector<int>(num_ranks) );
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

	// Find global indices of edges of coarse grid excluding the overlap
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
	*MpiManager::logout << "Grid size on rank " << my_rank << " is (";
	for (size_t d = 0; d < dims; d++) {
		*MpiManager::logout << "\t" << local_size[d];
	}
	*MpiManager::logout << "\t)" << std::endl;

	*MpiManager::logout << "Limits of the grid (indices) and (position) are (" <<
		global_edge_ind[0][my_rank] << "-" << global_edge_ind[1][my_rank] <<
		", " << global_edge_ind[2][my_rank] << "-" << global_edge_ind[3][my_rank] <<
		", " << global_edge_ind[4][my_rank] << "-" << global_edge_ind[5][my_rank] <<
		"), (" << global_edge_pos[0][my_rank] << "-" << global_edge_pos[1][my_rank] <<
		", " << global_edge_pos[2][my_rank] << "-" << global_edge_pos[3][my_rank] <<
		", " << global_edge_pos[4][my_rank] << "-" << global_edge_pos[5][my_rank] <<
		")" << std:: endl;
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

			std::cout << "Error: See Log File" << std::endl;
			*MpiManager::logout << "Error: Block sizes have been specified in the wrong order, faces do not line up. Exiting." << std::endl;

			// Tell user size it should be
			*MpiManager::logout <<
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

			MpiManager::logout->close();
			MPI_Finalize();
			exit(10000);

		 }

#else

	// 2D check
	if ( (	yRankSize[neighbour_rank[0]] != local_size[1]-2 || yRankSize[neighbour_rank[1]] != local_size[1]-2
		 ) || (
			xRankSize[neighbour_rank[4]] != local_size[0]-2 || xRankSize[neighbour_rank[5]] != local_size[0]-2
		 )
		) {

			std::cout << "Error: See Log File" << std::endl;
			*MpiManager::logout << "Error: Block sizes have been specified in the wrong order, faces do not line up. Exiting." << std::endl;

			// Tell user size it should be
			*MpiManager::logout <<
				" Y (left/right): " <<
				yRankSize[neighbour_rank[0]] << " needed " << local_size[1]-2 << ", " <<
				yRankSize[neighbour_rank[1]] << " needed " << local_size[1]-2 << ", " <<
				" X (up/down): " <<
				xRankSize[neighbour_rank[4]] << " needed " << local_size[0]-2 << ", " <<
				xRankSize[neighbour_rank[5]] << " needed " << local_size[0]-2;

			MpiManager::logout->close();
			MPI_Finalize();
			exit(10000);

		 }

#endif
#endif

}

// *************************************************************************************************** //
// Writes out f-buffers
void MpiManager::mpi_writeout_buf( std::string filename, int dir ) {

	std::ofstream rankout;
	rankout.open(filename.c_str(), std::ios::out);

	rankout << "f_buffer_send is of size " << f_buffer_send[dir].size() << " with values: " << std::endl;
	for (size_t v = 0; v < f_buffer_send[dir].size(); v++) {
		rankout << f_buffer_send[dir][v] << std::endl;
	}

	rankout << "f_buffer_recv is of size " << f_buffer_recv[dir].size() << " with values: " << std::endl;
	for (size_t v = 0; v < f_buffer_recv[dir].size(); v++) {
		rankout << f_buffer_recv[dir][v] << std::endl;
	}

	rankout.close();

	return;
}
// *************************************************************************************************** //
void MpiManager::mpi_communicate(int lev, int reg) {

	// Wall clock variables
	clock_t t_start, t_end, secs;

	// Tag
	int TAG;

	// Get grid object
	GridObj* Grid = NULL;
	GridUtils::getGrid(Grids, lev, reg,  Grid);


	///////////////////////
	// MPI Communication //
	///////////////////////

	/* Must be able to send data but not expect to receive with using MPI on multi-grid configurations
	* so cannot use MPI_Sendrecv_replace() and a single buffer as we did before. Use double buffer 
	* call instead. Also, not every rank will need to send-receive on every grid level as it might 
	* not have a sub-grid so need to handle sends and receives separately to ensure they only communicate 
	* with grids that are expecting it.
	*
	* IMPORTANT: MPI_Barrier() calls synchronise the entire topology. If these calls are made on MPI 
	* communications on a particular grid that only exists on some ranks, the sub-time steps on each rank 
	* will be out of sync. Need to allow the blocking nature of the send and receive calls to force correct 
	* synchronisation between processes and only call barriers outside the grid scope.
	*
	* For each sending direction, pack and load a message into the message queue 
	* for the destination rank with tag associated with direction.
	* Then for each receive direction, pull message with correct tag from the queue 
	* and unpack.
	*
	* In order to do this, need non-blocking send and receive calls and each needs
	* their own buffer to store the information which cannot be touched until the 
	* send is completed, hence this implementation carries a bigger memeory requirement
	* as buffer reuse through the direction loop is not possible.
	* Although MPI_Bsend() will do something similar it relies on creating and filling 
	* MPI background buffers which might have limited resources and which is slower so 
	* we use the MPI Manager class to hold the buffer in house. */

	// Start the clock
	t_start = clock();

	// Loop over directions in Cartesian topology
	for (int dir = 0; dir < MPI_dir; dir++) {

		/* Create a unique tag based on level (< 32), region (< 10) and direction (< 100).
		 * MPICH limits state that tag value cannot be greater than 32767 */
		TAG = ((Grid->level + 1) * 1000) + ((Grid->region_number + 1) * 100) + dir;

#ifdef MPI_VERBOSE
		*MpiManager::logout << "Processing Message with Tag --> " << TAG << std::endl;
#endif

		////////////////////////////
		// Resize and Pack Buffer //
		////////////////////////////

		// Adjust buffer size
		for (MpiManager::buffer_struct bufs : buffer_send_info) {
			if (bufs.level == Grid->level && bufs.region == Grid->region_number) {
				f_buffer_send[dir].resize(bufs.size[dir] * nVels);
			}
		}

		// Only pack and send if required
		if (f_buffer_send[dir].size()) {

			// Pass direction and Grid by reference and pack if required
			mpi_buffer_pack( dir, Grid );
		

			///////////////
			// Post Send //
			///////////////

#ifdef MPI_VERBOSE
			*MpiManager::logout << "L" << Grid->level << "R" << Grid->region_number << " -- Direction " << dir 
								<< " -->  Posting Send for " << f_buffer_send[dir].size() / nVels
								<< " sites to Rank " << neighbour_rank[dir] << "." << std::endl;
#endif
			// Post send message to message queue
			MPI_Isend( &f_buffer_send[dir].front(), f_buffer_send[dir].size(), MPI_DOUBLE, neighbour_rank[dir], TAG, my_comm, &request );

#ifdef MPI_VERBOSE
			*MpiManager::logout << "Direction " << dir << " --> Send Posted." << std::endl;
#endif

		}

		// Find opposite direction (neighbour it receives from)
		int opp_dir = mpi_getOpposite(dir);
		
		// Resize the receive buffer
		for (MpiManager::buffer_struct bufr : buffer_recv_info) {
			if (bufr.level == Grid->level && bufr.region == Grid->region_number) {
				f_buffer_recv[dir].resize(bufr.size[dir] * nVels);
			}
		}


		///////////////////
		// Fetch Message //
		///////////////////

		if (f_buffer_recv[dir].size()) {

#ifdef MPI_VERBOSE
			*MpiManager::logout << "L" << Grid->level << "R" << Grid->region_number << " -- Direction " << dir 
								<< " -->  Fetching message for " << f_buffer_recv[dir].size() / nVels	
								<< " sites from Rank " << neighbour_rank[opp_dir] << "." << std::endl;
#endif

			// Use a blocking receive call if required
			MPI_Recv( &f_buffer_recv[dir].front(), f_buffer_recv[dir].size(), MPI_DOUBLE, neighbour_rank[opp_dir], TAG,	my_comm, &stat );
		

#ifdef MPI_VERBOSE
			*MpiManager::logout << "Direction " << dir << " --> Received." << std::endl;

			// Write out buffers
			*MpiManager::logout << "L" << Grid->level << "R" << Grid->region_number << " -- Direction " << dir << " SUMMARY"
								<< " -- Sent " << f_buffer_send[dir].size() / nVels << " to " << neighbour_rank[dir] 
								<< ": Received " << f_buffer_recv[dir].size() / nVels << " from " << neighbour_rank[opp_dir] << std::endl;
			std::string filename = GridUtils::path_str + "/mpiBuffer_Rank" + std::to_string(MpiManager::my_rank) + "_Dir" + std::to_string(dir) + ".out";
			mpi_writeout_buf(filename, dir);
#endif

			///////////////////////////
			// Unpack Buffer to Grid //
			///////////////////////////

			// Pass direction and Grid by reference
			mpi_buffer_unpack( dir, Grid );

		}

	}


	// Print Time of MPI comms
	t_end = clock();
	secs = t_end - t_start;

	// Update average MPI overhead time for this particular grid
	Grid->timeav_mpi_overhead *= (Grid->t-1);
	Grid->timeav_mpi_overhead += ((double)secs)/CLOCKS_PER_SEC;
	Grid->timeav_mpi_overhead /= Grid->t;

#ifdef TEXTOUT
	if (Grid->t % out_every == 0) {
		*GridUtils::logfile << "Writing out to <Grids.out>" << std::endl;
		Grid->io_textout("POST MPI COMMS");
	}
#endif

	if (Grid->t % out_every == 0) {
		// Performance Data
		*GridUtils::logfile << "MPI overhead taking an average of " << Grid->timeav_mpi_overhead*1000 << "ms" << std::endl;
	}


}

// *************************************************************************************************** //
// Sets pointer to grid hierarchy and finds sizes of buffers for this rank across all grids -- 
// to be called post-initialisation.
void MpiManager::mpi_buffer_size() {

	*GridUtils::logfile << "Pre-computing buffer sizes for MPI..." << std::endl;

	/* For each grid in the hierarchy find communicating edges and store the buffer size.
	 * The data are arranged:
	 *
	 *		BufSize[direction][level][region]
	 *
	 * where direction is identified by the MPI labelling at the top of this file.
	 * A zero buffer size indicates that this edge does not communicate with a neighbour rank. */

	// Loop through levels and regions
	GridObj* g;	// Pointer to a GridObj
	for (int l = 0; l <= NumLev; l++) {
		for (int r = 0; r < NumReg; r++) {

			// Null the pointer
			g = NULL;

			// Get pointer to grid
			GridUtils::getGrid(Grids, l, r, g);

			// If the pointer is not updated then it doesn't exist on this rank
			if (g == NULL) {
				continue;	// Try find next grid in loop
			}

			// Expand buffer info arrays by one and add grid ID info
			buffer_send_info.emplace_back();		buffer_recv_info.emplace_back();
			buffer_send_info.back().level = l;		buffer_recv_info.back().level = l;
			buffer_send_info.back().region = r;		buffer_recv_info.back().region = r;

			/* To allow for the fact that not every edge needs to be communicated,
			 * we must account for 4 states in each Cartesian direction:
			 *
			 * 1. The x_min edge of the lattice needs communicating;
			 * 2. The x_max edge of the lattice needs communicating;
			 * 3. Both the x_min and the x_max edges of the lattice need communicating;
			 * 4. Neither the x_min nor the x_max edges of the lattice need communicating;
			 * 
			 * Therefore in 2D there are 4^2 = 16 permutations and in 3D there are 4^3 = 64 permutations.
			 *
			 * We can do this by looping over the appropriate region of the grid and checking whether 
			 * a site lies in the associated sender layer.
			 * If so increment the buffer size counter. */

			// Call send and recv size finding routines
			mpi_buffer_size_send(g);
			mpi_buffer_size_recv(g);


			// Write out buffer sizes
#ifdef MPI_VERBOSE
		*MpiManager::logout << "[" << buffer_send_info.back().level << "," << buffer_send_info.back().region << "]" << '\t';
		for (int i = 0; i < MPI_dir; i++) {
			*MpiManager::logout << buffer_send_info.back().size[i] << '\t';
		}
		*MpiManager::logout << std::endl;

		*MpiManager::logout << "[" << buffer_recv_info.back().level << "," << buffer_recv_info.back().region << "]" << '\t';
		for (int i = 0; i < MPI_dir; i++) {
			*MpiManager::logout << buffer_recv_info.back().size[i] << '\t';
		}
		*MpiManager::logout << std::endl;
#endif


		}
	}

	

}
// *************************************************************************************************** //
// Get opposite direction in MPI cartesian topology
int MpiManager::mpi_getOpposite(int direction) {

	/*	If direction is even, then opposite is direction+1.
		If direction is odd, then opposite is direction-1.
		e.g. direction 0 (+x direction) has opposite 1 (-x direction) --> +1
		however, direction 1 has opposite 0 --> -1
		Hence we can add (-1 ^ direction) so it alternates between +/-1
	*/

		return direction + (int)pow(-1,direction);

}