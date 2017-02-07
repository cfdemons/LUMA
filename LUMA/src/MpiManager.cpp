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

// Static declarations
MpiManager* MpiManager::me;

// Although I don't want this to be static VC++ 2013 doesn't support constructor initialiser lists
const int MpiManager::neighbour_vectors[3][26] = 
{
	{ 1, -1, 1, -1, 0, 0, -1, 1, 0, 0, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1, -1, 1, 0, 0, 1, -1 },
	{ 0, 0, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 0, 0, -1, 1, -1, 1, -1, 1 },
	{ 0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1 }
};

// ****************************************************************************
/// Default constructor
MpiManager::MpiManager()
{

#ifdef L_MPI_VERBOSE

	// Create stream for verbose logging
	logout = new std::ofstream();

#else

	logout = nullptr;

#endif

	// Resize buffer arrays based on number of MPI directions
	f_buffer_send.resize(L_MPI_DIRS, std::vector<double>(0));
	f_buffer_recv.resize(L_MPI_DIRS, std::vector<double>(0));

	// Assign the MPI topology size for decomposition
	size_t xSize, ySize, zSize;
	xSize = L_MPI_XCORES;
	ySize = L_MPI_YCORES;
#if (L_DIMS == 3)
	zSize = L_MPI_ZCORES;
#else
	zSize = 1;
#endif
	cRankSizeX.resize(xSize*ySize*zSize);
	cRankSizeY.resize(xSize*ySize*zSize);
	cRankSizeZ.resize(xSize*ySize*zSize);

	// Initialise the manager, grid information and topology
	mpi_init();
}

/// \brief	Default destructor.
///
///			Also closes the MPI logfile.
///
MpiManager::~MpiManager(void)
{
	// Close the logfile
	if (logout != nullptr)
	{
		logout->close();
		delete logout;
	}
}

/// Instance creator
MpiManager* MpiManager::getInstance() {

	if (!me) me = new MpiManager();	// Private construction
	return me;						// Return pointer to new object

}

/// Instance destroyer
void MpiManager::destroyInstance() {

	delete me;			// Delete pointer from static context (destructor will be called automatically)

}

// ************************************************************************* //
/// \brief	Initialisation routine.
///
///			Method is responsible for initialising the MPI topolgy and 
///			associated data. Must be called immediately after MPI_init().
///			For serial vuilds this gets called simply to intialise the 
///			MPIM with a basic set of grid information used by other methods.
void MpiManager::mpi_init()
{

	// Create communicator and topology
	int MPI_periodic[L_DIMS], MPI_reorder;
	MPI_reorder = true;
	dimensions[0] = L_MPI_XCORES; dimensions[1] = L_MPI_YCORES;
	MPI_periodic[0] = true;	MPI_periodic[1] = true;
#if (L_DIMS == 3)
	dimensions[2] = L_MPI_ZCORES; MPI_periodic[2] = true;
#endif

	MPI_Cart_create(MPI_COMM_WORLD, L_DIMS, dimensions, MPI_periodic, MPI_reorder, &world_comm);

	// Get Cartesian topology info
	MPI_Comm_rank( world_comm, &my_rank );
	MPI_Comm_size( world_comm, &num_ranks );

	// Store coordinates in the new topology
	MPI_Cart_coords(world_comm, my_rank, L_DIMS, rank_coords);

	// Output directory creation (only master rank)
	if (my_rank == 0) GridUtils::createOutputDirectory(GridUtils::path_str);

	// Buffer for passing path to other ranks
	char* path_buffer = const_cast<char*>(GridUtils::path_str.c_str());
	int path_buffer_size = static_cast<int>(GridUtils::path_str.size());

	// Broadcast directory name (acquire directory name if not rank 0)
	MPI_Bcast(path_buffer, path_buffer_size, MPI_CHAR, 0, world_comm);
	if (my_rank != 0) {
		std::string char_to_str(path_buffer);
		GridUtils::path_str = char_to_str;
	}

#ifdef L_MPI_VERBOSE
	
	// Open logfile now my_rank has been assigned
	logout->open( GridUtils::path_str + "/mpi_log_rank" + std::to_string(my_rank) + ".out", std::ios::out );

#else

	logout = nullptr;

#endif	


#ifdef L_MPI_VERBOSE
	// State my rank
	L_INFO("My rank is " + std::to_string(my_rank) + ". There are " + std::to_string(num_ranks) + " ranks.", logout);

	// Write out coordinates to application log
	std::string msg("Coordinates on rank " + std::to_string(my_rank) + " are (");
	for (size_t d = 0; d < L_DIMS; d++) {
		msg += "\t" + std::to_string(rank_coords[d]);
	}
	msg += "\t)";
	L_INFO(msg, logout);
#endif


	// Get neighbour ID //

	// Loop over grid direction
	for (int dir = 0; dir < L_MPI_DIRS; dir++) {

		MPI_Barrier(world_comm);

		// Get coordinates of neighbour (taking into account periodic structure)
		int coord_tmp[L_DIMS];
		for (size_t d = 0; d < L_DIMS; d++) {
			neighbour_coords[d][dir] = (rank_coords[d] + neighbour_vectors[d][dir] + dimensions[d]) % dimensions[d];
			coord_tmp[d] = neighbour_coords[d][dir];	// Store single vector for getting neighbour rank
		}

		// Get rank of neighbour and build vector of neighbour ranks
		int tmp;
		MPI_Cart_rank(world_comm, coord_tmp, &tmp);
		neighbour_rank[dir] = tmp;

		MPI_Barrier(world_comm);

#ifdef L_MPI_VERBOSE
		// Print out current neighbour coordinates and rank
		*logout << "Neighbour in direction " << dir << " of rank " << my_rank << " is (";
		for (size_t d = 0; d < L_DIMS; d++) {
			*logout << "\t" << neighbour_coords[d][dir];
		}
		*logout << "\t): Rank " << neighbour_rank[dir] << std::endl;
#endif

	}

	// End Initialisation //

	return;
}

// ************************************************************************* //
/// \brief	Domain decomposition.
///
///			Method to decompose the domain and identify local grid sizes.
///			Parameters defined here are used in GridObj construction.
///			Grid manager must have been initialised before calleing this hence the
///			requirement o have it as a non-null input parameter.
///
///	\param	grid_man	Pointer to an initialised grid manager.
void MpiManager::mpi_gridbuild(GridManager* const grid_man)
{
	
	// Check that when using MPI at least 2 cores have been specified as have assumed so in implementation
	if (
		L_MPI_XCORES < 2 || L_MPI_YCORES < 2
#if (L_DIMS == 3)
		|| L_MPI_ZCORES < 2
#endif
		)
	{
		L_ERROR("When using MPI must use at least 2 cores in each direction. Exiting.", GridUtils::logfile);
	}

	// If using custom sizes, user must set the L_MPI_ZCORES to 1
	if (L_DIMS == 2 && L_MPI_ZCORES != 1) {
		L_ERROR("L_MPI_ZCORES must be set to 1 when using custom MPI sizes in 2D. Exiting.", GridUtils::logfile);
	}
	
	// Global physical dimensions
	double Lx = L_BX;
	double dh = Lx / static_cast<double>(L_N);

	// Fill up the cRankSize arrays with the coordinates of each mpi region. 
	// Auxiliary variables
	int numCells[3];
	int numCores[3];
	numCells[0] = L_N;
	numCells[1] = L_M;
	numCells[2] = L_K;
	numCores[0] = L_MPI_XCORES;
	numCores[1] = L_MPI_YCORES;
#if (L_DIMS == 3)
	numCores[2] = L_MPI_ZCORES;
#else
	numCores[2] = 1;
#endif

	int cellsInDomains[3];
	int cellsLastDomain[3];

	// Set the last position to 0, just in case we are working with L_DIMS = 2
	cellsInDomains[2] = 0;
	cellsLastDomain[2] = 0;

	// Loop over dimensions
	for (size_t d = 0; d < L_DIMS; d++)
	{

		// Calculate the size of each domain
		cellsInDomains[d] = (int)ceil((float)numCells[d] / (float)numCores[d]);

		// Calculate the size of the last domain
		cellsLastDomain[d] = cellsInDomains[d] - (cellsInDomains[d]*numCores[d] - numCells[d]);

		// Throw an error if cellsInLast domain is <=0
		if (cellsLastDomain[d] <= 0)
		{
			L_ERROR("Last core in dir = " + std::to_string(d) + 
				" has 0 or fewer cells. Exiting. Change the number of cores in this direction or adjust resolution of problem.", logout);
		}
	}

	
	// Fill the cRankSize arrays
	int ind = 0;
	for (int i = 0; i < numCores[0]; i++)
	{
		int sX = i == (numCores[0] - 1) ? cellsLastDomain[0] : cellsInDomains[0];

		for (int j = 0; j < numCores[1]; j++)
		{
			int sY = j == (numCores[1] - 1) ? cellsLastDomain[1] : cellsInDomains[1];

			for (int k = 0; k < numCores[2]; k++)
			{
				int sZ = k == (numCores[2] - 1) ? cellsLastDomain[2] : cellsInDomains[2];
				cRankSizeX[ind] = sX;
				cRankSizeY[ind] = sY;
				cRankSizeZ[ind] = sZ;
				ind++;
			}
		}
	}

	
	// Compute required local grid size to pass to grid manager //
	std::vector<int> local_size;

	// Loop over dimensions
	for (size_t d = 0; d < L_DIMS; d++) {

		if (dimensions[d] == 1) {
			// If only 1 rank in this direction local grid is same size a global grid
			local_size.push_back( grid_man->global_size[d][0] );

		} else {

			// Else, get grids sizes from the arrays
			switch (d)
			{
			case 0:
				local_size.push_back( cRankSizeX[my_rank] + 2 );
				break;
			case 1:
				local_size.push_back( cRankSizeY[my_rank] + 2 );
				break;
			case 2:
				local_size.push_back( cRankSizeZ[my_rank] + 2 );
				break;

			default:
				break;
			}
		}
	}

	// Set sizes in grid manager
	grid_man->setLocalCoarseSize(local_size);

	// Resize the rank core edge array
	rank_core_edge.resize(6, std::vector<double>(num_ranks));

	// Find positions of edges of each grid excluding the overlap (recv layer)
	// If using custom sizing need to cumulatively establish how far from origin
	int adj_rank;
	
	// X edges
	rank_core_edge[eXMax][my_rank] = 0.0;
	for (int i = 0; i < rank_coords[0] + 1; i++) {
		// Find i-th rank
#if (L_DIMS == 3)
		int adj_coords[L_DIMS] = {i, rank_coords[1], rank_coords[2]};
#else
		int adj_coords[L_DIMS] = {i, rank_coords[1]};
#endif
		MPI_Cart_rank(world_comm, adj_coords, &adj_rank);
		// Add the width of this rank to total width
		rank_core_edge[eXMax][my_rank] += cRankSizeX[adj_rank] * dh;
	}
	rank_core_edge[eXMin][my_rank] = rank_core_edge[eXMax][my_rank] - cRankSizeX[my_rank] * dh;

	// Y edges
	rank_core_edge[eYMax][my_rank] = 0.0;
	for (int i = 0; i < rank_coords[1] + 1; i++) {
		// Find i-th rank
#if (L_DIMS == 3)
		int adj_coords[L_DIMS] = {rank_coords[0], i, rank_coords[2]};
#else
		int adj_coords[L_DIMS] = {rank_coords[0], i};
#endif
		MPI_Cart_rank(world_comm, adj_coords, &adj_rank);
		// Add the width of this rank to total width
		rank_core_edge[eYMax][my_rank] += cRankSizeY[adj_rank] * dh;
	}
	rank_core_edge[eYMin][my_rank] = rank_core_edge[eYMax][my_rank] - cRankSizeY[my_rank] * dh;

#if (L_DIMS == 3)
	// Z edges
	rank_core_edge[eZMax][my_rank] = 0.0;
	for (int i = 0; i < rank_coords[2] + 1; i++) {
		// Find i-th rank
		int adj_coords[L_DIMS] = {rank_coords[0], rank_coords[1], i};
		MPI_Cart_rank(world_comm, adj_coords, &adj_rank);
		// Add the width of this rank to total width
		rank_core_edge[eZMax][my_rank] += cRankSizeZ[adj_rank] * dh;
	}
	rank_core_edge[eZMin][my_rank] = rank_core_edge[eZMax][my_rank] - cRankSizeZ[my_rank] * dh;
#else
	rank_core_edge[eZMax][my_rank] = L_BZ;
	rank_core_edge[eZMin][my_rank] = 0.0;
#endif
	
	MPI_Barrier(world_comm);

#ifdef L_MPI_VERBOSE
	// Write out the Grid size vector
	*logout << "Grid size including halo on rank " << my_rank << " is (";
	for (size_t d = 0; d < L_DIMS; d++) {
		*logout << "\t" << local_size[d];
	}
	*logout << "\t)" << std::endl;

	*logout << "Limits of the grid core (position) are (" <<
				rank_core_edge[eXMin][my_rank] << "-" << rank_core_edge[eXMax][my_rank] <<
		", " << rank_core_edge[eYMin][my_rank] << "-" << rank_core_edge[eYMax][my_rank] <<
		", " << rank_core_edge[eZMin][my_rank] << "-" << rank_core_edge[eZMax][my_rank] <<
		")" << std:: endl;
#endif

}

// ************************************************************************* //
/// \brief	Buffer ASCII writer.
///
///			When verbose MPI logging is turned on this method will write out 
///			the communication buffer to an ASCII file.
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

// ************************************************************************* //
/// \brief	Communication routine.
///
///			This method implements the communication between grids of the same
///			level and region across MPI processes. Each call effects
///			communication in all valid directions for the grid of the supplied
///			level and region.
///
/// \param	lev	level of grid to communicate.
/// \param	reg	region number of grid to communicate.
void MpiManager::mpi_communicate(int lev, int reg) {

	// Wall clock variables
	clock_t t_start, t_end, secs;

	// Tag
	int TAG;
	int send_count = 0;

	// Get grid object
	GridObj* Grid = NULL;
	GridUtils::getGrid(GridManager::getInstance()->Grids, lev, reg,  Grid);


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
	for (int dir = 0; dir < L_MPI_DIRS; dir++) {

		/* Create a unique tag based on level (< 32), region (< 10) and direction (< 100).
		 * MPICH limits state that tag value cannot be greater than 32767 */
		TAG = ((Grid->level + 1) * 1000) + ((Grid->region_number + 1) * 100) + dir;

#ifdef L_MPI_VERBOSE
		*logout << "Processing Message with Tag --> " << TAG << std::endl;
#endif

		////////////////////////////
		// Resize and Pack Buffer //
		////////////////////////////

		// Adjust buffer size
		for (MpiManager::buffer_struct bufs : buffer_send_info) {
			if (bufs.level == Grid->level && bufs.region == Grid->region_number) {
				f_buffer_send[dir].resize(bufs.size[dir] * L_NUM_VELS);
			}
		}

		// Only pack and send if required
		if (f_buffer_send[dir].size()) {

			// Pass direction and Grid by reference and pack if required
			mpi_buffer_pack( dir, Grid );
		

			///////////////
			// Post Send //
			///////////////

			send_count++;

#ifdef L_MPI_VERBOSE
			*logout << "L" << Grid->level << "R" << Grid->region_number << " -- Direction " << dir 
								<< " -->  Posting Send for " << f_buffer_send[dir].size() / L_NUM_VELS
								<< " sites to Rank " << neighbour_rank[dir] << " with tag " << TAG << "." << std::endl;
#endif
			// Post send message to message queue and log request handle in array
			MPI_Isend( &f_buffer_send[dir].front(), static_cast<int>(f_buffer_send[dir].size()), MPI_DOUBLE, neighbour_rank[dir], 
				TAG, world_comm, &send_requests[send_count-1] );

#ifdef L_MPI_VERBOSE
			*logout << "Direction " << dir << " --> Send Posted." << std::endl;
#endif

		}

		// Find opposite direction (neighbour it receives from)
		int opp_dir = mpi_getOpposite(dir);
		
		// Resize the receive buffer
		for (MpiManager::buffer_struct bufr : buffer_recv_info) {
			if (bufr.level == Grid->level && bufr.region == Grid->region_number) {
				f_buffer_recv[dir].resize(bufr.size[dir] * L_NUM_VELS);
			}
		}


		///////////////////
		// Fetch Message //
		///////////////////

		if (f_buffer_recv[dir].size()) {

#ifdef L_MPI_VERBOSE
			*logout << "L" << Grid->level << "R" << Grid->region_number << " -- Direction " << dir 
								<< " -->  Fetching message for " << f_buffer_recv[dir].size() / L_NUM_VELS	
								<< " sites from Rank " << neighbour_rank[opp_dir] << " with tag " << TAG << "." << std::endl;
#endif

			// Use a blocking receive call if required
			MPI_Recv( &f_buffer_recv[dir].front(), static_cast<int>(f_buffer_recv[dir].size()), MPI_DOUBLE, neighbour_rank[opp_dir], 
				TAG, world_comm, &recv_stat );
		

#ifdef L_MPI_VERBOSE
			*logout << "Direction " << dir << " --> Received." << std::endl;

			// Write out buffers
			*logout << "SUMMARY for L" << Grid->level << "R" << Grid->region_number << " -- Direction " << dir 
								<< " -- Sent " << f_buffer_send[dir].size() / L_NUM_VELS << " to " << neighbour_rank[dir] 
								<< ": Received " << f_buffer_recv[dir].size() / L_NUM_VELS << " from " << neighbour_rank[opp_dir] << std::endl;
			std::string filename = GridUtils::path_str + "/mpiBuffer_Rank" + std::to_string(my_rank) + "_Dir" + std::to_string(dir) + ".out";
			mpi_writeout_buf(filename, dir);
#endif

			///////////////////////////
			// Unpack Buffer to Grid //
			///////////////////////////

			// Pass direction and Grid by reference
			mpi_buffer_unpack( dir, Grid );

		}

	}

#ifdef L_MPI_VERBOSE
	*logout << " *********************** Waiting for Send Completion  *********************** " << std::endl;
#endif

	/* Wait until other processes have handled all the sends from this rank
	 * Note that calls to this command destroy the handles once complete so
	 * do not need to clear the array afterward. */
	MPI_Waitall(send_count,send_requests,send_stat);


	// Print Time of MPI comms
	t_end = clock();
	secs = t_end - t_start;

	// Update average MPI overhead time for this particular grid
	Grid->timeav_mpi_overhead *= (Grid->t-1);
	Grid->timeav_mpi_overhead += ((double)secs)/CLOCKS_PER_SEC;
	Grid->timeav_mpi_overhead /= Grid->t;

#ifdef L_TEXTOUT
	if (Grid->t % L_OUT_EVERY == 0) {
		*GridUtils::logfile << "Writing out to <Grids.out>" << std::endl;
		Grid->io_textout("POST MPI COMMS");
	}
#endif

	if (Grid->t % L_OUT_EVERY == 0) {
		// Performance Data
		*GridUtils::logfile << "MPI overhead taking an average of " << Grid->timeav_mpi_overhead*1000 << "ms" << std::endl;
	}


}

// ************************************************************************* //
/// \brief	Pre-calcualtion of the buffer sizes.
///
///			Wrapper method for computing the buffer sizes for every grid on the
///			rank, both sender and receiver. Must be called post-initialisation.
void MpiManager::mpi_buffer_size() {

	*GridUtils::logfile << "Pre-computing buffer sizes for MPI...";

	/* For each grid in the hierarchy find communicating edges and store the buffer size.
	 * The data are arranged:
	 *
	 *		BufSize[direction][level][region]
	 *
	 * where direction is identified by the MPI labelling at the top of this file.
	 * A zero buffer size indicates that this edge does not communicate with a neighbour rank. */

	// Loop through levels and regions
	GridObj* g;	// Pointer to a GridObj
	for (int l = 0; l <= L_NUM_LEVELS; l++) {
		for (int r = 0; r < L_NUM_REGIONS; r++) {

			// Null the pointer
			g = NULL;

			// Get pointer to grid
			GridUtils::getGrid(GridManager::getInstance()->Grids, l, r, g);

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
#ifdef L_MPI_VERBOSE
		*logout << "Send buffer sizes are [L" << buffer_send_info.back().level << ",R" << buffer_send_info.back().region << "]" << '\t';
		for (int i = 0; i < L_MPI_DIRS; i++) {
			*logout << buffer_send_info.back().size[i] << '\t';
		}
		*logout << std::endl;

		*logout << "Recv buffer sizes are [L" << buffer_recv_info.back().level << ",R" << buffer_recv_info.back().region << "]" << '\t';
		for (int i = 0; i < L_MPI_DIRS; i++) {
			*logout << buffer_recv_info.back().size[i] << '\t';
		}
		*logout << std::endl;
#endif


		}
	}

	*GridUtils::logfile << "Complete." << std::endl;

}

// ************************************************************************* //
/// \brief	Helper method to find opposite direction in MPI topology.
///
///			The MPI directional vectors do not necessarily correspond to the 
///			lattice model direction. The MPI directional vectors are defined 
///			separately and hence there is a separate opposite finding method.
///
/// \param	direction	the outgoing direction whose opposite you wish to find.
int MpiManager::mpi_getOpposite(int direction) {

	/*	If direction is even, then opposite is direction+1.
		If direction is odd, then opposite is direction-1.
		e.g. direction 0 (+x direction) has opposite 1 (-x direction) --> +1
		however, direction 1 has opposite 0 --> -1
		Hence we can add (-1 ^ direction) so it alternates between +/-1
	*/

		return direction + (int)pow(-1,direction);

}
// ************************************************************************* //
/// \brief	Define writable sub-grid communicators.
///
///			When using HDF5 in parallel, collective IO operations require all
///			processes to write a non-zero amount of data to the same file.
///			This method examines availability of sub-grid and writable data on
///			the grid (if found) and ensures it is added to a new communicator. 
///			Must be called AFTER the grids and buffers have been initialised.
///
///	\param	grid_man	pointer to non-null grid manager.
int MpiManager::mpi_buildCommunicators(GridManager* const grid_man) {

	*GridUtils::logfile << "Creating infrastructure for HDF5...";

	// Declarations
	int status = 0;
	int colour;							// Colour indicates which new communicator this process belongs to
	int key = my_rank;					// Global rank as key (dictates numbering in new communicator)

	// Add L0 information (will always have an L0)
	grid_man->createWritableDataStore(grid_man->Grids);

	// Loop over the possible sub-grid combinations and add them to communicator if necessary
	for (int reg = 0; reg < L_NUM_REGIONS; reg++) {
		for (int lev = 1; lev <= L_NUM_LEVELS; lev++) {

#if defined L_HDF_DEBUG
			*logout << "Looking to add sub-grid L" << lev << " R" << reg << " to writable communicator...";
#endif

			// Try gain access to the sub-grid on this rank
			GridObj *targetGrid = NULL;
			GridUtils::getGrid(grid_man->Grids, lev, reg, targetGrid);

			// If sub-grid found
			if (targetGrid != NULL)
			{

#if defined L_HDF_DEBUG
				*logout << "Grid Found." << std::endl;
#endif

				if (grid_man->createWritableDataStore(targetGrid))
				{
					// Writable data found to include in communicator
					colour = 1;					
				}

				else
				{
					// No writable data on the grid so exclude from communicator
					colour = MPI_UNDEFINED;
				}

			}

			// Grid not on this rank
			else
			{
				// Grid not on this rank so exclude from communicator
				colour = MPI_UNDEFINED;

#if defined L_HDF_DEBUG
				*logout << "Grid not found." << std::endl;
#endif
			}


			// Define new communicator which contains writable sub-grids of this level and region
			// by splitting the global communicator using the flags provided.
			status = MPI_Comm_split(world_comm, colour, key, &subGrid_comm[(lev - 1) + reg * L_NUM_LEVELS]);

#if defined L_HDF_DEBUG

			// Only write out communicator info if this rank was added to the communicator
			if (colour != MPI_UNDEFINED)
			{
				*logout << "Grid has writable data of size " << grid_man->p_data.back().writable_data_count <<
					": Region " << reg << ", Level " << lev <<
					" has communicator value: " << subGrid_comm[lev + reg * L_NUM_LEVELS] <<
					" and colour " << colour << std::endl;
			}
#endif

		}
	}

	*GridUtils::logfile << "Communicator build complete. Status = " << status << std::endl;

	return status;
}
// ************************************************************************* //
/// \brief	Update the load balancing information stored in the MpiManager
///
///			This method is executed by all processes. Counts the ACTIVE cells on
///			the current rank and pushes the information to the master (rank 0)
///			which writes this information to an output file if required. Must be 
///			called after the grids have been built or will return zero.
///
///	\param	grid_man	pointer to non-null grid manager.
void MpiManager::mpi_updateLoadInfo(GridManager* const grid_man) {

	/* Loop over the grids and compute the number of ACTIVE cells on the rank.
	 * In other words exclude the cells covered by a sub-grid. However, does 
	 * include halo cells as these are part of the calculation. */
	long active_cell_count = 0;

	for (int lev = 0; lev < L_NUM_LEVELS + 1; ++lev) {
		for (int reg = 0; reg < L_NUM_REGIONS; ++reg) {

			// Get pointer to grid if available on this process
			GridObj *g = NULL;
			GridUtils::getGrid(grid_man->Grids, lev, reg, g);

			if (g != NULL) {

				// Get size (includes halo)
				long grid_cell_count = g->N_lim * g->M_lim * g->K_lim;

#ifdef L_MPI_VERBOSE
				*logout << "Level = " << lev << ", Region = " << reg << ", pointer = " << g << std::endl;
				*logout << "Local Size = " << g->N_lim << "x" << g->M_lim << "x" << g->K_lim << std::endl;
#endif

				// Update total
				if (lev == 0) {

					// Add cells to counter
					active_cell_count += grid_cell_count;
				}

				else {
					// Remove the coarse cells from the count
					active_cell_count -= grid_cell_count / 2;

					// Add the fine cells
					active_cell_count += grid_cell_count;
				}

			}
		}
	}

	// Pass active cell count to master
	long *cell_counts;
	if (my_rank == 0) {

		// Allocate space for receive buffer
		cell_counts = (long*)malloc(num_ranks * sizeof(long));
	}

	// Gather data (into root process)
	MPI_Gather(&active_cell_count, 1, MPI_LONG, cell_counts, 1, MPI_LONG, 0, world_comm);


	// Write information
#ifdef L_MPI_VERBOSE
	*logout << "Active Cell Count (this rank) = " << active_cell_count << std::endl;
#endif

	// If master, write special log file
#ifdef L_MPI_WRITE_LOAD_BALANCE
	
	if (my_rank == 0) {

		std::ofstream counts_out;
		counts_out.open(GridUtils::path_str + "/loadbalancing.out", std::ios::out);

		// Get max load
		double max_load = 0;
		for (int process = 0; process < num_ranks; ++process)
			if (cell_counts[process] > max_load) max_load = static_cast<double>(cell_counts[process]);

		for (int process = 0; process < num_ranks; ++process) {

			// Write process number, active count, load (as percentage of maximum)
			counts_out << process << "\t" << cell_counts[process] << "\t" << 
				(static_cast<double>(cell_counts[process]) * 100.0 / max_load) << std::endl;

		}

		counts_out.close();
	}

#endif

}
// ************************************************************************** //

