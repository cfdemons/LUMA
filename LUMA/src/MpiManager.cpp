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

	// Initialise the manager, grid information and topology
	mpi_init();

	// Resize the IBM-MPI helper classes for each grid level
	markerCommOwnerSide.resize(L_NUM_LEVELS+1);
	markerCommMarkerSide.resize(L_NUM_LEVELS+1);
	supportCommMarkerSide.resize(L_NUM_LEVELS+1);
	supportCommSupportSide.resize(L_NUM_LEVELS+1);
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
///			For serial builds this gets called simply to intialise the 
///			MPIM with a basic set of grid information used by other methods.
void MpiManager::mpi_init()
{

	// Create communicator and topology
	int MPI_periodic[3], MPI_reorder;
	MPI_reorder = true;
	dimensions[0] = L_MPI_XCORES;
	dimensions[1] = L_MPI_YCORES;
	dimensions[2] = L_MPI_ZCORES;
	MPI_periodic[0] = true;
	MPI_periodic[1] = true;
	MPI_periodic[2] = true;

	// Check total number of cores in the topology is the same as
	// the number of ranks in the initial MPI_COMM_WORLD
	// communicator
	int initial_num_ranks = -1;
	int initial_my_rank = -1;
	int total_cores = 1;
	MPI_Comm_size(MPI_COMM_WORLD, &initial_num_ranks);
	MPI_Comm_rank(MPI_COMM_WORLD, &initial_my_rank);
	for (int i = 0; i < L_DIMS; i++) {
	  total_cores *= dimensions[i];
	}
	if (total_cores != initial_num_ranks) {
	  if (initial_my_rank == 0) {
	    std::cerr << "ERROR: The total number of cores (" + std::to_string(total_cores)+") in the MPI topology defined by L_MPI_<d>CORES for a "+
	      std::to_string((int) L_DIMS)+"D simulation in definitions.h is not equal to the number of MPI ranks ("+
	      std::to_string(initial_num_ranks)+") specified in the mpirun command." << std::endl;
	  }
	  abort();
	}

	MPI_Cart_create(MPI_COMM_WORLD, L_DIMS, &dimensions[0], &MPI_periodic[0], MPI_reorder, &world_comm);

	// Get Cartesian topology info
	MPI_Comm_rank(world_comm, &my_rank);
	MPI_Comm_size(world_comm, &num_ranks);

	// Store coordinates in the new topology
	MPI_Cart_coords(world_comm, my_rank, L_DIMS, rank_coords);

	// Output directory creation (only master rank)
	if (my_rank == 0) GridUtils::createOutputDirectory(GridUtils::path_str);

	// Buffer for passing path to other ranks
	char* path_buffer = const_cast<char*>(GridUtils::path_str.c_str());
	int path_buffer_size = static_cast<int>(GridUtils::path_str.size());

	// Broadcast directory name (acquire directory name if not rank 0)
	MPI_Bcast(path_buffer, path_buffer_size, MPI_CHAR, 0, world_comm);
	if (my_rank != 0)
	{
		std::string char_to_str(path_buffer);
		GridUtils::path_str = char_to_str;
	}

#ifdef L_MPI_VERBOSE
	
	// Open logfile now my_rank has been assigned
	logout->open( GridUtils::path_str + "/mpi_log_rank" + std::to_string(my_rank) + ".log", std::ios::out );

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
///			requirement to have it as a non-null input parameter.
///
///	\param	grid_man	Pointer to an initialised grid manager.
void MpiManager::mpi_gridbuild(GridManager* const grid_man)
{
	// Auxiliary variables
	cRankSizeX.resize(num_ranks);
	cRankSizeY.resize(num_ranks);
	cRankSizeZ.resize(num_ranks);
	double dh = L_COARSE_SITE_WIDTH;
	int numCells[3];
	numCells[0] = L_N;
	numCells[1] = L_M;
	numCells[2] = L_K;

	// Compute block sizes based on chosen algorithm
#ifdef L_MPI_TOPOLOGY_REPORT
	mpi_reportOnDecomposition(dh);
#elif defined L_MPI_SMART_DECOMPOSE
	// Log use of SD
	L_INFO("Using Smart Decomposition...", GridUtils::logfile);
	mpi_smartDecompose(dh);
#else
	mpi_uniformDecompose(&numCells[0]);
#endif

#ifdef L_MPI_VERBOSE
	std::string msg("Rank Sizes in the X direction = ");
	for (int i = 0; i < num_ranks; i++) msg += std::to_string(cRankSizeX[i]) + " ";
	L_INFO(msg, logout); msg.clear();
	msg = "Rank Sizes in the Y direction = ";
	for (int i = 0; i < num_ranks; i++) msg += std::to_string(cRankSizeY[i]) + " ";
	L_INFO(msg, logout); msg.clear();
	msg = "Rank Sizes in the Z direction = ";
	for (int i = 0; i < num_ranks; i++) msg += std::to_string(cRankSizeZ[i]) + " ";
	L_INFO(msg, logout); msg.clear();
#endif

	// Compute required local grid size to pass to grid manager //
	std::vector<int> local_size;

	// Loop over dimensions
	for (size_t d = 0; d < L_DIMS; d++)
	{
		if (dimensions[d] == 1)
		{
			// If only 1 rank in this direction local grid is same size a global grid (no halo)
			local_size.push_back(grid_man->global_size[d][0]);

		}
		else
		{
			// Else, get grids sizes from the arrays and add halo
			switch (d)
			{
			case 0:
				local_size.push_back(cRankSizeX[my_rank] + 2);
				break;
			case 1:
				local_size.push_back(cRankSizeY[my_rank] + 2);
				break;
			case 2:
				local_size.push_back(cRankSizeZ[my_rank] + 2);
				break;

			default:
				break;
			}
		}
	}

	// Set sizes in grid manager
	grid_man->setLocalCoarseSize(local_size);


	/* Only rank 0 does the assignment of rank core edges and communicates the
	* result to the other ranks for consistency. */

	// Resize the rank core edge array
	rank_core_edge.resize(6, std::vector<double>(num_ranks));

	if (my_rank == 0)
	{
		// Variables
		int currentRank;
		int prevRank;
		std::vector<int> currentCoords(L_DIMS);
		std::vector<int> prevCoords(L_DIMS);

		// Loop over topology
#if (L_DIMS == 3)
		for (int k = 0; k < dimensions[eZDirection]; k++)
#endif

		{
#if (L_DIMS == 3)
			// Update Z edges
			currentCoords[eZDirection] = k;
			prevCoords[eZDirection] = GridUtils::upToZero(currentCoords[eZDirection] - 1);
#endif

			for (int j = 0; j < dimensions[eYDirection]; j++)
			{

				// Update Y edges
				currentCoords[eYDirection] = j;
				prevCoords[eYDirection] = GridUtils::upToZero(currentCoords[eYDirection] - 1);

				for (int i = 0; i < dimensions[eXDirection]; i++)
				{
					// Increment the coordinates
					currentCoords[eXDirection] = i;
					prevCoords[eXDirection] = GridUtils::upToZero(currentCoords[eXDirection] - 1);

					// Get rank numbers
					MPI_Cart_rank(world_comm, &currentCoords[0], &currentRank);
					MPI_Cart_rank(world_comm, &prevCoords[0], &prevRank);

					// Assign limits where minimum is the same as maximum from previous rank
					if (i == 0) rank_core_edge[eXMin][currentRank] = 0.0;
					else rank_core_edge[eXMin][currentRank] = rank_core_edge[eXMax][prevRank];
					rank_core_edge[eXMax][currentRank] = rank_core_edge[eXMin][currentRank] + (cRankSizeX[currentRank] * dh);

					if (j == 0) rank_core_edge[eYMin][currentRank] = 0.0;
					else rank_core_edge[eYMin][currentRank] = rank_core_edge[eYMax][prevRank];
					rank_core_edge[eYMax][currentRank] = rank_core_edge[eYMin][currentRank] + (cRankSizeY[currentRank] * dh);

#if (L_DIMS == 3)
					if (k == 0) rank_core_edge[eZMin][currentRank] = 0.0;
					else rank_core_edge[eZMin][currentRank] = rank_core_edge[eZMax][prevRank];
					rank_core_edge[eZMax][currentRank] = rank_core_edge[eZMin][currentRank] + (cRankSizeZ[currentRank] * dh);
#else
					rank_core_edge[eZMin][currentRank] = 0.0;
					rank_core_edge[eZMax][currentRank] = grid_man->global_edges[eZMax][0];
#endif

				}
			}
		}

	}

	// Communicate all the grid edges around the topology
	mpi_communicateBlockEdges();

	// Wait for all ranks
	MPI_Barrier(world_comm);

#ifdef L_MPI_VERBOSE
	// Write out the Grid size vector
	*logout << "Grid size including halo on rank " << my_rank << " is (";
	for (size_t d = 0; d < L_DIMS; d++) {
		*logout << "\t" << local_size[d];
	}
	*logout << "\t)" << std::endl;

	L_INFO("Limits of the grid core (position) are ("
		+ std::to_string(rank_core_edge[eXMin][my_rank]) + "-" + std::to_string(rank_core_edge[eXMax][my_rank]) + ", "
		+ std::to_string(rank_core_edge[eYMin][my_rank]) + "-" + std::to_string(rank_core_edge[eYMax][my_rank]) + ", "
		+ std::to_string(rank_core_edge[eZMin][my_rank]) + "-" + std::to_string(rank_core_edge[eZMax][my_rank]) +
		")", logout);
#endif

	/* Now update the halo regions taking into account periodicity.
	 *
	 * Note that if we don't want to have a halo in a given direction due to there
	 * only being one core in that direction and need to avoid Issue #157 then we
	 * need to set the sender and receiver layers to be outside the world space
	 * used by LUMA, i.e. < 0.0. */

	// X
	if (dimensions[eXDirection] == 1)
	{
		sender_layer_pos.X[eLeftMin] = -1.0;
		sender_layer_pos.X[eLeftMax] = -1.0;
		sender_layer_pos.X[eRightMin] = -1.0;
		sender_layer_pos.X[eRightMax] = -1.0;
		recv_layer_pos.X[eLeftMin] = -1.0;
		recv_layer_pos.X[eLeftMax] = -1.0;
		recv_layer_pos.X[eRightMin] = -1.0;
		recv_layer_pos.X[eRightMax] = -1.0;
	}
	else
	{
		sender_layer_pos.X[eLeftMin] = rank_core_edge[eXMin][my_rank];
		sender_layer_pos.X[eLeftMax] = rank_core_edge[eXMin][my_rank] + dh;
		sender_layer_pos.X[eRightMin] = rank_core_edge[eXMax][my_rank] - dh;
		sender_layer_pos.X[eRightMax] = rank_core_edge[eXMax][my_rank];

		// Check if need to wrap or not
		if (sender_layer_pos.X[eLeftMin] - dh < 0.0)
		{
			recv_layer_pos.X[eLeftMin] = grid_man->global_edges[eXMax][0] - dh;
			recv_layer_pos.X[eLeftMax] = grid_man->global_edges[eXMax][0];
		}
		else
		{
			recv_layer_pos.X[eLeftMin] = sender_layer_pos.X[eLeftMin] - dh;
			recv_layer_pos.X[eLeftMax] = sender_layer_pos.X[eLeftMin];
		}
		if (sender_layer_pos.X[eRightMax] + dh > grid_man->global_edges[eXMax][0])
		{
			recv_layer_pos.X[eRightMin] = 0.0;
			recv_layer_pos.X[eRightMax] = dh;
		}
		else
		{
			recv_layer_pos.X[eRightMin] = sender_layer_pos.X[eRightMax];
			recv_layer_pos.X[eRightMax] = sender_layer_pos.X[eRightMax] + dh;
		}
	}

	// Y
	if (dimensions[eYDirection] == 1)
	{
		sender_layer_pos.Y[eLeftMin] = -1.0;
		sender_layer_pos.Y[eLeftMax] = -1.0;
		sender_layer_pos.Y[eRightMin] = -1.0;
		sender_layer_pos.Y[eRightMax] = -1.0;
		recv_layer_pos.Y[eLeftMin] = -1.0;
		recv_layer_pos.Y[eLeftMax] = -1.0;
		recv_layer_pos.Y[eRightMin] = -1.0;
		recv_layer_pos.Y[eRightMax] = -1.0;
	}
	else
	{
		sender_layer_pos.Y[eLeftMin] = rank_core_edge[eYMin][my_rank];
		sender_layer_pos.Y[eLeftMax] = rank_core_edge[eYMin][my_rank] + dh;
		sender_layer_pos.Y[eRightMin] = rank_core_edge[eYMax][my_rank] - dh;
		sender_layer_pos.Y[eRightMax] = rank_core_edge[eYMax][my_rank];
		if (sender_layer_pos.Y[eLeftMin] - dh < 0.0)
		{
			recv_layer_pos.Y[eLeftMin] = grid_man->global_edges[eYMax][0] - dh;
			recv_layer_pos.Y[eLeftMax] = grid_man->global_edges[eYMax][0];
		}
		else
		{
			recv_layer_pos.Y[eLeftMin] = sender_layer_pos.Y[eLeftMin] - dh;
			recv_layer_pos.Y[eLeftMax] = sender_layer_pos.Y[eLeftMin];
		}
		if (sender_layer_pos.Y[eRightMax] + dh > grid_man->global_edges[eYMax][0])
		{
			recv_layer_pos.Y[eRightMin] = 0.0;
			recv_layer_pos.Y[eRightMax] = dh;
		}
		else
		{
			recv_layer_pos.Y[eRightMin] = sender_layer_pos.Y[eRightMax];
			recv_layer_pos.Y[eRightMax] = sender_layer_pos.Y[eRightMax] + dh;
		}
	}

	// Z
#if (L_DIMS == 3)
	if (dimensions[eZDirection] == 1)
	{
		sender_layer_pos.Z[eLeftMin] = -1.0;
		sender_layer_pos.Z[eLeftMax] = -1.0;
		sender_layer_pos.Z[eRightMin] = -1.0;
		sender_layer_pos.Z[eRightMax] = -1.0;
		recv_layer_pos.Z[eLeftMin] = -1.0;
		recv_layer_pos.Z[eLeftMax] = -1.0;
		recv_layer_pos.Z[eRightMin] = -1.0;
		recv_layer_pos.Z[eRightMax] = -1.0;
	}
	else
	{
		sender_layer_pos.Z[eLeftMin] = rank_core_edge[eZMin][my_rank];
		sender_layer_pos.Z[eLeftMax] = rank_core_edge[eZMin][my_rank] + dh;
		sender_layer_pos.Z[eRightMin] = rank_core_edge[eZMax][my_rank] - dh;
		sender_layer_pos.Z[eRightMax] = rank_core_edge[eZMax][my_rank];
		if (sender_layer_pos.Z[eLeftMin] - dh < 0.0)
		{
			recv_layer_pos.Z[eLeftMin] = grid_man->global_edges[eZMax][0] - dh;
			recv_layer_pos.Z[eLeftMax] = grid_man->global_edges[eZMax][0];
		}
		else
		{
			recv_layer_pos.Z[eLeftMin] = sender_layer_pos.Z[eLeftMin] - dh;
			recv_layer_pos.Z[eLeftMax] = sender_layer_pos.Z[eLeftMin];
		}
		if (sender_layer_pos.Z[eRightMax] + dh > grid_man->global_edges[eZMax][0])
		{
			recv_layer_pos.Z[eRightMin] = 0.0;
			recv_layer_pos.Z[eRightMax] = dh;
		}
		else
		{
			recv_layer_pos.Z[eRightMin] = sender_layer_pos.Z[eRightMax];
			recv_layer_pos.Z[eRightMax] = sender_layer_pos.Z[eRightMax] + dh;
		}
}
#endif

#ifdef L_MPI_VERBOSE
	if (dimensions[eXDirection] == 1 || dimensions[eYDirection] == 1
#if (L_DIMS == 3)
		|| dimensions[eZDirection] == 1
#endif
		)
		L_WARN("One of the MPI dimensions is 1 so sender and receiver layers in this direction will be set to -1 as they have no meaning in this context.", logout);

	L_INFO("X sender layers are: " +
		std::to_string(sender_layer_pos.X[eLeftMin]) + " -- " + std::to_string(sender_layer_pos.X[eLeftMax]) + " (min edge) , " +
		std::to_string(sender_layer_pos.X[eRightMin]) + " -- " + std::to_string(sender_layer_pos.X[eRightMax]) + " (max edge)",
		logout);
	L_INFO("X recv layers are: " +
		std::to_string(recv_layer_pos.X[eLeftMin]) + " -- " + std::to_string(recv_layer_pos.X[eLeftMax]) + " (min edge) , " +
		std::to_string(recv_layer_pos.X[eRightMin]) + " -- " + std::to_string(recv_layer_pos.X[eRightMax]) + " (max edge)",
		logout);

	L_INFO("Y sender layers are: " +
		std::to_string(sender_layer_pos.Y[eLeftMin]) + " -- " + std::to_string(sender_layer_pos.Y[eLeftMax]) + " (min edge) , " +
		std::to_string(sender_layer_pos.Y[eRightMin]) + " -- " + std::to_string(sender_layer_pos.Y[eRightMax]) + " (max edge)",
		logout);
	L_INFO("Y recv layers are: " +
		std::to_string(recv_layer_pos.Y[eLeftMin]) + " -- " + std::to_string(recv_layer_pos.Y[eLeftMax]) + " (min edge) , " +
		std::to_string(recv_layer_pos.Y[eRightMin]) + " -- " + std::to_string(recv_layer_pos.Y[eRightMax]) + " (max edge)",
		logout);

#if (L_DIMS == 3)
	L_INFO("Z sender layers are: " +
		std::to_string(sender_layer_pos.Z[eLeftMin]) + " -- " + std::to_string(sender_layer_pos.Z[eLeftMax]) + " (min edge) , " +
		std::to_string(sender_layer_pos.Z[eRightMin]) + " -- " + std::to_string(sender_layer_pos.Z[eRightMax]) + " (max edge)",
		logout);
	L_INFO("Z recv layers are: " +
		std::to_string(recv_layer_pos.Z[eLeftMin]) + " -- " + std::to_string(recv_layer_pos.Z[eLeftMax]) + " (min edge) , " +
		std::to_string(recv_layer_pos.Z[eRightMin]) + " -- " + std::to_string(recv_layer_pos.Z[eRightMax]) + " (max edge)",
		logout);
#endif
#endif

}


// ************************************************************************* //
/// \brief	Communication of all rank core edges around the topology.
///
///			After this all ranks will have positional limits of all rank cores.
void MpiManager::mpi_communicateBlockEdges()
{

	// Declare exchange buffer
	int numElements = static_cast<int>(rank_core_edge.size() * rank_core_edge[0].size());
	std::vector<double> buffer(numElements, 0.0);

	if (my_rank == 0)
	{
		// Populate buffer
		for (int edge = 0; edge < rank_core_edge.size(); edge++)
		{
			for (int rank = 0; rank < num_ranks; rank++)
			{
				buffer[edge + rank_core_edge.size() * rank] = rank_core_edge[edge][rank];
			}
		}
	}

	// Broadcast to other ranks
	MPI_Bcast(&buffer[0], numElements, MPI_DOUBLE, 0, world_comm);

	// Unpack the buffer
	for (int edge = 0; edge < rank_core_edge.size(); edge++)
	{
		for (int rank = 0; rank < num_ranks; rank++)
		{
			rank_core_edge[edge][rank] = buffer[edge + rank_core_edge.size() * rank];
		}
	}
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
	for (int dir = 0; dir < L_MPI_DIRS; dir++)
	{

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
		for (MpiManager::BufferSizeStruct bufs : buffer_send_info) {
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
		for (MpiManager::BufferSizeStruct bufr : buffer_recv_info) {
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
#endif

			///////////////////////////
			// Unpack Buffer to Grid //
			///////////////////////////

			// Pass direction and Grid by reference
			mpi_buffer_unpack( dir, Grid );

		}

#ifdef L_MPI_VERBOSE

		*logout << "SUMMARY for L" << Grid->level << "R" << Grid->region_number << " -- Direction " << dir
			<< " -- Sent " << f_buffer_send[dir].size() / L_NUM_VELS << " to " << neighbour_rank[dir]
			<< ": Received " << f_buffer_recv[dir].size() / L_NUM_VELS << " from " << neighbour_rank[opp_dir] << std::endl;

		// Write out buffers
		std::string filename = GridUtils::path_str + "/mpiBuffer_Rank" + std::to_string(my_rank) + "_Dir" + std::to_string(dir) + ".out";
		mpi_writeout_buf(filename, dir);
#endif

	}

#ifdef L_MPI_VERBOSE
	*logout << " *********************** Waiting for Sends to be Received on L" + 
		std::to_string(lev) + "R" + std::to_string(reg) + 
		" *********************** " << std::endl;
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
	if (Grid->t % L_GRID_OUT_FREQ == 0) {
		*GridUtils::logfile << "Writing out to <Grids.out>" << std::endl;
		Grid->io_textout("POST MPI COMMS");
	}
#endif

	if (Grid->t % L_GRID_OUT_FREQ == 0) {
		// Performance Data
		L_INFO("MPI overhead taking an average of " + 
			std::to_string(Grid->timeav_mpi_overhead * 1000) + "ms", GridUtils::logfile);
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
			if (g == NULL)
			{
				continue;	// Try find next grid in loop
			}

			// Expand buffer info arrays by one and add grid ID info
			buffer_send_info.emplace_back(l, r);
			buffer_recv_info.emplace_back(l, r);

			/* Not every edge of a grid needs to be communicated.
			 * We therefore loop over the appropriate region of the grid and 
			 * check whether a site lies in the associated sender / receiver 
			 * layer. If so increment the buffer size counter. */

			// Call send and recv size finding routines
			mpi_buffer_size_send(g);
			mpi_buffer_size_recv(g);


#ifdef L_MPI_VERBOSE
			// Write out buffer sizes for reference
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

#ifdef L_MPI_VERBOSE
	/* Historically, there have been cases of MPI hangs due to buffer
	* size inconsistencies between ranks. These are difficult to debug
	* as it involves going through the log files and cross checking
	* each send and receive by hand. To speed up the process we
	* implement a debugging activity here to check the send and receive
	* buffer sizes tally up across the topology. */
	int * const bufSendSizesOut = new int[L_NUM_LEVELS * L_NUM_REGIONS + 1];
	int * const bufSendSizesIn = new int[L_NUM_LEVELS * L_NUM_REGIONS + 1];
	bool bInfoFound = false;

	// Loop over the MPI directions
	for (int i = 0; i < L_MPI_DIRS; i++)
	{
		int counter = 0;

		// Loop over grids and sub-grids and pack buffer one at a time
		for (int l = 0; l <= L_NUM_LEVELS; l++)
		{
			for (int r = 0; r < L_NUM_REGIONS; r++)
			{
				bInfoFound = false;
				if (l == 0 && r != 0) continue;		// L0 can only be R0

				// Try retireve the buffer size info
				for (MpiManager::BufferSizeStruct bufs : buffer_send_info)
				{
					if (bufs.level == l && bufs.region == r)
					{
						bufSendSizesOut[counter] = bufs.size[i];
						bInfoFound = true;
					}
				}

				// If no info found then grid doesn't exist on this rank so mark as zero
				if (!bInfoFound) bufSendSizesOut[counter] = 0;

				// Increment counter
				counter++;
			}
		}

		L_INFO("Send sizes packed for direction " + std::to_string(i) + ". Sending to rank " + std::to_string(neighbour_rank[i]) + "...", logout);

		// Once buffer complete, post send
		MPI_Isend(bufSendSizesOut, L_NUM_LEVELS * L_NUM_REGIONS + 1, MPI_INT, neighbour_rank[i],
			my_rank, world_comm, &send_requests[i]);

		L_INFO("Sent.", logout);

		// Get opposite direction
		int opp_i = mpi_getOpposite(i);
		L_INFO("Receiving send sizes from rank " + std::to_string(neighbour_rank[opp_i]) + " in direction " + std::to_string(opp_i) + "...", logout);

		// Receive the send sizes from opposite direction
		MPI_Recv(bufSendSizesIn, L_NUM_LEVELS * L_NUM_REGIONS + 1, MPI_INT, neighbour_rank[opp_i],
			neighbour_rank[opp_i], world_comm, &recv_stat);

		L_INFO("Received.", logout);

		// Cross check that the send sizes match the receive sizes
		counter = 0;
		// Loop over grids and sub-grids and unpack buffer one at a time
		for (int l = 0; l <= L_NUM_LEVELS; l++)
		{
			for (int r = 0; r < L_NUM_REGIONS; r++)
			{
				bInfoFound = false;
				if (l == 0 && r != 0) continue;		// L0 can only be R0

				// Try retireve the buffer size info
				for (MpiManager::BufferSizeStruct bufr : buffer_recv_info)
				{
					if (bufr.level == l && bufr.region == r)
					{
						bInfoFound = true;
						if (bufSendSizesIn[counter] != bufr.size[i])
						{
							L_ERROR("L" + std::to_string(l) + "R" + std::to_string(r) + 
								" -- Mismatch in buffer sizes! Rank " + std::to_string(neighbour_rank[opp_i]) +
								" is sending " + std::to_string(bufSendSizesIn[counter]) + 
								". Expected to receive " + std::to_string(bufr.size[i]), GridUtils::logfile);
						}						
					}
				}

				// If no info found then grid doesn't exist on this rank -- make sure buffer size is zero
				if (!bInfoFound)
				{
					if (bufSendSizesIn[counter] != 0)
					{
						L_ERROR("L" + std::to_string(l) + "R" + std::to_string(r) + 
							" -- Mismatch in buffer sizes! Rank " + std::to_string(neighbour_rank[opp_i]) +
							" is sending " + std::to_string(bufSendSizesIn[counter]) +
							". Expected to receive 0 as no grid of this L+R on this rank.", GridUtils::logfile);
					}
				}

				// Increment counter
				counter++;
			}
		}		
	}
	delete[] bufSendSizesOut;
	delete[] bufSendSizesIn;

#endif

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
	
	L_INFO("Creating infrastructure for HDF5...", GridUtils::logfile);
	
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
			L_INFO("MpiManager: Looking to add sub-grid L" + std::to_string(lev) +
				" R" + std::to_string(reg) + " to writable communicator...", GridUtils::logfile);
#endif

			// Try gain access to the sub-grid on this rank
			GridObj *targetGrid = NULL;
			GridUtils::getGrid(grid_man->Grids, lev, reg, targetGrid);

			// If sub-grid found
			if (targetGrid != NULL)
			{

#if defined L_HDF_DEBUG
				L_INFO("Grid Found.", GridUtils::logfile);
#endif

				if (grid_man->createWritableDataStore(targetGrid))
				{
					// Writable data found to include in communicator
					colour = lev + reg * L_NUM_LEVELS;
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
				L_INFO("Grid NOT found.", GridUtils::logfile);
#endif
			}


			// Define new communicator which contains writable sub-grids of this level and region
			// by splitting the global communicator using the flags provided.
			status = MPI_Comm_split(world_comm, colour, key, &subGrid_comm[(lev - 1) + reg * L_NUM_LEVELS]);
			if (status != MPI_SUCCESS) L_ERROR("Sub-grid comm split was unsuccessful.", GridUtils::logfile);

#if defined L_HDF_DEBUG

			// Only write out communicator info if this rank was added to the communicator
			if (colour != MPI_UNDEFINED)
			{
				L_INFO("Grid has writable data of size " + std::to_string(grid_man->p_data.back().writable_data_count)
					+ ": Region " + std::to_string(reg) + ", Level " + std::to_string(lev)
					+ " has communicator value: " + std::to_string(subGrid_comm[(lev - 1) + reg * L_NUM_LEVELS])
					+ " and colour " + std::to_string(colour), GridUtils::logfile);
			}
			else
			{
				L_INFO("Grid has no writable data", GridUtils::logfile);
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
///			This method is executed by all processes. Counts the ACTIVE cell 
///			operations on the current rank and pushes the to the master (rank 0)
///			which writes this information to an output file if required. Must be 
///			called after the grids have been built or will return zero.
///			Active cell count takes into account the number of operations 
///			each cell requires per timestep by scaling the count by the 
///			refinement ratio. i.e. in 2D a coarse site refined to level 1 is 
///			represented	by 2^2 cells which each do 2^1 operations per coarse 
///			time step. Operation count is therefore 2^3 = 8 operations.
///
///	\param	grid_man	pointer to non-null grid manager.
void MpiManager::mpi_updateLoadInfo(GridManager* const grid_man) {

	/* Loop over the grids and compute the number of ACTIVE cells on the rank.
	 * In other words exclude the cells covered by a sub-grid. However, does 
	 * include halo cells as these are part of the calculation. Count up ops
	 * per coarse time step. */
	long ops_count = 0;

	for (int lev = 0; lev < L_NUM_LEVELS + 1; ++lev)
	{
		for (int reg = 0; reg < L_NUM_REGIONS; ++reg)
		{

			// Get pointer to grid if available on this process
			GridObj *g = nullptr;
			GridUtils::getGrid(lev, reg, g);

			if (g)
			{

				// Get size (includes halo)
				long grid_cell_count = g->N_lim * g->M_lim * g->K_lim;

#ifdef L_MPI_VERBOSE
				*logout << "Level = " << lev << ", Region = " << reg << ", pointer = " << g << std::endl;
				*logout << "Local Size = " << g->N_lim << "x" << g->M_lim << "x" << g->K_lim << std::endl;
#endif

				// Update total operations (cells * sub-cycles)
				ops_count += grid_cell_count * static_cast<long>(pow(2, lev));

				if (lev != 0)
				{
					// Remove the coarse ops from the count
					ops_count -= grid_cell_count * static_cast<long>(pow(2, lev - 1 - L_DIMS));
				}

			}
		}
	}

	// Gather ops count (into root process)
	long *ops_count_buffer = nullptr;
	if (my_rank == 0)
	{
		// Allocate space for receive buffer
		ops_count_buffer = new long[num_ranks];
		MPI_Gather(&ops_count, 1, MPI_LONG, ops_count_buffer, 1, MPI_LONG, 0, world_comm);
	}
	else
	{
		// No buffer need for non-root processes
		MPI_Gather(&ops_count, 1, MPI_LONG, NULL, 1, MPI_LONG, 0, world_comm);
	}


	// Write information
#ifdef L_MPI_VERBOSE
	*logout << "Lattice updates on this rank = " << ops_count << std::endl;
#endif

	// If master, write special log file
#ifdef L_MPI_WRITE_LOAD_BALANCE
	
	if (my_rank == 0)
	{

		std::ofstream counts_out;
		counts_out.open(GridUtils::path_str + "/loadbalancing.out", std::ios::out);

		// Get max load
		double max_load = 0;
		for (int process = 0; process < num_ranks; ++process)
			if (ops_count_buffer[process] > max_load) max_load = static_cast<double>(ops_count_buffer[process]);

		for (int process = 0; process < num_ranks; ++process)
		{
			// Write process number, active count, load (as percentage of maximum)
			counts_out << process << "\t" << ops_count_buffer[process] << "\t" << 
				(static_cast<double>(ops_count_buffer[process]) * 100.0 / max_load) << std::endl;
		}

		counts_out.close();
	}

#endif

	if (ops_count_buffer) delete[] ops_count_buffer;

}
// ************************************************************************* //
/// \brief	Populate the rank size arrays based on uniform decomposition logic
void MpiManager::mpi_uniformDecompose(int *numCells)
{
	// Auxiliary variables
	int cellsInDomains[3];
	int cellsLastDomain[3];

	// Set the last position to 0, just in case we are working with L_DIMS = 2
	cellsInDomains[2] = 0;
	cellsLastDomain[2] = 0;

	// Compute domain size based on uniform decomposition
	for (size_t d = 0; d < L_DIMS; d++)
	{
		// Calculate the size of each domain in dimension d
		cellsInDomains[d] = static_cast<int>(std::ceil(static_cast<double>(numCells[d]) / static_cast<double>(dimensions[d])));

		// Calculate the size of the last domain
		cellsLastDomain[d] = cellsInDomains[d] - (cellsInDomains[d] * dimensions[d] - numCells[d]);

		// If cellsInLast domain is <=0
		if (cellsLastDomain[d] <= 0)
		{
			// Try using floor instead of ceil
			cellsInDomains[d] = static_cast<int>(std::floor(static_cast<double>(numCells[d]) / static_cast<double>(dimensions[d])));
			cellsLastDomain[d] = cellsInDomains[d] - (cellsInDomains[d] * dimensions[d] - numCells[d]);

			// If still na issue then error
			if (cellsLastDomain[d] <= 0)
			{
				L_ERROR("Last core in dir = " + std::to_string(d) +
					" has 0 or fewer cells. Exiting. Change the number of cores in this direction or adjust resolution of problem.", logout);
			}
		}
	}


	// Fill the cRankSize arrays for each MPI block
	int ind = 0;
	for (int i = 0; i < dimensions[eXDirection]; i++)
	{
		int sX = i == (dimensions[eXDirection] - 1) ? cellsLastDomain[eXDirection] : cellsInDomains[eXDirection];

		for (int j = 0; j < dimensions[eYDirection]; j++)
		{
			int sY = j == (dimensions[eYDirection] - 1) ? cellsLastDomain[eYDirection] : cellsInDomains[eYDirection];

			for (int k = 0; k < dimensions[eZDirection]; k++)
			{
				int sZ = k == (dimensions[eZDirection] - 1) ? cellsLastDomain[eZDirection] : cellsInDomains[eZDirection];
				cRankSizeX[ind] = sX;
				cRankSizeY[ind] = sY;
				cRankSizeZ[ind] = sZ;
				ind++;
			}
		}
	}
}

// ************************************************************************* //
/// \brief	Checks whether the perturbation vector is valid and adjusts delta 
///			if necessary.
///
///			Also reconstructs the solution vector at the same time.
///
///	\param		solutionData	structure to hold SD information.
///	\param		dh				coarsest cell width.
///	\param		numCores		reference to vector holding core topology.
///	\returns					indication whether the perturbation was valid or not.
bool MpiManager::mpi_SDCheckDelta(SDData& solutionData, double dh, std::vector<int>& numCores)
{

	// Declarations
	bool bAdjustComplete = false;
	bool bAdjusted = false;
	unsigned int c = 0;

	// Create perturbed variable vector using current delta
	for (size_t p = 0; p < solutionData.theta.size(); ++p)
		solutionData.thetaNew[p] = solutionData.theta[p] + solutionData.delta[p];

	// Reconstruct the solution vector
	mpi_SDReconstructSolution(solutionData, numCores);

	// Loop
	while (!bAdjustComplete)
	{
		// Reset the adjusted flag
		bAdjusted = false;		

		// Check validity of edges and reset the delta if causes a problem
		for (int i = 0; i < numCores[eXDirection]; ++i)
		{
			for (int j = 0; j < numCores[eYDirection]; ++j)
			{
				for (int k = 0; k < numCores[eZDirection]; ++k)
				{
					// Nothing to change if on the last block
					if (
						i == numCores[eXDirection] - 1 || j == numCores[eYDirection] - 1
#if (L_DIMS == 3)
						|| k == numCores[eZDirection] - 1
#endif
						) continue;

					// Perturbation cannot cause a block to have a zero size
					if ((solutionData.XSol[i + 1] - solutionData.XSol[i]) < dh)
					{
						// Adjust delta such that size is equal to dh
						c = i;
						solutionData.delta[c] = solutionData.XSol[i] - solutionData.theta[c] + dh;
						solutionData.thetaNew[c] = solutionData.theta[c] + solutionData.delta[c];
						solutionData.XSol[i + 1] = solutionData.thetaNew[c];
						bAdjusted = true;
					}
					if ((solutionData.YSol[j + 1] - solutionData.YSol[j]) < dh)
					{
						c = numCores[eXDirection] - 1 + j;
						solutionData.delta[c] = solutionData.YSol[j] - solutionData.theta[c] + dh;
						solutionData.thetaNew[c] = solutionData.theta[c] + solutionData.delta[c];
						solutionData.YSol[j + 1] = solutionData.thetaNew[c];
						bAdjusted = true;
					}
#if (L_DIMS == 3)
					if ((solutionData.ZSol[k + 1] - solutionData.ZSol[k]) < dh)
					{
						c = numCores[eXDirection] + numCores[eYDirection] - 2 + k;
						solutionData.delta[c] = solutionData.ZSol[k] - solutionData.theta[c] + dh;
						solutionData.thetaNew[c] = solutionData.theta[c] + solutionData.delta[c];
						solutionData.ZSol[k + 1] = solutionData.thetaNew[c];
						bAdjusted = true;
					}
#endif
				}
			}
		}
		bAdjustComplete = true;
	}

	return !bAdjusted;

}

// ************************************************************************* //
/// \brief	Reassembles the solutions vectors given the vector of variable block edges.
///
///	\param	solutionData	structure to hold SD information.
///	\param	numCores		reference to vector holding core topology.
void MpiManager::mpi_SDReconstructSolution(SDData& solutionData, std::vector<int>& numCores)
{
	// 1D index for theta
	int c = 0;
	int i = 0;
	for (i = 1; i < numCores[eXDirection]; ++i)
	{
		solutionData.XSol[i] = solutionData.theta[c];
		c++;
	}
	for (i = 1; i < numCores[eYDirection]; ++i)
	{
		solutionData.YSol[i] = solutionData.theta[c];
		c++;
	}
	for (i = 1; i < numCores[eZDirection]; ++i)
	{
		solutionData.ZSol[i] = solutionData.theta[c];
		c++;
	}

}

// ************************************************************************* //
/// \brief	Populate the active cell counts for the given block configuration.
///
///			Imbalance measured as the difference between the heaviest and 
///			lightest block as a percentage of the heaviest. The load imbalance
///			information is used to update the structure provided.
///
///	\param[out]	load			load imbalance information structure.
///	\param		solutionData	structure to hold SD information.
///	\param		numCores		reference to vector holding core topology.					.
void MpiManager::mpi_SDComputeImbalance(LoadImbalanceData& load,
	SDData& solutionData, std::vector<int>& numCores)
{
	size_t count = 0;
	size_t countMax = 0;
	size_t countMin = std::numeric_limits<size_t>::max();

	// Construct bounds for each block and then find active cell count from grid manager
	double bounds[6];
	for (int i = 0; i < numCores[eXDirection]; ++i)
	{
		for (int j = 0; j < numCores[eYDirection]; ++j)
		{
			for (int k = 0; k < numCores[eZDirection]; ++k)
			{
				// Set bounds from solution vector
				bounds[eXMin] = solutionData.XSol[i];
				bounds[eXMax] = solutionData.XSol[i + 1];
				bounds[eYMin] = solutionData.YSol[j];
				bounds[eYMax] = solutionData.YSol[j + 1];
				bounds[eZMin] = solutionData.ZSol[k];
				bounds[eZMax] = solutionData.ZSol[k + 1];

				// Get active operation count
				count = GridManager::getInstance()->getActiveCellCount(&bounds[0], true);

				// Update the extremes
				if (count > countMax)
				{
					countMax = count;
					load.heaviestBlock[eXDirection] = i;
					load.heaviestBlock[eYDirection] = j;
					load.heaviestBlock[eZDirection] = k;
				}
				if (count < countMin) countMin = count;

			}
		}
	}

	// Update load imbalance
	load.loadImbalance = 
		std::abs(static_cast<double>(countMax) - static_cast<double>(countMin)) * 100.0 / static_cast<double>(countMax);
	load.heaviestOps = countMax;

}

// ************************************************************************* //
/// \brief	Populate the rank size arrays based on an algorithm that seeks to 
///			load balance.
///
///			This method is independent of the topology in use. Topologies of custom
///			dimensions can be passed through the optional argument.
///
///	\param	reqDims		pointer to a vector containing the desired MPI dimensions.
///	\param	dh			size of a voxel on the coarsest grid.
///	\returns			structure containing imbalance information.
MpiManager::LoadImbalanceData MpiManager::mpi_smartDecompose(double dh, std::vector<int> reqDims)
{
	// Make a suitable copy of the information depending on the argument passed in
	std::vector<int> numCores(3);
	if (!reqDims.size())
	{
		numCores[eXDirection] = L_MPI_XCORES;
		numCores[eYDirection] = L_MPI_YCORES;
		numCores[eZDirection] = L_MPI_ZCORES;
	}
	else
	{
		numCores = reqDims;
	}

	// Set up solution vector for variables
	SDData solutionData;
	solutionData.XSol.resize(numCores[eXDirection] + 1, -1 * dh);
	solutionData.YSol.resize(numCores[eYDirection] + 1, -1 * dh);
	solutionData.ZSol.resize(numCores[eZDirection] + 1, -1 * dh);

	// Create imbalance structure
	LoadImbalanceData load;

	// Only rank 0 does the calculation and result is pushed to other ranks
	if (my_rank == 0)
	{
		// Data
		int p = (numCores[eXDirection] + numCores[eYDirection] + numCores[eZDirection]) - 3;	// Number of unknowns
		int i = 0;
		int c = 0;
		std::vector<int> domainSize(3);
		domainSize[eXDirection] = L_N;
		domainSize[eYDirection] = L_M;
		domainSize[eZDirection] = L_K;

		// Fix the edges as we know where they are
		solutionData.XSol[0] = 0.0;
		solutionData.YSol[0] = 0.0;
		solutionData.ZSol[0] = 0.0;
		solutionData.XSol.back() = GridManager::getInstance()->global_edges[eXMax][0];
		solutionData.YSol.back() = GridManager::getInstance()->global_edges[eYMax][0];
		solutionData.ZSol.back() = GridManager::getInstance()->global_edges[eZMax][0];

		// Handle the 1, 1, 1 case
		if (p == 0)
		{
			// Communicate information around topology if not performing a report
			if (!reqDims.size()) mpi_SDCommunicateSolution(solutionData, load.loadImbalance, dh);

			// Return as no need to perform the iteration
			return load;
		}

		// Create theta vector with initial guesses (uniform decomposition)
		solutionData.theta.resize(p, 0.0);
		solutionData.thetaNew = solutionData.theta;
		for (int d = 0; d < 3; ++d)
		{
			// Coarse sites in a block if decomposed uniformly
			int uniSpace = static_cast<int>(std::round(domainSize[d] / numCores[d]));

			// Upper edge of block is a variable
			for (i = 0; i < numCores[d] - 1; ++i)
			{
				solutionData.theta[c] = (i + 1) * uniSpace * dh;
				c++;
			}
		}

		// Perturbation vector
		solutionData.delta.resize(p, 0);

		// Populate initial solution vectors and imbalance from uniform decomposition
		mpi_SDCheckDelta(solutionData, dh, numCores);
		mpi_SDComputeImbalance(load, solutionData, numCores);

		// Update uniform decomposition quantity
		load.uniImbalance = load.loadImbalance;
#ifndef L_MPI_TOPOLOGY_REPORT
		L_INFO("Uniform decomposition produces an imbalance of " + std::to_string(load.uniImbalance) + "%.", GridUtils::logfile);
#endif

		// Temporaries
		SDData tempData(solutionData);		// Make a copy
		LoadImbalanceData tmpLoad(load);	// Make a copy

		// Start iteration
		int k = 0;
		double midHeavyBlockX, midHeavyBlockY, midHeavyBlockZ;
		double midCurrentBlockX, midCurrentBlockY, midCurrentBlockZ;
		double dirX, dirY, dirZ;
		while (k < L_MPI_SD_MAX_ITER)
		{

			// Set perturbation directions by driving towards heaviest block
			midHeavyBlockX =
				(tempData.XSol[tmpLoad.heaviestBlock[eXDirection] + 1] + tempData.XSol[tmpLoad.heaviestBlock[eXDirection]]) / 2.0;
			midHeavyBlockY =
				(tempData.YSol[tmpLoad.heaviestBlock[eYDirection] + 1] + tempData.YSol[tmpLoad.heaviestBlock[eYDirection]]) / 2.0;
			midHeavyBlockZ =
				(tempData.ZSol[tmpLoad.heaviestBlock[eZDirection] + 1] + tempData.ZSol[tmpLoad.heaviestBlock[eZDirection]]) / 2.0;

			for (int i = 0; i < numCores[eXDirection]; i++)
			{
				for (int j = 0; j < numCores[eYDirection]; j++)
				{
					for (int k = 0; k < numCores[eZDirection]; k++)
					{
						if (
							i == numCores[eXDirection] - 1 || j == numCores[eYDirection] - 1
#if (L_DIMS == 3)
							|| k == numCores[eZDirection] - 1
#endif
							) continue;

						// Compute middle of current block
						midCurrentBlockX = (tempData.XSol[i + 1] + tempData.XSol[i]) / 2.0;
						midCurrentBlockY = (tempData.YSol[j + 1] + tempData.YSol[j]) / 2.0;
						midCurrentBlockZ = (tempData.ZSol[k + 1] + tempData.ZSol[k]) / 2.0;

						// Compute direction to heaviest block
						dirX = midHeavyBlockX - midCurrentBlockX;
						dirY = midHeavyBlockY - midCurrentBlockY;
						dirZ = midHeavyBlockZ - midCurrentBlockZ;

						// Set deltas					
						if (dirX == 0)
							tempData.delta[i] = -dh;
						else
							tempData.delta[i] = (dirX / std::fabs(dirX)) * dh;

						if (dirY == 0)
							tempData.delta[numCores[eXDirection] - 1 + j] = -dh;
						else
							tempData.delta[numCores[eXDirection] - 1 + j] = (dirY / std::fabs(dirY)) * dh;

#if (L_DIMS == 3)
						if (dirZ == 0)
							tempData.delta[numCores[eXDirection] + numCores[eYDirection] - 2 + k] = -dh;
						else
							tempData.delta[numCores[eXDirection] + numCores[eYDirection] - 2 + k] = (dirZ / std::fabs(dirZ)) * dh;
#endif
					}
				}
			}

			// Check and adjust delta if necessary
			mpi_SDCheckDelta(tempData, dh, numCores);

			// Obtain new imbalance under the adjusted delta
			mpi_SDComputeImbalance(tmpLoad, tempData, numCores);

			// If better than current solution, update
			if (tmpLoad.loadImbalance <= load.loadImbalance)
			{
				load.loadImbalance = tmpLoad.loadImbalance;
				solutionData.XSol = tempData.XSol;
				solutionData.YSol = tempData.YSol;
				solutionData.ZSol = tempData.ZSol;
			}

			// Update theta
			tempData.theta = tempData.thetaNew;

			// Increment k
			k++;
		}
	}

	// Communicate information around topology if not performing a report
	if (!reqDims.size()) mpi_SDCommunicateSolution(solutionData, load.loadImbalance, dh);
	return load;

}
// ************************************************************************** //
/// \brief	Communicate the decomposition result between ranks.
///
///			Uses a broadcast throughout the topolgy to communicate the 
///			information and hence is called by all ranks.
///
///	\param	solutionData	structure to hold SD information.
///	\param	imbalance		percentage of imbalance assocaited with solution.
///	\param	dh				coarse cell spacing.
void MpiManager::mpi_SDCommunicateSolution(SDData& solutionData, double imbalance, double dh)
{

	// Buffer creation for sending/receiving sizes
	int *bufRankSize = new int[num_ranks * L_DIMS];

	// Stash constant pointer to start of buffer
	int * const bufRankSizeStart = bufRankSize;

	if (my_rank == 0)
	{
		// Pack information to send
		for (int i = 0; i < dimensions[eXDirection]; i++)
		{
			for (int j = 0; j < dimensions[eYDirection]; j++)
			{
				for (int k = 0; k < dimensions[eZDirection]; k++)
				{
					// Snap to nearest complete cell
					*bufRankSize++ = static_cast<int>(std::round((solutionData.XSol[i + 1] - solutionData.XSol[i]) / dh));
					*bufRankSize++ = static_cast<int>(std::round((solutionData.YSol[j + 1] - solutionData.YSol[j]) / dh));
#if (L_DIMS == 3)
					*bufRankSize++ = static_cast<int>(std::round((solutionData.ZSol[k + 1] - solutionData.ZSol[k]) / dh));
#endif
				}
			}
		}

		L_INFO("Smart decomposition finished. Imbalance of " + std::to_string(imbalance) + "%.", GridUtils::logfile);

	}

	// Update ranks sizes from solution
	MPI_Bcast(bufRankSizeStart, L_MPI_XCORES * L_MPI_YCORES * L_MPI_ZCORES * L_DIMS, MPI_INT, 0, world_comm);
	int blockIdx = 0;
	bufRankSize = bufRankSizeStart;
	for (int i = 0; i < dimensions[eXDirection]; i++)
	{
		for (int j = 0; j < dimensions[eYDirection]; j++)
		{
			for (int k = 0; k < dimensions[eZDirection]; k++)
			{
				// Update vectors
				cRankSizeX[blockIdx] = *bufRankSize++;
				cRankSizeY[blockIdx] = *bufRankSize++;
#if (L_DIMS == 3)
				cRankSizeZ[blockIdx] = *bufRankSize++;
#endif
				blockIdx++;
			}
		}
	}

	delete[] bufRankSizeStart;

}
// ************************************************************************** //
/// \brief	Writes a report on imbalances from different decomposition topologies.
///
///			This method terminates the application on completion. Only compatible
///			with smart decomposition at present. Uses the values of L_MPI_?CORES
///			as the upper threshold for options.
///
///	\param	dh			coarse cell spacing.
void MpiManager::mpi_reportOnDecomposition(double dh)
{

#if (L_DIMS != 3)
	#undef L_MPI_TOP_ZCORES
	#define L_MPI_TOP_ZCORES 1
#endif

	L_WARN("Topology report mode enabled. No simulation will take place.", GridUtils::logfile);

	if (my_rank == 0)
	{
		// Declarations
		std::vector<int> coreCombo(3);
		std::string logInfo;
		std::ofstream reportFile;
		reportFile.open(GridUtils::path_str + "/topologyreport.out", std::ios::out);
		if (!reportFile.is_open()) L_ERROR("Could not open topology report file. Exiting.", GridUtils::logfile);
		L_INFO("Writing report...", GridUtils::logfile);

		// Write header
		reportFile << "Case\tXCORES\tYCORES\tZCORES\tTotalCore\tImbalance\tUniform\tHeaviestOps\t" << std::endl;

		// Loop over each case
		for (int i = 1; i < L_MPI_TOP_XCORES + 1; ++i)
		{
			for (int j = 1; j < L_MPI_TOP_YCORES + 1; ++j)
			{
				for (int k = 1; k < L_MPI_TOP_ZCORES + 1; ++k)
				{
					coreCombo[eXDirection] = i;
					coreCombo[eYDirection] = j;
					coreCombo[eZDirection] = k;
					LoadImbalanceData load;
					load = mpi_smartDecompose(dh, coreCombo);

					// Log information
					reportFile << std::to_string(
						(k - 1) + (j - 1) * L_MPI_TOP_ZCORES + (i - 1) * L_MPI_TOP_ZCORES * L_MPI_TOP_YCORES
						) + "\t";
					reportFile << std::to_string(i) + "\t";
					reportFile << std::to_string(j) + "\t";
					reportFile << std::to_string(k) + "\t";
					reportFile << std::to_string(i * j * k) + "\t";
					reportFile << std::to_string(load.loadImbalance) + "\t";
					reportFile << std::to_string(load.uniImbalance) + "\t";
					reportFile << std::to_string(load.heaviestOps);
					reportFile << std::endl;
				}
			}
		}

		reportFile.close();
	}

	// Synchronise before trying to exit to ensure report is complete
	MPI_Barrier(world_comm);

	// Exit
	MPI_Finalize();
	exit(EXIT_SUCCESS);
}

// ************************************************************************** //
/// \brief Sets the number of accessible sub-grids on this rank.
///
///			Called by the initialiser once the grid hierarchy has been populated.
void MpiManager::mpi_setSubGridDepth()
{

	// Highest level this rank owns
	int highestLevel = 0;

	// Grid pointer
	GridObj* g = NULL;

	// Loop through all subgrids on this rank
	for (int lev = 1; lev < (L_NUM_LEVELS + 1); lev++)
	{
		// Check if this rank has this level
		g = NULL;
		GridUtils::getGrid(GridManager::getInstance()->Grids, lev, 0, g);

		// Check if g is NULL
		if (g == NULL) break;

		// If it makes it here then it has this level
		highestLevel = lev;
	}

	// Resize vector
	rankGrids.resize(num_ranks, 0);

	// Now pass this value to all ranks
	MPI_Allgather(&highestLevel, 1, MPI_INT, &rankGrids.front(), 1, MPI_INT, world_comm);

	// Split the communicators
	lev_comm.resize(L_NUM_LEVELS + 1);
	int key = my_rank;

	// Now loop through all levels and split communicators
	for (int lev = 0; lev < (L_NUM_LEVELS + 1); lev++)
	{
		// Set colour to undefined
		int colour = MPI_UNDEFINED;

		// If this rank has access to this grid level then it is part of the communicator
		if (rankGrids[my_rank] >= lev)
			colour = 1;

		// Split the communicator
		MPI_Comm_split(world_comm, colour, key, &lev_comm[lev]);
	}
}

// *****************************************************************************
///	\brief	Maps rank numbers from level communicator to world communcator
///
///	\param	level	current grid level
///	\returns		vector of ranks mappings
std::vector<int> MpiManager::mpi_mapRankLevelToWorld(int level) {

	// Declare vector
	std::vector<int> mapping;

	// Loop through all ranks
	for (int rank = 0; rank < num_ranks; rank++) {

		// Check if this rank has this level
		if (rankGrids[rank] >= level)
			mapping.push_back(rank);
	}

	// Return
	return mapping;
}

// *****************************************************************************
///	\brief	Map ranks from global communicator to level rank ID
///
///	\param	level	current grid level
///	\returns		vector of ranks mappings
std::vector<int> MpiManager::mpi_mapRankWorldToLevel(int level) {

	// Declare vector
	std::vector<int> mapping(num_ranks, MPI_UNDEFINED);
	int count = 0;

	// Loop through all ranks
	for (int rank = 0; rank < num_ranks; rank++) {

		// Check if this rank has this level
		if (rankGrids[rank] >= level) {
			mapping[rank] = count;
			count++;
		}
	}

	// Return
	return mapping;
}
