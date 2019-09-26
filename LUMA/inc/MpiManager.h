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

#ifndef MPIMAN_H
#define MPIMAN_H

#include "stdafx.h"
#include "HDFstruct.h"
#include "IBInfo.h"
class GridObj;
class GridManager;
class IBBody;


// Define the loop expressions required to inspect the overlap regions of a grid for ease of coding
#define range_i_left	i = 0; i < GridUtils::downToLimit((int)pow(2, g->level + 1), N_lim); i++	///< For loop definition for left halo
#define range_j_down	j = 0; j < GridUtils::downToLimit((int)pow(2, g->level + 1), M_lim); j++	///< For loop definition for bottom halo
#define range_k_front	k = 0; k < GridUtils::downToLimit((int)pow(2, g->level + 1), K_lim); k++	///< For loop definition for front halo
#define range_i_right	i = GridUtils::upToZero(N_lim - (int)pow(2, g->level + 1)); i < N_lim; i++	///< For loop definition for right halo
#define range_j_up		j = GridUtils::upToZero(M_lim - (int)pow(2, g->level + 1)); j < M_lim; j++	///< For loop definition for top halo
#define range_k_back	k = GridUtils::upToZero(K_lim - (int)pow(2, g->level + 1)); k < K_lim; k++	///< For loop definition for back halo

/// \brief	MPI Manager class.
///
///			Class to manage all MPI apsects of the code.
class MpiManager
{
	///	Class to hold smart decomposition data
	class SDData
	{
	public:

		SDData() {};
		~SDData() {};
		// Copy constructor
		SDData(SDData& other)
			: XSol(other.XSol), YSol(other.YSol), ZSol(other.ZSol),
			theta(other.theta), thetaNew(other.thetaNew), delta(other.delta)
		{ };

		// Solution vectors
		std::vector<double> XSol;	///< Solution vector for X variables in SD problem.
		std::vector<double> YSol;	///< Solution vector for Y variables in SD problem.
		std::vector<double> ZSol;	///< Solution vector for Z variables in SD problem.

		// Iteration vectors
		std::vector<double> theta;		///< Complete vector of variable block edges.
		std::vector<double> delta;		///< Perturbation vector.
		std::vector<double> thetaNew;	///< Perturbed vector of block edges.
	};

	/// class to hold imbalance information
	class LoadImbalanceData
	{
	public:
		LoadImbalanceData()
		{
			heaviestBlock.resize(3);
		};
		~LoadImbalanceData() {};

		// Copy constructor
		LoadImbalanceData(LoadImbalanceData& other)
			: loadImbalance(other.loadImbalance), uniImbalance(other.uniImbalance),
			heaviestBlock(other.heaviestBlock), heaviestOps(other.heaviestOps)
		{ };

		double loadImbalance;		///< Imbalance assocaited with smart decomposition.
		double uniImbalance;		///< Imbalance assocaited with uniform decomposition.
		size_t heaviestOps;			///< Number of operations on heaviest rank.
		std::vector<unsigned int> heaviestBlock;	///< MPI indices of the heaviest block.
	};


private :
	MpiManager();			///< Private constructor
	~MpiManager();			///< Private destructor
	static MpiManager* me;	///< Pointer to self

public :

	/************** Member Data **************/

	// MPI world data (all public)
	MPI_Comm world_comm;	///< Global MPI communicator

	/// \brief	Cartesian unit vectors pointing to each neighbour in Cartesian topology.
	///
	///			Define 3D such that first 8 mimic the 2D ones. Opposites are 
	///			simply the next or previous column in the array. MSVC 2013 does
	///			not support initialiser lists tagged onto the constructor although
	///			it is valid C++ so I have had to make it static even though it goes
	///			against the idea of the singleton design.
	static const int neighbour_vectors[3][26];
	int dimensions[3];							///< Size of MPI Cartesian topology
	int neighbour_rank[L_MPI_DIRS];				///< Neighbour rank number for each direction in Cartesian topology
	int neighbour_coords[L_DIMS][L_MPI_DIRS];	///< Coordinates in MPI topology of neighbour ranks

	// Sizes of each of the MPI domains
	/// Number of sites in X direction for each custom rank.
	std::vector<int> cRankSizeX;
	/// Number of sites in Y direction for each custom rank.
	std::vector<int> cRankSizeY;
	/// Number of sites in Z direction for each custom rank.
	std::vector<int> cRankSizeZ;
	
	/// Communicators for sub-grid / region combinations
#if (L_NUM_LEVELS > 0)
	MPI_Comm subGrid_comm[L_NUM_LEVELS * L_NUM_REGIONS];	
#else
	MPI_Comm subGrid_comm[1];	// Default to size = 1
#endif
	
	// Communicators for IBM-level specific communications
	std::vector<MPI_Comm> lev_comm;

	// Commonly used properties of the rank / topology
	int my_rank;				///< Rank number
	int num_ranks;				///< Total number of ranks in MPI Cartesian topology
	int rank_coords[L_DIMS];	///< Coordinates in MPI Cartesian topology


	/// \brief	Absolute positions of edges of the core region represented on this rank.
	///
	///			Excludes outer overlapping layer (recv layer). 
	///			Rows are x,y,z start and end pairs and columns are rank number.
	///			Access the rows using the eCartMinMax enumeration.
	std::vector< std::vector<double> > rank_core_edge;

	/// Vector of size num_ranks which indicates how many sub-grids each rank has access to
	std::vector<int> rankGrids;

	/// \struct HaloEdgeStruct
	/// \brief	Structure containing absolute positions of the edges of halos.
	///
	///			Sender (inner) and receiver (outer) parts of halo are located 
	///			using the convention [left_min left_max right_min right_max] 
	///			for X and similar for Y and Z. Access using the enumerator 
	///			eEdgeMinMax.
	struct HaloEdgeStruct {
		double X[4];	///< X limits
		double Y[4];	///< Y limits
		double Z[4];	///< Z limits
	};
	HaloEdgeStruct sender_layer_pos;	///< Structure containing sender layer edge positions.
	HaloEdgeStruct recv_layer_pos;		///< Structure containing receiver layer edge positions.
	

	// Buffer data
	std::vector< std::vector<double>> f_buffer_send;	///< Array of resizeable outgoing buffers used for data transfer
	std::vector< std::vector<double>> f_buffer_recv;	///< Array of resizeable incoming buffers used for data transfer
	MPI_Status recv_stat;					///< Status structure for Receive return information
	MPI_Request send_requests[L_MPI_DIRS];	///< Array of request structures for handles to posted ISends
	MPI_Status send_stat[L_MPI_DIRS];		///< Array of statuses for each ISend

	/// \struct BufferSizeStruct
	/// \brief	Structure storing buffers sizes in each direction for particular grid.
	struct BufferSizeStruct
	{
		int size[L_MPI_DIRS];	///< Buffer sizes for each direction
		int level;				///< Grid level
		int region;				///< Region number

		BufferSizeStruct(int l, int r) 
			: level(l), region(r){};
	};
	std::vector<BufferSizeStruct> buffer_send_info;	///< Vectors of buffer_info structures holding sender layer size info.
	std::vector<BufferSizeStruct> buffer_recv_info;	///< Vectors of buffer_info structures holding receiver layer size info.

	/// Logfile handle
	std::ofstream* logout;

	// IBM comm classes
	std::vector<std::vector<MarkerCommOwnerSideClass>> markerCommOwnerSide;			///< Owner-side marker-owner comm
	std::vector<std::vector<MarkerCommMarkerSideClass>> markerCommMarkerSide;		///< Marker-side marker-owner comm
	std::vector<std::vector<SupportCommMarkerSideClass>> supportCommMarkerSide;		///< Marker-side marker-support comm
	std::vector<std::vector<SupportCommSupportSideClass>> supportCommSupportSide;	///< Support-side marker-support comm



	/************** Member Methods **************/

	// Singleton design
	static MpiManager* getInstance();	// Get the pointer to the singleton instance (create it if necessary)
	static void destroyInstance();

	// Initialisation
	void mpi_init();												// Initialisation of MpiManager & Cartesian topology
	void mpi_gridbuild(GridManager* const grid_man);				// Do domain decomposition to build local grid dimensions
	void mpi_communicateBlockEdges();								// Get the positional limits of all ranks
	int mpi_buildCommunicators(GridManager* const grid_man);		// Create a new communicator for each sub-grid and region combo
	void mpi_updateLoadInfo(GridManager* const grid_man);			// Method to compute the number of active cells on the rank and pass to master
	void mpi_uniformDecompose(int *numCells);						// Method to perform uniform decomposition into MPI blocks
	LoadImbalanceData mpi_smartDecompose(double dh,
		std::vector<int> combo = std::vector<int>(0));				// Method to perform load-balanced decomposition into MPI blocks
	void mpi_reportOnDecomposition(double dh);						// Method to provide a report on decomposition options
	void mpi_SDReconstructSolution(SDData& solutionData, std::vector<int>& numCores);
	void mpi_SDComputeImbalance(LoadImbalanceData& load, SDData& solutionData, std::vector<int>& numCores);
	bool mpi_SDCheckDelta(SDData& solutionData, double dh, std::vector<int>& numCores);
	void mpi_SDCommunicateSolution(SDData& solutionData, double imbalance, double dh);
	void mpi_setSubGridDepth();										// Method to initialise the rankGrids variable

	// Helper functions
	std::vector<int> mpi_mapRankLevelToWorld(int level);			// Map rank numbers from level communicator to world communcator
	std::vector<int> mpi_mapRankWorldToLevel(int level);			// Map rank numbers from world communicator to level communicator

	// Buffer methods
	void mpi_buffer_pack(int dir, GridObj* const g);		// Pack the buffer ready for data transfer on the supplied grid in specified direction
	void mpi_buffer_unpack(int dir, GridObj* const g);		// Unpack the buffer back to the grid given
	void mpi_buffer_size();									// Set buffer size information for grids in hierarchy given and 
															// set pointer to hierarchy for subsequent access
	void mpi_buffer_size_send( GridObj* const g );			// Routine to find the size of the sending buffer on supplied grid
	void mpi_buffer_size_recv( GridObj* const g );			// Routine to find the size of the receiving buffer on supplied grid

	// IO
	void mpi_writeout_buf(std::string filename, int dir);		// Write out the buffers of direction dir to file

	// Comms
	void mpi_communicate( int level, int regnum );		// Wrapper routine for communication between grids of given level/region
	int mpi_getOpposite(int direction);					// Version of GridUtils::getOpposite for MPI_directions rather than lattice directions

	// IBM
	void mpi_buildMarkerComms(int level);												// Build comms required for epsilon calculation
	void mpi_buildSupportComms(int level);												// Build comms required for support communication
	void mpi_epsilonCommGather(int level);												// Do communication required for epsilon calculation
	void mpi_epsilonCommScatter(int level);												// Do communication required for epsilon calculation
	void mpi_uniEpsilonCommGather(int level, int rootRank, IBBody &iBodyTmp);			// Do communication required for universal epsilon calculation
	void mpi_uniEpsilonCommScatter(int level, int rootRank, IBBody &iBodyTmp);			// Do communication required for universal epsilon calculation
	void mpi_interpolateComm(int level, std::vector<std::vector<double>> &interpVels);	// Do communication required for velocity interpolation
	void mpi_spreadComm(int level, std::vector<std::vector<double>> &spreadForces);		// Do communication required for force spreading
	void mpi_dsCommScatter(int level);													// Spread the ds values from owner to other ranks
	void mpi_ptCloudMarkerGather(IBBody *iBody, std::vector<double> &recvPositionBuffer, std::vector<int> &recvIDBuffer, std::vector<int> &recvSizeBuffer, std::vector<int> &recvDisps);		// Gather in info for pt cloud sorter
	void mpi_ptCloudMarkerScatter(IBBody *iBody, std::vector<int> &recvIDBuffer, std::vector<int> &recvSizeBuffer, std::vector<int> &recvDisps);	// Scatter info for pt cloud sorter

	// FEM
	void mpi_forceCommGather(int level);
	void mpi_spreadNewMarkers(int level, std::vector<std::vector<int>> &markerIDs, std::vector<std::vector<std::vector<double>>> &positions, std::vector<std::vector<std::vector<double>>> &vels);
};

#endif
