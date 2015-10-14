#pragma once

#include "definitions.h"
#include "GridObj.h"

class MPI_manager
{

	// Default constructor and destructor
public:
	MPI_manager(void);
	~MPI_manager(void);

	/*	
	***************************************************************************************************************
	********************************************* Member Data *****************************************************
	***************************************************************************************************************
	*/

	// MPI world data (all public)

	MPI_Comm my_comm;						// MPI communicator
	static const int MPI_cartlab[3][26];	// Cartesian unit vectors pointing to each neighbour in Cartesian topology
	int MPI_dims[dims];						// Size of MPI Cartesian topology
	int my_rank;							// Rank number
	int num_ranks;							// Total number of ranks in MPI Cartesian topology
	int MPI_coords[dims];					// Coordinates in MPI Cartesian topolgy
	int neighbour_rank[MPI_dir];			// Neighbour rank number for each direction in Cartesian topology
	int neighbour_coords[dims][MPI_dir];	// Coordinates in MPI topology of neighbour ranks


	// Grid data
	int global_dims[3];						// Dimensions of problem lattice
	std::vector<unsigned int> local_size;	// Dimensions of lattice represented on this rank (includes inner and outer overlapping layers)
	// Global indices of lattice represented on this rank (excluding outer overlapping layer)
	std::vector< std::vector<unsigned int> > global_edge_ind;	// Rows are x,y,z start and end pairs and columns are rank number
	// Global positions of lattice represented on this rank (excluding outer overlapping layer)
	std::vector< std::vector<double> > global_edge_pos;	// Rows are x,y,z start and end pairs and columns are rank number
	
	// Buffer data
	ivector<double> f_buffer;				// Resizeable buffer used for data transfer
	MPI_Status stat;						// Status structure for Send-Receive return information

	/*	
	***************************************************************************************************************
	********************************************* Member Methods **************************************************
	***************************************************************************************************************
	*/

	// Initialisation
	void mpi_init( );		// Initialisation of MPI_manager & Cartesian topology
	void mpi_gridbuild( );	// Do domain decomposition to build local grid dimensions

	// Buffer methods
	void mpi_buffer_pack( int dir, GridObj& Grids );	// Pack the buffer ready for data transfer in specified direction
	void mpi_buffer_unpack( int dir, GridObj& Grids );	// Unpack the buffer back to the grid

	// IO
	void writeout_buf( std::string filename );		// Write out the buffer to file
};

