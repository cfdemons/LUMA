#include "stdafx.h"
#include <mpi.h>
#include "definitions.h"
#include <iostream>
#include <fstream>
#include "MPI_manager.h"
#include "GridObj.h"


void MPI_manager::mpi_buffer_pack( int dir, GridObj& Grids ) {

	// Imagine every grid has an inner layer with complete information post-stream 
	// and an outer layer with incomplete information post-stream.
	// The inner layers need copying from one grid to the outer layer of its neighbour on
	// the opposite side of the grid.
	// To start the process we copy the inner values to the f_buffer (intermediate buffer).
	// This buffer will also be used to receive the new values using MPI_Sendrecv_replace().
	unsigned int count, i, j , k, v;

	MPI_Barrier(my_comm);

#ifdef MPI_VERBOSE
	std::ofstream logout;
	logout.open( "./Output/mpiLog_Rank_" + std::to_string(my_rank) + ".out", std::ios::out | std::ios::app );
	logout << "Buffering direction " << dir << std::endl;
	logout.close();
#endif

	// Copy outgoing information from inner layer to f_buffer
	count = 0;
	switch (dir)
	{

	case 0:
		// Right

#if (dims != 3)
		// 2D version //
		// Resize buffer based on y-dimension
		f_buffer.resize( (local_size[1]-2)*nVels );

			
		// Populate the buffer
		i = local_size[0]-2;
		for (j = 1; j < local_size[1]-1; j++) {
			k = 0;
			for (v = 0; v < nVels; v++) {

				f_buffer[count]
					= Grids.f( i, j, k, v, local_size[1], 1, nVels );

				count++;

			}
		}
		break;
#else

		// 3D version //
		// Resize buffer based on y-dimension and z-dimension
		f_buffer.resize( (local_size[1]-2) * (local_size[2]-2) * nVels );

			
		// Populate the buffer
		i = local_size[0]-2;
		for (j = 1; j < local_size[1]-1; j++) {
			for (k = 1; k < local_size[2]-1; k++) {
				for (v = 0; v < nVels; v++) {

					f_buffer[count]
						= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

					count++;
				}
			}
		}
		break;
#endif

			
	case 1:
		// Left
#if (dims != 3)
		// 2D version //
		// Resize buffer based on y-dimension
		f_buffer.resize((local_size[1]-2)*nVels);

		// Populate the buffer
		i = 1;
		for (j = 1; j < local_size[1]-1; j++) {
			k = 0;
			for (v = 0; v < nVels; v++) {

				f_buffer[count]
					= Grids.f( i, j, k, v, local_size[1], 1, nVels);

				count++;

			}
		}
		break;

#else

		// 3D version //
		// Resize buffer based on y-dimension and z-dimension
		f_buffer.resize( (local_size[1]-2) * (local_size[2]-2) * nVels );

			
		// Populate the buffer
		i = 1;
		for (j = 1; j < local_size[1]-1; j++) {
			for (k = 1; k < local_size[2]-1; k++) {
				for (v = 0; v < nVels; v++) {

					f_buffer[count]
						= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

					count++;
				}
			}
		}
		break;
#endif

	case 2:
		// Right-Up
#if (dims != 3)
		// 2D version //
		// Resize buffer to single element of grid
		f_buffer.resize(nVels);

		// Populate the buffer
		i = local_size[0]-2;
		j = 1;
		k = 0;
		for (v = 0; v < nVels; v++) {

			f_buffer[count] 
				= Grids.f( i, j, k, v, local_size[1], 1, nVels);
				
			count++;
		}
		break;
#else
		// 3D version //
		// Resize buffer to single vector of grid
		f_buffer.resize( (local_size[2]-2) * nVels );

		// Populate the buffer
		i = local_size[0]-2;
		j = 1;
		for (k = 1; k < local_size[2]-1; k++) {
			for (v = 0; v < nVels; v++) {

				f_buffer[count] 
					= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels);
				
				count++;
			}
		}
		break;
#endif

	case 3:
		// Left-Down
#if (dims != 3)
		// 2D version //
		// Resize buffer to single element of grid
		f_buffer.resize(nVels);

		// Populate the buffer
		i = 1;
		j = local_size[1]-2;
		k = 0;
		for (v = 0; v < nVels; v++) {

			f_buffer[count] 
				= Grids.f( i, j, k, v, local_size[1], 1, nVels);
				
			count++;
		}
		break;
#else
		// 3D version //
		// Resize buffer to single vector of grid
		f_buffer.resize( (local_size[2]-2) * nVels );

		// Populate the buffer
		i = 1;
		j = local_size[1]-2;
		for (k = 1; k < local_size[2]-1; k++) {
			for (v = 0; v < nVels; v++) {

				f_buffer[count] 
					= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels);
				
				count++;
			}
		}
		break;
#endif

	case 4:
		// Up
#if (dims != 3)
		// 2D version //
		// Resize buffer based on x-dimension
		f_buffer.resize((local_size[0]-2)*nVels);

		// Populate the buffer
		for (i = 1; i < local_size[0]-1; i++) {
			j = 1;
			k = 0;
			for (v = 0; v < nVels; v++) {

				f_buffer[count]
					= Grids.f( i, j, k, v, local_size[1], 1, nVels);

				count++;

			}
		}
		break;
#else

		// 3D version //
		// Resize buffer based on y-dimension and z-dimension
		f_buffer.resize( (local_size[0]-2) * (local_size[2]-2) * nVels );

			
		// Populate the buffer
		for (i = 1; i < local_size[0]-1; i++) {	
			j = 1;
			for (k = 1; k < local_size[2]-1; k++) {
				for (v = 0; v < nVels; v++) {

					f_buffer[count]
						= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

					count++;
				}
			}
		}
		break;
#endif

	case 5:
		// Down
#if (dims != 3)
		// 2D version //
		// Resize buffer based on x-dimension
		f_buffer.resize((local_size[0]-2)*nVels);

		// Populate the buffer
		for (i = 1; i < local_size[0]-1; i++) {
			j = local_size[1]-2;
			k = 0;
			for (v = 0; v < nVels; v++) {

				f_buffer[count]
					= Grids.f( i, j, k, v, local_size[1], 1, nVels);

				count++;

			}
		}
		break;
#else

		// 3D version //
		// Resize buffer based on y-dimension and z-dimension
		f_buffer.resize( (local_size[0]-2) * (local_size[2]-2) * nVels );

			
		// Populate the buffer
		for (i = 1; i < local_size[0]-1; i++) {	
			j = local_size[1]-2;
			for (k = 1; k < local_size[2]-1; k++) {
				for (v = 0; v < nVels; v++) {

					f_buffer[count]
						= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

					count++;
				}
			}
		}
		break;
#endif

	case 6:
		// Left-Up
#if (dims != 3)
		// 2D version //
		// Resize buffer to single element of grid
		f_buffer.resize(nVels);

		// Populate the buffer
		i = 1;
		j = 1;
		k = 0;
		for (v = 0; v < nVels; v++) {

			f_buffer[count] 
				= Grids.f( i, j, k, v, local_size[1], 1, nVels);
				
			count++;
		}

		break;
#else
		// 3D version //
		// Resize buffer to single vector of grid
		f_buffer.resize( (local_size[2]-2) * nVels );

		// Populate the buffer
		i = 1;
		j = 1;
		for (k = 1; k < local_size[2]-1; k++) {
			for (v = 0; v < nVels; v++) {

				f_buffer[count] 
					= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels);
				
				count++;
			}
		}
		break;
#endif

	case 7:
		// Right-Down
#if (dims != 3)
		// 2D version //
		// Resize buffer to single element of grid
		f_buffer.resize(nVels);

		// Populate the buffer
		i = local_size[0]-2;
		j = local_size[1]-2;
		k = 0;
		for (v = 0; v < nVels; v++) {

			f_buffer[count] 
				= Grids.f( i, j, k, v, local_size[1], 1, nVels);
				
			count++;
		}
		break;
#else
		// 3D version //
		// Resize buffer to single vector of grid
		f_buffer.resize( (local_size[2]-2) * nVels );

		// Populate the buffer
		i = local_size[0]-2;
		j = local_size[1]-2;
		for (k = 1; k < local_size[2]-1; k++) {
			for (v = 0; v < nVels; v++) {

				f_buffer[count] 
					= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels);
				
				count++;
			}
		}
		break;
#endif


	///////////////////////
	// 3D-specific cases //
	///////////////////////

	case 8:
		// Back

		// Resize buffer
		f_buffer.resize( (local_size[0]-2) * (local_size[1]-2) * nVels );

		// Populate the buffer
		for (i = 1; i < local_size[0]-1; i++) {	
			for (j = 1; j < local_size[1]-1; j++) {
			k = local_size[2]-2;
				for (v = 0; v < nVels; v++) {

					f_buffer[count]
						= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

					count++;
				}
			}
		}
		break;

	case 9:
		// Front

		// Resize buffer
		f_buffer.resize( (local_size[0]-2) * (local_size[1]-2) * nVels );

		// Populate the buffer
		for (i = 1; i < local_size[0]-1; i++) {	
			for (j = 1; j < local_size[1]-1; j++) {
			k = 1;
				for (v = 0; v < nVels; v++) {

					f_buffer[count]
						= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

					count++;
				}
			}
		}
		break;

	case 10:
		// Right-Back

		// Resize buffer
		f_buffer.resize( (local_size[1]-2) * nVels );

		// Populate the buffer
		i = local_size[0]-2;
		for (j = 1; j < local_size[1]-1; j++) {
			k = local_size[2]-2;
			for (v = 0; v < nVels; v++) {

				f_buffer[count]
					= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

				count++;
			}
		}
		break;

	case 11:
		// Left-Front

		// Resize buffer
		f_buffer.resize( (local_size[1]-2) * nVels );

		// Populate the buffer
		i = 1;
		for (j = 1; j < local_size[1]-1; j++) {
			k = 1;
			for (v = 0; v < nVels; v++) {

				f_buffer[count]
					= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

				count++;
			}
		}
		break;

	case 12:
		// Right-Up-Back

		// Resize buffer
		f_buffer.resize( nVels );

		// Populate the buffer
		i = local_size[0]-2;
		j = 1;
		k = local_size[2]-2;
		for (v = 0; v < nVels; v++) {

			f_buffer[count]
				= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

			count++;
		}
		break;

	case 13:
		// Left-Down-Front

		// Resize buffer
		f_buffer.resize( nVels );

		// Populate the buffer
		i = 1;
		j = local_size[1]-2;
		k = 1;
		for (v = 0; v < nVels; v++) {

			f_buffer[count]
				= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

			count++;
		}
		break;
		
	case 14:
		// Up-Back

		// Resize buffer
		f_buffer.resize( (local_size[0]-2) * nVels );

		// Populate the buffer
		for (i = 1; i < local_size[0]-1; i++) {
			j = 1;
			k = local_size[2]-2;
			for (v = 0; v < nVels; v++) {

				f_buffer[count]
					= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

				count++;
			}
		}
		break;

	case 15:
		// Down-Front

		// Resize buffer
		f_buffer.resize( (local_size[0]-2) * nVels );

		// Populate the buffer
		for (i = 1; i < local_size[0]-1; i++) {
			j = local_size[1]-2;
			k = 1;
			for (v = 0; v < nVels; v++) {

				f_buffer[count]
					= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

				count++;
			}
		}
		break;

	case 16:
		// Left-Up-Back

		// Resize buffer
		f_buffer.resize( nVels );

		// Populate the buffer
		i = 1;
		j = 1;
		k = local_size[2]-2;
		for (v = 0; v < nVels; v++) {

			f_buffer[count]
				= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

			count++;
		}
		break;

	case 17:
		// Right-Down-Front

		// Resize buffer
		f_buffer.resize( nVels );

		// Populate the buffer
		i = local_size[0]-2;
		j = local_size[1]-2;
		k = 1;
		for (v = 0; v < nVels; v++) {

			f_buffer[count]
				= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

			count++;
		}
		break;

	case 18:
		// Left-Back

		// Resize buffer
		f_buffer.resize( (local_size[1]-2) * nVels );

		// Populate the buffer
		i = 1;
		for (j = 1; j < local_size[1]-1; j++) {
			k = local_size[2]-2;
			for (v = 0; v < nVels; v++) {

				f_buffer[count]
					= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

				count++;
			}
		}
		break;

	case 19:
		// Right-Front

		// Resize buffer
		f_buffer.resize( (local_size[1]-2) * nVels );

		// Populate the buffer
		i = local_size[0]-2;
		for (j = 1; j < local_size[1]-1; j++) {
			k = 1;
			for (v = 0; v < nVels; v++) {

				f_buffer[count]
					= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

				count++;
			}
		}
		break;

	case 20:
		// Left-Down-Back

		// Resize buffer
		f_buffer.resize( nVels );

		// Populate the buffer
		i = 1;
		j = local_size[1]-2;
		k = local_size[2]-2;
		for (v = 0; v < nVels; v++) {

			f_buffer[count]
				= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

			count++;
		}
		break;

	case 21:
		// Right-Up-Front

		// Resize buffer
		f_buffer.resize( nVels );

		// Populate the buffer
		i = local_size[0]-2;
		j = 1;
		k = 1;
		for (v = 0; v < nVels; v++) {

			f_buffer[count]
				= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

			count++;
		}
		break;

	case 22:
		// Down-Back

		// Resize buffer
		f_buffer.resize( (local_size[0]-2) * nVels );

		// Populate the buffer
		for (i = 1; i < local_size[0]-1; i++) {
			j = local_size[1]-2;
			k = local_size[2]-2;
			for (v = 0; v < nVels; v++) {

				f_buffer[count]
					= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

				count++;
			}
		}
		break;

	case 23:
		// Up-Front

		// Resize buffer
		f_buffer.resize( (local_size[0]-2) * nVels );

		// Populate the buffer
		for (i = 1; i < local_size[0]-1; i++) {
			j = 1;
			k = 1;
			for (v = 0; v < nVels; v++) {

				f_buffer[count]
					= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

				count++;
			}
		}
		break;

	case 24:
		// Right-Down-Back

		// Resize buffer
		f_buffer.resize( nVels );

		// Populate the buffer
		i = local_size[0]-2;
		j = local_size[1]-2;
		k = local_size[2]-2;
		for (v = 0; v < nVels; v++) {

			f_buffer[count]
				= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

			count++;
		}
		break;

	case 25:
		// Left-Up-Front

		// Resize buffer
		f_buffer.resize( nVels );

		// Populate the buffer
		i = 1;
		j = 1;
		k = 1;
		for (v = 0; v < nVels; v++) {

			f_buffer[count]
				= Grids.f( i, j, k, v, local_size[1], local_size[2], nVels );

			count++;
		}
		break;

	default:
		std::cout << "Unidentified direction: " << dir << std::endl;
		break;
	}

}