#include "../inc/stdafx.h"
#include <mpi.h>
#include "../inc/definitions.h"
#include <iostream>
#include <fstream>
#include "../inc/MPI_manager.h"
#include "../inc/GridObj.h"

// ********************************************************************************************************

// Routine to unpack the buffer and update the macroscopic quantities at that site based on the new values.
void MPI_manager::mpi_buffer_unpack( int dir, GridObj& Grids ) {

	// Copy received information back to grid using the EXACT
	// reverse algorithm of the copying procedure
	unsigned int count, i, j , k, v;

	MPI_Barrier(my_comm);

#ifdef MPI_VERBOSE
	std::ofstream logout;
	logout.open( timeout_str + "/mpiLog_Rank_" + std::to_string(my_rank) + ".out", std::ios::out | std::ios::app );
	logout << "Unpacking direction " << dir << std::endl;
	logout.close();
#endif

	// Copy outgoing information from f_buffer to outer layer
	count = 0;
	switch (dir)
	{

	case 0:
		// Right

#if (dims != 3)
		// 2D version //
		// Read the buffer (into the left-hand outer layer of the grid)
		i = 0;
		for (j = 1; j < local_size[1]-1; j++) {
			k = 0;
			for (v = 0; v < nVels; v++) {

				Grids.f( i, j, k, v, local_size[1], 1, nVels ) = f_buffer[count];

				count++;

			}
			// Update macroscopic (but not time-averaged quantities)
			Grids.LBM_macro(i,j,k);
		}
		break;
#else
		// 3D version //
		// Read the buffer (into the left-hand outer layer of the grid)
		i = 0;
		for (j = 1; j < local_size[1]-1; j++) {
			for (k = 1; k < local_size[2]-1; k++) {
				for (v = 0; v < nVels; v++) {

					Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

					count++;
				}
				// Update macroscopic (but not time-averaged quantities)
				Grids.LBM_macro(i,j,k);
			}
		}
		break;
#endif

	case 1:
		// Left
#if (dims != 3)
		// 2D version //
		// Read the buffer
		i = local_size[0]-1;
		for (j = 1; j < local_size[1]-1; j++) {
			k = 0;
			for (v = 0; v < nVels; v++) {

				Grids.f( i, j, k, v, local_size[1], 1, nVels) = f_buffer[count];

				count++;

			}
			// Update macroscopic (but not time-averaged quantities)
			Grids.LBM_macro(i,j,k);
		}
		break;
#else
		// 3D version //
		i = local_size[0]-1;
		for (j = 1; j < local_size[1]-1; j++) {
			for (k = 1; k < local_size[2]-1; k++) {
				for (v = 0; v < nVels; v++) {

					Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

					count++;
				}
				// Update macroscopic (but not time-averaged quantities)
				Grids.LBM_macro(i,j,k);
			}
		}
		break;
#endif


	case 2:
		// Right-Up
#if (dims != 3)
		// 2D version //
		// Read the buffer
		i = 0;
		j = local_size[1]-1;
		k = 0;
		for (v = 0; v < nVels; v++) {

			Grids.f( i, j, k, v, local_size[1], 1, nVels) = f_buffer[count];

			count++;
		}
		// Update macroscopic (but not time-averaged quantities)
		Grids.LBM_macro(i,j,k);
		break;
#else
		// 3D version //
		// Read the buffer
		i = 0;
		j = local_size[1]-1;
		for (k = 1; k < local_size[2]-1; k++) {
			for (v = 0; v < nVels; v++) {

				Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

				count++;
			}
			// Update macroscopic (but not time-averaged quantities)
			Grids.LBM_macro(i,j,k);
		}
		break;

#endif


	case 3:
		// Left-Down
#if (dims != 3)
		// 2D version //
		// Read the buffer
		i = local_size[0]-1;
		j = 0;
		k = 0;
		for (v = 0; v < nVels; v++) {

			Grids.f( i, j, k, v, local_size[1], 1, nVels) = f_buffer[count];

			count++;
		}
		// Update macroscopic (but not time-averaged quantities)
		Grids.LBM_macro(i,j,k);
		break;
#else
		// 3D version //
		// Read the buffer
		i = local_size[0]-1;
		j = 0;
		for (k = 1; k < local_size[2]-1; k++) {
			for (v = 0; v < nVels; v++) {

				Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

				count++;
			}
			// Update macroscopic (but not time-averaged quantities)
			Grids.LBM_macro(i,j,k);
		}
		break;
#endif


	case 4:
		// Up
#if (dims != 3)
		// 2D version //
		// Read the buffer
		for (i = 1; i < local_size[0]-1; i++) {
			j = local_size[1]-1;
			k = 0;
			for (v = 0; v < nVels; v++) {

				Grids.f( i, j, k, v, local_size[1], 1, nVels) = f_buffer[count];

				count++;

			}
			// Update macroscopic (but not time-averaged quantities)
			Grids.LBM_macro(i,j,k);
		}
		break;
#else
		// 3D version //
		// Read the buffer
		for (i = 1; i < local_size[0]-1; i++) {
			j = local_size[1]-1;
			for (k = 1; k < local_size[2]-1; k++) {
				for (v = 0; v < nVels; v++) {

					Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

					count++;
				}
				// Update macroscopic (but not time-averaged quantities)
				Grids.LBM_macro(i,j,k);
			}
		}
		break;
#endif


	case 5:
		// Down
#if (dims != 3)
		// 2D version //
		// Read the buffer
		for (i = 1; i < local_size[0]-1; i++) {
			j = 0;
			k = 0;
			for (v = 0; v < nVels; v++) {

				Grids.f( i, j, k, v, local_size[1], 1, nVels) = f_buffer[count];

				count++;

			}
			// Update macroscopic (but not time-averaged quantities)
			Grids.LBM_macro(i,j,k);
		}
		break;
#else
		// 3D version //
		// Read the buffer
		for (i = 1; i < local_size[0]-1; i++) {
			j = 0;
			for (k = 1; k < local_size[2]-1; k++) {
				for (v = 0; v < nVels; v++) {

					Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

					count++;
				}
				// Update macroscopic (but not time-averaged quantities)
				Grids.LBM_macro(i,j,k);
			}
		}
		break;

#endif


	case 6:
		// Left-Up
#if (dims != 3)
		// 2D version //
		// Read the buffer
		i = local_size[0]-1;
		j = local_size[1]-1;
		k = 0;
		for (v = 0; v < nVels; v++) {

			Grids.f( i, j, k, v, local_size[1], 1, nVels) = f_buffer[count];

			count++;
		}
		// Update macroscopic (but not time-averaged quantities)
		Grids.LBM_macro(i,j,k);
		break;
#else
		// 3D version //
		// Read the buffer
		i = local_size[0]-1;
		j = local_size[1]-1;
		for (k = 1; k < local_size[2]-1; k++) {
			for (v = 0; v < nVels; v++) {

				Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

				count++;
			}
			// Update macroscopic (but not time-averaged quantities)
			Grids.LBM_macro(i,j,k);
		}
		break;

#endif


	case 7:
		// Right-Down
#if (dims != 3)
		// 2D version //
		// Read the buffer
		i = 0;
		j = 0;
		k = 0;
		for (v = 0; v < nVels; v++) {

			Grids.f( i, j, k, v, local_size[1], 1, nVels) = f_buffer[count];

			count++;
		}
		// Update macroscopic (but not time-averaged quantities)
		Grids.LBM_macro(i,j,k);
		break;
#else
		// 3D version //
		// Read the buffer
		i = 0;
		j = 0;
		for (k = 1; k < local_size[2]-1; k++) {
			for (v = 0; v < nVels; v++) {

				Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

				count++;
			}
			// Update macroscopic (but not time-averaged quantities)
			Grids.LBM_macro(i,j,k);
		}
		break;

#endif

	///////////////////////
	// 3D-specific cases //
	///////////////////////

	case 8:
		// Back
		for (i = 1; i < local_size[0]-1; i++) {
			for (j = 1; j < local_size[1]-1; j++) {
			k = 0;
				for (v = 0; v < nVels; v++) {

					Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

					count++;
				}
				// Update macroscopic (but not time-averaged quantities)
				Grids.LBM_macro(i,j,k);
			}
		}
		break;

	case 9:
		// Front
		for (i = 1; i < local_size[0]-1; i++) {
			for (j = 1; j < local_size[1]-1; j++) {
			k = local_size[2]-1;
				for (v = 0; v < nVels; v++) {

					Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

					count++;
				}
				// Update macroscopic (but not time-averaged quantities)
				Grids.LBM_macro(i,j,k);
			}
		}
		break;

	case 10:
		// Right-Back
		i = 0;
		for (j = 1; j < local_size[1]-1; j++) {
			k = 0;
			for (v = 0; v < nVels; v++) {

				Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

				count++;
			}
			// Update macroscopic (but not time-averaged quantities)
			Grids.LBM_macro(i,j,k);
		}
		break;

	case 11:
		// Left-Front
		i = local_size[0]-1;
		for (j = 1; j < local_size[1]-1; j++) {
			k = local_size[2]-1;
			for (v = 0; v < nVels; v++) {

				Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

				count++;
			}
			// Update macroscopic (but not time-averaged quantities)
			Grids.LBM_macro(i,j,k);
		}
		break;

	case 12:
		// Right-Up-Back
		i = 0;
		j = local_size[1]-1;
		k = 0;
		for (v = 0; v < nVels; v++) {

			Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

			count++;
		}
		// Update macroscopic (but not time-averaged quantities)
		Grids.LBM_macro(i,j,k);
		break;

	case 13:
		// Left-Down-Front
		i = local_size[0]-1;
		j = 0;
		k = local_size[2]-1;
		for (v = 0; v < nVels; v++) {

			Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

			count++;
		}
		// Update macroscopic (but not time-averaged quantities)
		Grids.LBM_macro(i,j,k);
		break;

	case 14:
		// Up-Back
		for (i = 1; i < local_size[0]-1; i++) {
			j = local_size[1]-1;
			k = 0;
			for (v = 0; v < nVels; v++) {

				Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

				count++;
			}
			// Update macroscopic (but not time-averaged quantities)
			Grids.LBM_macro(i,j,k);
		}
		break;

	case 15:
		// Down-Front
		for (i = 1; i < local_size[0]-1; i++) {
			j = 0;
			k = local_size[2]-1;
			for (v = 0; v < nVels; v++) {

				Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

				count++;
			}
			// Update macroscopic (but not time-averaged quantities)
			Grids.LBM_macro(i,j,k);
		}
		break;

	case 16:
		// Left-Up-Back
		i = local_size[0]-1;
		j = local_size[1]-1;
		k = 0;
		for (v = 0; v < nVels; v++) {

			Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

			count++;
		}
		// Update macroscopic (but not time-averaged quantities)
		Grids.LBM_macro(i,j,k);
		break;

	case 17:
		// Right-Down-Front
		i = 0;
		j = 0;
		k = local_size[2]-1;
		for (v = 0; v < nVels; v++) {

			Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

			count++;
		}
		// Update macroscopic (but not time-averaged quantities)
		Grids.LBM_macro(i,j,k);
		break;

	case 18:
		// Left-Back
		i = local_size[0]-1;
		for (j = 1; j < local_size[1]-1; j++) {
			k = 0;
			for (v = 0; v < nVels; v++) {

				Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

				count++;
			}
			// Update macroscopic (but not time-averaged quantities)
			Grids.LBM_macro(i,j,k);
		}
		break;

	case 19:
		// Right-Front
		i = 0;
		for (j = 1; j < local_size[1]-1; j++) {
			k = local_size[2]-1;
			for (v = 0; v < nVels; v++) {

				Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

				count++;
			}
			// Update macroscopic (but not time-averaged quantities)
			Grids.LBM_macro(i,j,k);
		}
		break;

	case 20:
		// Left-Down-Back
		i = local_size[0]-1;
		j = 0;
		k = 0;
		for (v = 0; v < nVels; v++) {

			Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

			count++;
		}
		// Update macroscopic (but not time-averaged quantities)
		Grids.LBM_macro(i,j,k);
		break;

	case 21:
		// Right-Up-Front
		i = 0;
		j = local_size[1]-1;
		k = local_size[2]-1;
		for (v = 0; v < nVels; v++) {

			Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

			count++;
		}
		// Update macroscopic (but not time-averaged quantities)
		Grids.LBM_macro(i,j,k);
		break;

	case 22:
		// Down-Back
		for (i = 1; i < local_size[0]-1; i++) {
			j = 0;
			k = 0;
			for (v = 0; v < nVels; v++) {

				Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

				count++;
			}
			// Update macroscopic (but not time-averaged quantities)
			Grids.LBM_macro(i,j,k);
		}
		break;

	case 23:
		// Up-Front
		for (i = 1; i < local_size[0]-1; i++) {
			j = local_size[1]-1;
			k = local_size[2]-1;
			for (v = 0; v < nVels; v++) {

				Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

				count++;
			}
			// Update macroscopic (but not time-averaged quantities)
			Grids.LBM_macro(i,j,k);
		}
		break;

	case 24:
		// Right-Down-Back
		i = 0;
		j = 0;
		k = 0;
		for (v = 0; v < nVels; v++) {

			Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

			count++;
		}
		// Update macroscopic (but not time-averaged quantities) (but not time-averaged quantities)
		Grids.LBM_macro(i,j,k);
		break;

	case 25:
		// Left-Up-Front
		i = local_size[0]-1;
		j = local_size[1]-1;
		k = local_size[2]-1;
		for (v = 0; v < nVels; v++) {

			Grids.f( i, j, k, v, local_size[1], local_size[2], nVels ) = f_buffer[count];

			count++;
		}
		// Update macroscopic (but not time-averaged quantities) (but not time-averaged quantities)
		Grids.LBM_macro(i,j,k);
		break;

	default:
		std::cout << "Unidentified direction: " << dir << std::endl;
		break;
	}

}
