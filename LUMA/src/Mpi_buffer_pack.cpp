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
#include <mpi.h>
#include "../inc/definitions.h"
#include <iostream>
#include <fstream>
#include "../inc/MpiManager.h"
#include "../inc/GridObj.h"


// Routine to pack the buffer on the supplied grid to be passed in the indicated direction.
void MpiManager::mpi_buffer_pack( int dir, GridObj* g ) {
	
	/* Imagine every grid overlap has an inner region with complete information post-stream
	 * and an outer region with incomplete information post-stream.
	 * In the case of the lower level grids the layers will increase in thickness by a 
	 * factor of 2 with each refinement.
	 * At every exchange, the inner layers need copying from one grid to the outer layer 
	 * of its neighbour on the opposite side of the grid.
	 * To start the process we copy the inner values to the f_buffer_send (intermediate buffer).
	 * This buffer will also be used to receive the new values using MPI_Sendrecv_replace(). */
	int idx, i, j , k, v;
	int N_lim = static_cast<int>(g->XInd.size()), M_lim = static_cast<int>(g->YInd.size())		// Local grid sizes for read/writing arrays
#if (L_dims == 3)
		, K_lim = static_cast<int>(g->ZInd.size());
#else
		, K_lim = 1;
#endif

#ifdef L_MPI_VERBOSE
	*MpiManager::logout << "Packing direction " << dir << std::endl;
#endif

	/* Copy outgoing information from inner layers to f_buffer_send using 
	 * same idenifying logic as the buffer_size routine. */
	idx = 0;


	// Switch based on direction
	switch (dir)
	{

	case 0:
		// Right
		
		// Examine possible inner and outer buffer locations
		for (i_right) {
			for (j = 0; j < M_lim; j++) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if (  GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMaximum) && 
							(!GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum))
#if (L_dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum))
#endif
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 1:
		// Left
		
		for (i_left) {
			for (j = 0; j < M_lim; j++) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if (  GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMinimum) && 
							(!GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum))
#if (L_dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum))
#endif
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 2:
		// Right-Up		

		for (i_right) {
			for (j_up) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if (  GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMaximum) && 
							GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMaximum)
#if (L_dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum))
#endif
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 3:
		// Left-Down		

		for (i_left) {
			for (j_down) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if (  GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMinimum) && 
							GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMinimum)
#if (L_dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum))
#endif
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 4:
		// Up		

		for (i = 0; i < N_lim; i++) {
			for (j_up) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if (  (!GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum)) &&
							GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMaximum)
#if (L_dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum))
#endif
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;


	case 5:
		// Down		

		for (i = 0; i < N_lim; i++) {
			for (j_down) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if (  (!GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum)) &&
							GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMinimum)
#if (L_dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum))
#endif
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 6:
		// Left-Up		

		for (i_left) {
			for (j_up) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if (  GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMinimum) && 
							GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMaximum)
#if (L_dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum))
#endif
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 7:
		// Right-Down		

		for (i_right) {
			for (j_down) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if (  GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMaximum) && 
							GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMinimum)
#if (L_dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum))
#endif
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;


	///////////////////////
	// 3D-specific cases //
	///////////////////////

	case 8:
		// Back		

		for (i = 0; i < N_lim; i++) {
			for (j = 0; j < M_lim; j++) {
				for (k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 9:
		// Front		

		for (i = 0; i < N_lim; i++) {
			for (j = 0; j < M_lim; j++) {
				for (k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 10:
		// Right-Back		

		for (i_right) {
			for (j = 0; j < M_lim; j++) {
				for (k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 11:
		// Left-Front		

		for (i_left) {
			for (j = 0; j < M_lim; j++) {
				for (k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 12:
		// Right-Up-Back		

		for (i_right) {
			for (j_up) {
				for (k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 13:
		// Left-Down-Front		

		for (i_left) {
			for (j_down) {
				for (k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 14:
		// Up-Back		

		for (i = 0; i < N_lim; i++) {
			for (j_up) {
				for (k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 15:
		// Down-Front		

		for (i = 0; i < N_lim; i++) {
			for (j_down) {
				for (k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 16:
		// Left-Up-Back		

		for (i_left) {
			for (j_up) {
				for (k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 17:
		// Right-Down-Front		

		for (i_right) {
			for (j_down) {
				for (k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 18:
		// Left-Back		

		for (i_left) {
			for (j = 0; j < M_lim; j++) {
				for (k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 19:
		// Right-Front		

		for (i_right) {
			for (j = 0; j < M_lim; j++) {
				for (k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 20:
		// Left-Down-Back		

		for (i_left) {
			for (j_down) {
				for (k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 21:
		// Right-Up-Front		

		for (i_right) {
			for (j_up) {
				for (k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 22:
		// Down-Back		

		for (i = 0; i < N_lim; i++) {
			for (j_down) {
				for (k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 23:
		// Up-Front		

		for (i = 0; i < N_lim; i++) {
			for (j_up) {
				for (k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 24:
		// Right-Down-Back		

		for (i_right) {
			for (j_down) {
				for (k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	case 25:
		// Left-Up-Front		

		for (i_left) {
			for (j_up) {
				for (k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be a site to send
							for (v = 0; v < L_nVels; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_nVels);
								idx++;
							}
						}
					}

				}
			}
		}
		break;

	default:
		std::cout << "Unidentified direction: " << dir << std::endl;
		break;
	}

#ifdef L_MPI_VERBOSE
	*MpiManager::logout << "Packing direction " << dir << " complete." << std::endl;
#endif

}