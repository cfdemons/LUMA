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


// ****************************************************************************
/// \brief	Method to pack the communication buffer.
///
///			Communication buffer is packed with distribution values from the 
///			supplied grid. Amount of information is dictated by the direction 
///			of the communication being prepared.
///
/// \param	dir	communication direction.
/// \param	g	grid from which information is being sent during the communication.
void MpiManager::mpi_buffer_pack(int dir, GridObj* const g) {
	
	/* Imagine every grid overlap has an inner region with complete information post-stream
	 * and an outer region with incomplete information post-stream.
	 * In the case of the lower level grids the layers will increase in thickness by a 
	 * factor of 2 with each refinement.
	 * At every exchange, the inner layers need copying from one grid to the outer layer 
	 * of its neighbour on the opposite side of the grid.
	 * To start the process we copy the inner values to the f_buffer_send (intermediate buffer).
	 * This buffer will also be used to receive the new values using MPI_Sendrecv_replace(). */
	int idx, i, j , k, v;
	// Local grid sizes for read/writing arrays
	int N_lim = static_cast<int>(g->N_lim), M_lim = static_cast<int>(g->M_lim)
#if (L_DIMS == 3)
		, K_lim = static_cast<int>(g->K_lim);
#else
		, K_lim = 1;
#endif

#ifdef L_MPI_VERBOSE
	*logout << "Packing direction " << dir << std::endl;
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
		for (range_i_right) {
			for (j = 0; j < M_lim; j++) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if (  GridUtils::isOnSenderLayer(g->XPos[i],eXMax) && 
							(!GridUtils::isOnRecvLayer(g->YPos[j],eYMax) && !GridUtils::isOnRecvLayer(g->YPos[j],eYMin))
#if (L_DIMS == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZMax) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZMin))
#endif
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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
		
		for (range_i_left) {
			for (j = 0; j < M_lim; j++) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if (  GridUtils::isOnSenderLayer(g->XPos[i],eXMin) && 
							(!GridUtils::isOnRecvLayer(g->YPos[j],eYMax) && !GridUtils::isOnRecvLayer(g->YPos[j],eYMin))
#if (L_DIMS == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZMax) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZMin))
#endif
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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

		for (range_i_right) {
			for (range_j_up) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if (  GridUtils::isOnSenderLayer(g->XPos[i],eXMax) && 
							GridUtils::isOnSenderLayer(g->YPos[j],eYMax)
#if (L_DIMS == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZMax) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZMin))
#endif
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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

		for (range_i_left) {
			for (range_j_down) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if (  GridUtils::isOnSenderLayer(g->XPos[i],eXMin) && 
							GridUtils::isOnSenderLayer(g->YPos[j],eYMin)
#if (L_DIMS == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZMax) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZMin))
#endif
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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
			for (range_j_up) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if (  (!GridUtils::isOnRecvLayer(g->XPos[i],eXMin) && !GridUtils::isOnRecvLayer(g->XPos[i],eXMax)) &&
							GridUtils::isOnSenderLayer(g->YPos[j],eYMax)
#if (L_DIMS == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZMax) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZMin))
#endif
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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
			for (range_j_down) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if (  (!GridUtils::isOnRecvLayer(g->XPos[i],eXMin) && !GridUtils::isOnRecvLayer(g->XPos[i],eXMax)) &&
							GridUtils::isOnSenderLayer(g->YPos[j],eYMin)
#if (L_DIMS == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZMax) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZMin))
#endif
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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

		for (range_i_left) {
			for (range_j_up) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if (  GridUtils::isOnSenderLayer(g->XPos[i],eXMin) && 
							GridUtils::isOnSenderLayer(g->YPos[j],eYMax)
#if (L_DIMS == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZMax) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZMin))
#endif
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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

		for (range_i_right) {
			for (range_j_down) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if (  GridUtils::isOnSenderLayer(g->XPos[i],eXMax) && 
							GridUtils::isOnSenderLayer(g->YPos[j],eYMin)
#if (L_DIMS == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZMax) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZMin))
#endif
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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
				for (range_k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXMin) && !GridUtils::isOnRecvLayer(g->XPos[i],eXMax)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYMin) && !GridUtils::isOnRecvLayer(g->YPos[j],eYMax)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMax))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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
				for (range_k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXMin) && !GridUtils::isOnRecvLayer(g->XPos[i],eXMax)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYMin) && !GridUtils::isOnRecvLayer(g->YPos[j],eYMax)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMin))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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

		for (range_i_right) {
			for (j = 0; j < M_lim; j++) {
				for (range_k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXMax)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYMin) && !GridUtils::isOnRecvLayer(g->YPos[j],eYMax)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMax))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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

		for (range_i_left) {
			for (j = 0; j < M_lim; j++) {
				for (range_k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXMin)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYMin) && !GridUtils::isOnRecvLayer(g->YPos[j],eYMax)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMin))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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

		for (range_i_right) {
			for (range_j_up) {
				for (range_k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXMax)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYMax)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMax))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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

		for (range_i_left) {
			for (range_j_down) {
				for (range_k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXMin)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYMin)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMin))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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
			for (range_j_up) {
				for (range_k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXMin) && !GridUtils::isOnRecvLayer(g->XPos[i],eXMax)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYMax)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMax))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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
			for (range_j_down) {
				for (range_k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXMin) && !GridUtils::isOnRecvLayer(g->XPos[i],eXMax)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYMin)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMin))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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

		for (range_i_left) {
			for (range_j_up) {
				for (range_k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXMin)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYMax)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMax))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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

		for (range_i_right) {
			for (range_j_down) {
				for (range_k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXMax)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYMin)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMin))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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

		for (range_i_left) {
			for (j = 0; j < M_lim; j++) {
				for (range_k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXMin)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYMin) && !GridUtils::isOnRecvLayer(g->YPos[j],eYMax)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMax))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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

		for (range_i_right) {
			for (j = 0; j < M_lim; j++) {
				for (range_k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXMax)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYMin) && !GridUtils::isOnRecvLayer(g->YPos[j],eYMax)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMin))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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

		for (range_i_left) {
			for (range_j_down) {
				for (range_k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXMin)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYMin)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMax))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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

		for (range_i_right) {
			for (range_j_up) {
				for (range_k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXMax)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYMax)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMin))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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
			for (range_j_down) {
				for (range_k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXMin) && !GridUtils::isOnRecvLayer(g->XPos[i],eXMax)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYMin)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMax))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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
			for (range_j_up) {
				for (range_k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXMin) && !GridUtils::isOnRecvLayer(g->XPos[i],eXMax)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYMax)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMin))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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

		for (range_i_right) {
			for (range_j_down) {
				for (range_k_back) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXMax)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYMin)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMax))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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

		for (range_i_left) {
			for (range_j_up) {
				for (range_k_front) {

					// Check conditions for sender
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Do not pass refined sites as zero anyway
					{
						if ( (GridUtils::isOnSenderLayer(g->XPos[i],eXMin)) && 
								(GridUtils::isOnSenderLayer(g->YPos[j],eYMax)) &&
								(GridUtils::isOnSenderLayer(g->ZPos[k],eZMin))
						) {
							// Must be a site to send
							for (v = 0; v < L_NUM_VELS; v++) {
								f_buffer_send[dir][idx] = g->f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
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
	*logout << "Packing direction " << dir << " complete." << std::endl;
#endif

}