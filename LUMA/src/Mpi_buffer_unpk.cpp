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
#include "../inc/MpiManager.h"
#include "../inc/GridObj.h"
#include "../inc/GridUtils.h"

// ****************************************************************************
/// \brief	Method to unpack the communication buffer.
///
///			Communication buffer is unpacked onto the supplied grid. Amount and
///			region of unpacking is dictated by the direction of the communication 
///			taking place.
///
/// \param	dir	communication direction.
/// \param	g	grid doing the communication.
void MpiManager::mpi_buffer_unpack( int dir, GridObj* g ) {
	
	// Copy received information back to grid using the EXACT
	// reverse algorithm of the copying procedure
	int i, j , k, v, idx;
	// Local grid sizes for read/writing arrays
	int N_lim = static_cast<int>(g->XInd.size()), M_lim = static_cast<int>(g->YInd.size())
#if (L_dims == 3)
		, K_lim = static_cast<int>(g->ZInd.size());
#else
		, K_lim = 1;
#endif

#ifdef L_MPI_VERBOSE
	*MpiManager::logout << "Unpacking direction " << dir << std::endl;
#endif

	// Copy outgoing information from f_buffer_recv to outer layers
	idx = 0;

	// Switch based on direction
	switch (dir)
	{

	case 0:
		// Right
		
		// Examine possible inner and outer buffer locations
		for (range_i_left) {
			for (j = 0; j < M_lim; j++) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if (  GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum) && 
							(!GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum))
#if (L_dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum))
#endif
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 1:
		// Left
		

		for (range_i_right) {
			for (j = 0; j < M_lim; j++) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if (  GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum) && 
							(!GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum))
#if (L_dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum))
#endif
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}
					
				}
			}
		}
		break;

	case 2:
		// Right-Up
		

		for (range_i_left) {
			for (range_j_down) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if (  GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum) && 
							GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum)
#if (L_dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum))
#endif
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 3:
		// Left-Down
		

		for (range_i_right) {
			for (range_j_up) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if (  GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum) && 
							GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum)
#if (L_dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum))
#endif
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 4:
		// Up
		

		for (i = 0; i < N_lim; i++) {
			for (range_j_down) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if (  (!GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum)) &&
							GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum)
#if (L_dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum))
#endif
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;


	case 5:
		// Down
		

		for (i = 0; i < N_lim; i++) {
			for (range_j_up) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if (  (!GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum)) &&
							GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum)
#if (L_dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum))
#endif
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 6:
		// Left-Up
		

		for (range_i_right) {
			for (range_j_down) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if (  GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum) && 
							GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum)
#if (L_dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum))
#endif
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 7:
		// Right-Down
		

		for (range_i_left) {
			for (range_j_up) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if (  GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum) && 
							GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum)
#if (L_dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum) && !GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum))
#endif
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
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
				for (range_k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
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
				for (range_k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 10:
		// Right-Back
		

		for (range_i_left) {
			for (j = 0; j < M_lim; j++) {
				for (range_k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 11:
		// Left-Front
		

		for (range_i_right) {
			for (j = 0; j < M_lim; j++) {
				for (range_k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 12:
		// Right-Up-Back
		

		for (range_i_left) {
			for (range_j_down) {
				for (range_k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 13:
		// Left-Down-Front
		

		for (range_i_right) {
			for (range_j_up) {
				for (range_k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 14:
		// Up-Back
		

		for (i = 0; i < N_lim; i++) {
			for (range_j_down) {
				for (range_k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 15:
		// Down-Front
		

		for (i = 0; i < N_lim; i++) {
			for (range_j_up) {
				for (range_k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 16:
		// Left-Up-Back
		

		for (range_i_right) {
			for (range_j_down) {
				for (range_k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 17:
		// Right-Down-Front
		

		for (range_i_left) {
			for (range_j_up) {
				for (range_k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 18:
		// Left-Back
		

		for (range_i_right) {
			for (j = 0; j < M_lim; j++) {
				for (range_k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 19:
		// Right-Front
		

		for (range_i_left) {
			for (j = 0; j < M_lim; j++) {
				for (range_k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 20:
		// Left-Down-Back
		

		for (range_i_right) {
			for (range_j_up) {
				for (range_k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 21:
		// Right-Up-Front
		

		for (range_i_left) {
			for (range_j_down) {
				for (range_k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 22:
		// Down-Back
		

		for (i = 0; i < N_lim; i++) {
			for (range_j_up) {
				for (range_k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 23:
		// Up-Front
		

		for (i = 0; i < N_lim; i++) {
			for (range_j_down) {
				for (range_k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum) && !GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 24:
		// Right-Down-Back
		

		for (range_i_left) {
			for (range_j_up) {
				for (range_k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMinimum)) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMaximum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMinimum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
						}
					}

				}
			}
		}
		break;

	case 25:
		// Left-Up-Front
		

		for (range_i_right) {
			for (range_j_down) {
				for (range_k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != eRefined)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],eXDirection,eMaximum)) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],eYDirection,eMinimum)) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],eZDirection,eMaximum))
						) {
							// Must be suitable receiver site
							for (v = 0; v < L_nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,L_nVels) = f_buffer_recv[dir][idx];
								idx++;
							}
							// Update macroscopic (but not time-averaged quantities)
							g->LBM_macro(i,j,k);
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
	*MpiManager::logout << "Unpacking direction " << dir << " complete." << std::endl;
#endif

}
