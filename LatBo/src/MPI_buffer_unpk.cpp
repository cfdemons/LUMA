#include "../inc/stdafx.h"
#include <mpi.h>
#include "../inc/definitions.h"
#include <iostream>
#include <fstream>
#include "../inc/MpiManager.h"
#include "../inc/GridObj.h"

// ********************************************************************************************************

// Routine to unpack the buffer and update the macroscopic quantities at that site based on the new values.
void MpiManager::mpi_buffer_unpack( int dir, GridObj* g ) {
	
	// Copy received information back to grid using the EXACT
	// reverse algorithm of the copying procedure
	int i, j , k, v, idx;
	int N_lim = g->XInd.size(), M_lim = g->YInd.size()		// Local grid sizes for read/writing arrays
#if (dims == 3)
		, K_lim = g->ZInd.size();
#else
		, K_lim = 1;
#endif

#ifdef MPI_VERBOSE
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
		for (i_left) {
			for (j = 0; j < M_lim; j++) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if (  GridUtils::isOnRecvLayer(g->XPos[i],'x',"min") && 
							(!GridUtils::isOnRecvLayer(g->YPos[j],'y',"min") && !GridUtils::isOnRecvLayer(g->YPos[j],'y',"max"))
#if (dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],'z',"min") && !GridUtils::isOnRecvLayer(g->ZPos[k],'z',"max"))
#endif
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
		

		for (i_right) {
			for (j = 0; j < M_lim; j++) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if (  GridUtils::isOnRecvLayer(g->XPos[i],'x',"max") && 
							(!GridUtils::isOnRecvLayer(g->YPos[j],'y',"min") && !GridUtils::isOnRecvLayer(g->YPos[j],'y',"max"))
#if (dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],'z',"min") && !GridUtils::isOnRecvLayer(g->ZPos[k],'z',"max"))
#endif
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
		

		for (i_left) {
			for (j_down) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if (  GridUtils::isOnRecvLayer(g->XPos[i],'x',"min") && 
							GridUtils::isOnRecvLayer(g->YPos[j],'y',"min")
#if (dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],'z',"min") && !GridUtils::isOnRecvLayer(g->ZPos[k],'z',"max"))
#endif
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
		

		for (i_right) {
			for (j_up) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if (  GridUtils::isOnRecvLayer(g->XPos[i],'x',"max") && 
							GridUtils::isOnRecvLayer(g->YPos[j],'y',"max")
#if (dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],'z',"min") && !GridUtils::isOnRecvLayer(g->ZPos[k],'z',"max"))
#endif
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
			for (j_down) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if (  (!GridUtils::isOnRecvLayer(g->XPos[i],'x',"max") && !GridUtils::isOnRecvLayer(g->XPos[i],'x',"min")) &&
							GridUtils::isOnRecvLayer(g->YPos[j],'y',"min")
#if (dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],'z',"min") && !GridUtils::isOnRecvLayer(g->ZPos[k],'z',"max"))
#endif
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
			for (j_up) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if (  (!GridUtils::isOnRecvLayer(g->XPos[i],'x',"max") && !GridUtils::isOnRecvLayer(g->XPos[i],'x',"min")) &&
							GridUtils::isOnRecvLayer(g->YPos[j],'y',"max")
#if (dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],'z',"min") && !GridUtils::isOnRecvLayer(g->ZPos[k],'z',"max"))
#endif
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
		

		for (i_right) {
			for (j_down) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if (  GridUtils::isOnRecvLayer(g->XPos[i],'x',"max") && 
							GridUtils::isOnRecvLayer(g->YPos[j],'y',"min")
#if (dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],'z',"min") && !GridUtils::isOnRecvLayer(g->ZPos[k],'z',"max"))
#endif
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
		

		for (i_left) {
			for (j_up) {
				for (k = 0; k < K_lim; k++) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if (  GridUtils::isOnRecvLayer(g->XPos[i],'x',"min") && 
							GridUtils::isOnRecvLayer(g->YPos[j],'y',"max")
#if (dims == 3)
							&&
							(!GridUtils::isOnRecvLayer(g->ZPos[k],'z',"min") && !GridUtils::isOnRecvLayer(g->ZPos[k],'z',"max"))
#endif
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
				for (k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],'x',"max") && !GridUtils::isOnRecvLayer(g->XPos[i],'x',"min")) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],'y',"max") && !GridUtils::isOnRecvLayer(g->YPos[j],'y',"min")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"min"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
				for (k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],'x',"max") && !GridUtils::isOnRecvLayer(g->XPos[i],'x',"min")) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],'y',"max") && !GridUtils::isOnRecvLayer(g->YPos[j],'y',"min")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"max"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
		

		for (i_left) {
			for (j = 0; j < M_lim; j++) {
				for (k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],'x',"min")) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],'y',"max") && !GridUtils::isOnRecvLayer(g->YPos[j],'y',"min")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"min"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
		

		for (i_right) {
			for (j = 0; j < M_lim; j++) {
				for (k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],'x',"max")) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],'y',"max") && !GridUtils::isOnRecvLayer(g->YPos[j],'y',"min")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"max"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
		

		for (i_left) {
			for (j_down) {
				for (k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],'x',"min")) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],'y',"min")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"min"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
		

		for (i_right) {
			for (j_up) {
				for (k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],'x',"max")) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],'y',"max")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"max"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
			for (j_down) {
				for (k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],'x',"max") && !GridUtils::isOnRecvLayer(g->XPos[i],'x',"min")) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],'y',"min")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"min"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
			for (j_up) {
				for (k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],'x',"max") && !GridUtils::isOnRecvLayer(g->XPos[i],'x',"min")) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],'y',"max")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"max"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
		

		for (i_right) {
			for (j_down) {
				for (k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],'x',"max")) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],'y',"min")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"min"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
		

		for (i_left) {
			for (j_up) {
				for (k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],'x',"min")) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],'y',"max")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"max"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
		

		for (i_right) {
			for (j = 0; j < M_lim; j++) {
				for (k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],'x',"max")) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],'y',"max") && !GridUtils::isOnRecvLayer(g->YPos[j],'y',"min")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"min"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
		

		for (i_left) {
			for (j = 0; j < M_lim; j++) {
				for (k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],'x',"min")) && 
								(!GridUtils::isOnRecvLayer(g->YPos[j],'y',"max") && !GridUtils::isOnRecvLayer(g->YPos[j],'y',"min")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"max"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
		

		for (i_right) {
			for (j_up) {
				for (k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],'x',"max")) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],'y',"max")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"min"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
		

		for (i_left) {
			for (j_down) {
				for (k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],'x',"min")) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],'y',"min")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"max"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
			for (j_up) {
				for (k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],'x',"max") && !GridUtils::isOnRecvLayer(g->XPos[i],'x',"min")) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],'y',"max")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"min"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
			for (j_down) {
				for (k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (!GridUtils::isOnRecvLayer(g->XPos[i],'x',"max") && !GridUtils::isOnRecvLayer(g->XPos[i],'x',"min")) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],'y',"min")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"max"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
		

		for (i_left) {
			for (j_up) {
				for (k_front) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],'x',"min")) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],'y',"max")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"min"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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
		

		for (i_right) {
			for (j_down) {
				for (k_back) {

					// Check conditions for receiver
					if (g->LatTyp(i,j,k,M_lim,K_lim) != 2)	// Refined sites are not passed
					{
						if ( (GridUtils::isOnRecvLayer(g->XPos[i],'x',"max")) && 
								(GridUtils::isOnRecvLayer(g->YPos[j],'y',"min")) &&
								(GridUtils::isOnRecvLayer(g->ZPos[k],'z',"max"))
						) {
							// Must be suitable receiver site
							for (v = 0; v < nVels; v++) {
								g->f(i,j,k,v,M_lim,K_lim,nVels) = f_buffer_recv[dir][idx];
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

#ifdef MPI_VERBOSE
	*MpiManager::logout << "Unpacking direction " << dir << " complete." << std::endl;
#endif

}
