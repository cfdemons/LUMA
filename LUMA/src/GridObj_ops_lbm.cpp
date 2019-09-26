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

/* This file holds all the code for the core LBM operations including collision,
streaming and macroscopic calulcation.
*/

#include "../inc/stdafx.h"
#include "../inc/GridObj.h"
#include "../inc/IVector.h"
#include "../inc/ObjectManager.h"

using namespace std;

// *****************************************************************************
/// \brief	KBC collision operator.
///
///			Applies KBC collision operator using the KBC-N4 and KBC-D models in 
///			3D and 2D, respectively.
///
/// \param i		i-index of lattice site. 
/// \param j		j-index of lattice site.
/// \param k		k-index of lattice site.
/// \param f_new	reference to the temporary, post-collision grid.
void GridObj::LBM_kbcCollide( int i, int j, int k, IVector<double>& f_new ) {
	
	// Declarations
	double ds[L_NUM_VELS], dh[L_NUM_VELS], gamma;

	// Compute required moments and equilibrium moments
#if (L_DIMS == 3)
		
	// Most moments are required in 3D for the KBC-N4 model

	// Stress (second order)
	double M200 = 0.0, M200eq = 0.0;
	double M020 = 0.0, M020eq = 0.0;
	double M002 = 0.0, M002eq = 0.0;
	// Pis (second order)
	double M110 = 0.0, M110eq = 0.0;
	double M101 = 0.0, M101eq = 0.0;
	double M011 = 0.0, M011eq = 0.0;
	// Qs (third order)
	double M111 = 0.0, M111eq = 0.0;
	double M102 = 0.0, M102eq = 0.0;
	double M210 = 0.0, M210eq = 0.0;
	double M021 = 0.0, M021eq = 0.0;
	double M201 = 0.0, M201eq = 0.0;
	double M120 = 0.0, M120eq = 0.0;
	double M012 = 0.0, M012eq = 0.0;

	for (int v = 0; v < L_NUM_VELS; v++) {
		
		// Update feq
		feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) = _LBM_equilibrium_opt(k + j * K_lim + i * K_lim * M_lim, v);

		// These are actually rho * MXXX but no point in dividing to multiply later
		M200 += f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[0][v]);
		M020 += f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[1][v] * c[1][v]);
		M002 += f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[2][v] * c[2][v]);
		M110 += f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[1][v]);
		M101 += f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[2][v]);
		M011 += f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[1][v] * c[2][v]);
		M111 += f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[1][v] * c[2][v]);
		M102 += f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[2][v] * c[2][v]);
		M210 += f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[0][v] * c[1][v]);
		M021 += f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[1][v] * c[1][v] * c[2][v]);
		M201 += f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[0][v] * c[2][v]);
		M120 += f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[1][v] * c[1][v]);
		M012 += f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[1][v] * c[2][v] * c[2][v]);

		M200eq += feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[0][v]);
		M020eq += feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[1][v] * c[1][v]);
		M002eq += feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[2][v] * c[2][v]);
		M110eq += feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[1][v]);
		M101eq += feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[2][v]);
		M011eq += feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[1][v] * c[2][v]);
		M111eq += feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[1][v] * c[2][v]);
		M102eq += feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[2][v] * c[2][v]);
		M210eq += feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[0][v] * c[1][v]);
		M021eq += feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[1][v] * c[1][v] * c[2][v]);
		M201eq += feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[0][v] * c[2][v]);
		M120eq += feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[1][v] * c[1][v]);
		M012eq += feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[1][v] * c[2][v] * c[2][v]);
	}

	// Compute ds
	for (int v = 0; v < L_NUM_VELS; v++) {

		// s part dictated by KBC model choice and directions
		if (c[0][v] == 0 && c[1][v] == 0 && c[2][v] == 0) {

			// First family
			ds[v] = ( -(M200 + M020 + M002) ) - 
					( -(M200eq + M020eq + M002eq) );

		} else if (c[0][v] != 0 && c[1][v] == 0 && c[2][v] == 0) {

			// Second family
			ds[v] =	( (2 * (M200 - M002) - (M020 - M002)) / 6 + (M200 + M020 + M002) / 6  - c[0][v] * 0.5 * (M120 + M102) ) - 
					( (2 * (M200eq - M002eq) - (M020eq - M002eq)) / 6 + (M200eq + M020eq + M002eq) / 6 - c[0][v] * 0.5 * (M120eq + M102eq));

		} else if (c[0][v] == 0 && c[1][v] != 0 && c[2][v] == 0) {

			// Third family
			ds[v] =	( (-(M200 - M002) + 2 * (M020 - M002)) / 6 + (M200 + M020 + M002) / 6 - c[1][v] * 0.5 * (M210 + M012) ) - 
					( (-(M200eq - M002eq) + 2 * (M020eq - M002eq)) / 6 + (M200eq + M020eq + M002eq) / 6 - c[1][v] * 0.5 * (M210eq + M012eq));

		} else if (c[0][v] == 0 && c[1][v] == 0 && c[2][v] != 0) {

			// Fourth family
			ds[v] =	( (-(M200 - M002) - (M020 - M002)) / 6 + (M200 + M020 + M002) / 6 - c[2][v] * 0.5 * (M201 + M021) ) - 
					( (-(M200eq - M002eq) + (M020eq - M002eq)) / 6 + (M200eq + M020eq + M002eq) / 6 - c[2][v] * 0.5 * (M201eq + M021eq) );

		} else if (c[0][v] != 0 && c[1][v] != 0 && c[2][v] == 0) {

			// Fifth family
			ds[v] =	( c[0][v] * c[1][v] * 0.25 * M110 + (c[1][v] * 0.25 * M210 + c[0][v] * 0.25 * M120) ) - 
					( c[0][v] * c[1][v] * 0.25 * M110eq + (c[1][v] * 0.25 * M210eq + c[0][v] * 0.25 * M120eq) );

		} else if (c[0][v] != 0 && c[1][v] == 0 && c[2][v] != 0) {

			// Sixth family
			ds[v] =	( c[0][v] * c[2][v] * 0.25 * M101 + (c[2][v] * 0.25 * M201 + c[0][v] * 0.25 * M102) ) - 
					( c[0][v] * c[2][v] * 0.25 * M101eq + (c[2][v] * 0.25 * M201eq + c[0][v] * 0.25 * M102eq) );

		} else if (c[0][v] == 0 && c[1][v] != 0 && c[2][v] != 0) {

			// Seventh family
			ds[v] =	( c[1][v] * c[2][v] * 0.25 * M011 + (c[2][v] * 0.25 * M021 + c[1][v] * 0.25 * M012) ) - 
					( c[1][v] * c[2][v] * 0.25 * M011eq + (c[2][v] * 0.25 * M021eq + c[1][v] * 0.25 * M012eq) );

		} else if (c[0][v] != 0 && c[1][v] != 0 && c[2][v] != 0) {

			// Eighth family
			ds[v] =	( c[0][v] * c[1][v] * c[2][v] * M111 / 8 ) - 
					( c[0][v] * c[1][v] * c[2][v] * M111eq / 8 );

		}


		// Compute dh
		dh[v] = f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) - feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) - ds[v];

	}



#else

	// Only need to compute these 3 in 2D for KBC-D model
	double M20 = 0.0, M20eq = 0.0;
	double M02 = 0.0, M02eq = 0.0;
	double M11 = 0.0, M11eq = 0.0;

	for (int v = 0; v < L_NUM_VELS; v++) {
		
		// Update feq
		feq(i, j, k, v, M_lim, K_lim, L_NUM_VELS) = _LBM_equilibrium_opt(k + j * K_lim + i * M_lim * K_lim, v);
		
		// These are actually rho * MXX but no point in dividing to multiply later
		M20 += f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[0][v]);
		M02 += f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[1][v] * c[1][v]);
		M11 += f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[1][v]);

		M20eq += feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[0][v]);
		M02eq += feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[1][v] * c[1][v]);
		M11eq += feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) * (c[0][v] * c[1][v]);
	}

	// Compute ds
	for (int v = 0; v < L_NUM_VELS; v++) {

		// s part dictated by KBC model choice and directions
		if (c[0][v] == 0 && c[1][v] == 0) {

			// First family
			ds[v] = 0.0;

		} else if (c[0][v] != 0 && c[1][v] == 0) {

			// Second family
			ds[v] =	( 0.25 * (M20 - M02) ) - 
					( 0.25 * (M20eq - M02eq) );

		} else if (c[0][v] == 0 && c[1][v] != 0) {

			// Third family
			ds[v] =	( -0.25 * (M20 - M02) ) -
					( -0.25 * (M20eq - M02eq) );

		} else {

			// Fourth family
			ds[v] =	( 0.25 * c[0][v] * c[1][v] * M11 ) - 
					( 0.25 * c[0][v] * c[1][v] * M11eq );

		}


		// Compute dh
		dh[v] = f(i,j,k,v,M_lim,K_lim,L_NUM_VELS) - feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS) - ds[v];

	}

#endif

	// Once all dh and ds have been computed, compute the products
	double top_prod = 0.0, bot_prod = 0.0;
	for (int v = 0; v < L_NUM_VELS; v++) {

		// Compute scalar products
		top_prod += ds[v] * dh[v] / feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
		bot_prod += dh[v] * dh[v] / feq(i,j,k,v,M_lim,K_lim,L_NUM_VELS);

	}
	
	// Compute gamma
	if (bot_prod == 0.0) gamma = (2/omega);
	else gamma = (2/omega) - ( 2 - (2/omega) ) * (top_prod / bot_prod);

	// Finally perform collision
	for (int v = 0; v < L_NUM_VELS; v++) {

		// Perform collision
		f_new(i, j, k, v, M_lim, K_lim, L_NUM_VELS) =
			f(i, j, k, v, M_lim, K_lim, L_NUM_VELS) -
			(omega / 2) * (2 * ds[v] + gamma * dh[v])

#if (defined L_GRAVITY_ON || defined L_IBM_ON)
			+ force_i(i,j,k,v,M_lim,K_lim,L_NUM_VELS)
#endif
			;
	}


}


// *****************************************************************************
/// \brief	Site-specific macroscopic update.
///
///			Overload of macroscopic quantity calculation to allow it to be 
///			applied to a single site as used by the MPI unpacking routine to 
///			update the values for the next collision step. This routine does not
///			update the time-averaged quantities.
///
/// \param i		i-index of lattice site. 
/// \param j		j-index of lattice site.
/// \param k		k-index of lattice site.
void GridObj::LBM_macro( int i, int j, int k ) {

	// Declarations
	double rho_temp = 0.0;
	double fux_temp = 0.0;
	double fuy_temp = 0.0;
	double fuz_temp = 0.0;


	if (LatTyp(i,j,k,M_lim,K_lim) == eRefined) {

		// Refined site so set both density and velocity to zero
		rho(i,j,k,M_lim,K_lim) = 0.0;
		u(i,j,k,0,M_lim,K_lim,L_DIMS) = 0.0;
		u(i,j,k,1,M_lim,K_lim,L_DIMS) = 0.0;
#if (L_DIMS == 3)
		u(i,j,k,2,M_lim,K_lim,L_DIMS) = 0.0;
#endif

	} else if (LatTyp(i,j,k,M_lim,K_lim) == eSolid) {

		// Solid site so do not update density but set velocity to zero
		rho(i,j,k,M_lim,K_lim) = 1.0;
		u(i,j,k,0,M_lim,K_lim,L_DIMS) = 0.0;
		u(i,j,k,1,M_lim,K_lim,L_DIMS) = 0.0;
#if (L_DIMS == 3)
		u(i,j,k,2,M_lim,K_lim,L_DIMS) = 0.0;
#endif

	} else if (LatTyp(i,j,k,M_lim,K_lim) == eVelocity) {

		// Velocity BC update themselves prior to collision

	} else {

		// Any other of type of site compute both density and velocity from populations
		rho_temp = 0.0; fux_temp = 0.0; fuy_temp = 0.0; fuz_temp = 0.0;

		for (int v = 0; v < L_NUM_VELS; v++) {

			// Sum up to find mass flux
			fux_temp += (double)c[0][v] * f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
			fuy_temp += (double)c[1][v] * f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);
			fuz_temp += (double)c[2][v] * f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);

			// Sum up to find density
			rho_temp += f(i,j,k,v,M_lim,K_lim,L_NUM_VELS);

		}

		// Assign density
		rho(i,j,k,M_lim,K_lim) = rho_temp;

#if (defined L_GRAVITY_ON || defined L_IBM_ON)
		// Add forces to momentum (rho * time step * 0.5 * force -- eqn 19 in Favier 2014)
		fux_temp += 0.5 * force_xyz(i,j,k,0,M_lim,K_lim,L_DIMS);
		fuy_temp += 0.5 * force_xyz(i,j,k,1,M_lim,K_lim,L_DIMS);
#if (L_DIMS == 3)
		fuz_temp += 0.5 * force_xyz(i,j,k,2,M_lim,K_lim,L_DIMS);
#endif
#endif

		// Assign velocity
		u(i,j,k,0,M_lim,K_lim,L_DIMS) = fux_temp / rho_temp;
		u(i,j,k,1,M_lim,K_lim,L_DIMS) = fuy_temp / rho_temp;
#if (L_DIMS == 3)
		u(i,j,k,2,M_lim,K_lim,L_DIMS) = fuz_temp / rho_temp;
#endif

	}

}
