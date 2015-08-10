#include "stdafx.h"
#include "GridObj.h"
#include "definitions.h"
#include "generic_ops.h"
#include <cmath>



// ***************************************************************************************************
// Method to initialise array of IB_bodies
void GridObj::ibm_initialise() {

	// Loop over the number of bodies in the iBody array
	for (size_t ib = 0; ib < iBody.size(); ib++) {

#ifdef IBM_DEBUG
		// DEBUG -- write out marker coordinates
		std::ofstream bodyout;
		bodyout.open("./Output/IBbody_" + std::to_string(ib) + ".out");
		bodyout << "x\ty\tz" << std::endl;
		for (size_t i = 0; i < iBody[ib].markers.size(); i++) {
			bodyout << iBody[ib].markers[i].position[0] << "\t" << iBody[ib].markers[i].position[1] << "\t" << iBody[ib].markers[i].position[2] << std::endl;
		}
		bodyout.close();
#endif

		// Compute support for each marker
		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
			ibm_findsupport(ib, m);	// Pass body ID and marker ID
		}

#ifdef IBM_DEBUG
		// DEBUG -- write out support coordinates
		std::ofstream suppout;
		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
			suppout.open("./Output/Supp_" + std::to_string(ib) + "_" + std::to_string(m) + ".out");
			suppout << "x\ty\tz" << std::endl;
			for (size_t i = 0; i < iBody[ib].markers[m].supp_i.size(); i++) {
				if (dims == 3) {
					suppout << XPos[iBody[ib].markers[m].supp_i[i]] << "\t" << YPos[iBody[ib].markers[m].supp_j[i]] << "\t" << ZPos[iBody[ib].markers[m].supp_k[i]] << std::endl;
				} else {
					suppout << XPos[iBody[ib].markers[m].supp_i[i]] << "\t" << YPos[iBody[ib].markers[m].supp_j[i]] << "\t" << 0.0 << std::endl;
				}
			}
			suppout.close();
		}
#endif

		// Find epsilon for the body
		ibm_findepsilon(ib);

#ifdef IBM_DEBUG
		// DEBUG -- write out epsilon values
		std::ofstream epout;
		epout.open("./Output/Epsilon_" + std::to_string(ib) + ".out");
		for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
			epout << iBody[ib].markers[m].epsilon << std::endl;
		}
		epout.close();
#endif

	}
		
}

// ***************************************************************************************************
// Method to evaluate delta kernel at supplied location
// Both radius and dilation are expected in lattice units
double GridObj::ibm_deltakernel(double radius, double dilation) {

	double mag_r, value;

	// Absolute value of radius
	mag_r = fabs(radius) / dilation;
	
	// Piecemeal function evaluation
	if (mag_r > 1.5) {
		value = 0.0;
	} else if (mag_r > 0.5) {
		value = (5.0 - (3.0 * mag_r) - sqrt(-3.0 * pow( 1.0 - mag_r, 2) + 1.0)) / 6.0;
	} else {
		value = (1.0 + sqrt(1.0 - 3.0 * pow(mag_r, 2) ) ) / 3.0;
	}

	return value;
}

// ***************************************************************************************************
// Method to find the support points in the Eulerian fluid of Lagrange marker m
void GridObj::ibm_findsupport(unsigned int ib, unsigned int m) {

	// Declarations
	unsigned int inear, jnear;					// Nearest node indices
	double dist_x, dist_y, delta_x, delta_y;	// Distances and deltas
#if (dims == 3)
	// Extras for 3D
	double dist_z, delta_z;	
	size_t knear;
#endif


#ifdef CHEAP_NEAREST_NODE_DETECTION
	inear = (unsigned int)floor( iBody[ib].markers[m].position[0]/dx );
	jnear = (unsigned int)floor( iBody[ib].markers[m].position[1]/dy );

#if (dims == 3)
	knear = (unsigned int)floor( iBody[ib].markers[m].position[2]/dz );
#endif


#else

	// Declarations
	double r_min, radius;			// Minimum radius limit and actual radius from centre of kernel

	// Following (Pinelli et al. 2010, JCP) we find closest node
	r_min = abs(b_x-a_x) / dx;	// Initial minimum distance (in lu) taken from problem geometry

	// Loop over the grid to find closest node first
	for (size_t i = 0; i < XPos.size(); i++) {
		for (size_t j = 0; j < YPos.size(); j++) {
			for (size_t k = 0; k < ZPos.size(); k++) {

#if (dims == 3)
				// Find r = sqrt(x^2 + y^2 + z^2)
				radius = vecnorm( XPos[i]-iBody[ib].markers[m].position[0], 
							YPos[j]-iBody[ib].markers[m].position[1],
							ZPos[k]-iBody[ib].markers[m].position[2]
						) / dx;
#else
				// Find r = sqrt(x^2 + y^2)
				radius = vecnorm ( XPos[i]-iBody[ib].markers[m].position[0],
							YPos[j]-iBody[ib].markers[m].position[1]
						) / dx;
				
#endif
				// Check that radius is valid otherwise jacowire must have failed
				if ( _finite(radius) == false ) {
					std::cout << "Jacowire calculation of new position has failed. Exiting." << std::endl;
					system("pause");
					exit(EXIT_FAILURE);					
				}


				if (radius < r_min) {
					r_min = radius;	// This node is closer than the last one so update criterion
					// Store nearest node position
					inear = i;
					jnear = j;
					knear = k;

				}

			}
		}
	}
#endif

	// Define limits of support region (only have to do one since latice is uniformly spaced with dx = dy = dz)
	// Following protocol for arbitrary grid spacing each marker should have at least 3 support nodes in each direction
	double h_plus = std::max( abs(XPos[inear + 1] - XPos[inear]), abs(XPos[inear] - XPos[inear-1]) ) /dx;
	double h_minus = std::min( abs(XPos[inear + 1] - XPos[inear]), abs(XPos[inear] - XPos[inear-1]) ) /dx;

	// Side length of support region defined as 3 x dilation paramter which is found from:
	iBody[ib].markers[m].dilation = (5.0/6.0) * h_plus + (1.0/6.0) * h_minus 
		+ ( (1.0/9.0) * (1 / pow(2,level)) );	// This last term is a small fraction of the local grid spacing in lattice units

	
	// Test to see if required support nodes are available
#if (dims == 3)

	if ( (int)inear - 5 < 0 || inear + 5 >= XPos.size() ) {
		std::cout << "IB body " << std::to_string(ib) << " is too near the X boundary of the grid so support cannot be guaranteed. Exiting." << std::endl;
		system("pause");
		exit(EXIT_FAILURE);
	
	} else if ( (int)jnear - 5 < 0 || jnear + 5 >= YPos.size() ) {
		std::cout << "IB body " << std::to_string(ib) << " is too near the Y boundary of the grid so support cannot be guaranteed. Exiting." << std::endl;
		system("pause");
		exit(EXIT_FAILURE);
		
	} else if ( (int)knear - 5 < 0 || knear + 5 >= ZPos.size() ) {
		std::cout << "IB body " << std::to_string(ib) << " is too near the Z boundary of the grid so support cannot be guaranteed. Exiting." << std::endl;
		system("pause");
		exit(EXIT_FAILURE);
	}


	// Loop over surrounding 5 nodes to find support nodes
	for (size_t i = inear - 5; i <= inear + 5; i++) {
		for (size_t j = jnear - 5; j <= jnear + 5; j++) {
			for (size_t k = knear - 5; k <= knear + 5; k++) {

				// Find distance between Lagrange marker and possible support node and decide whether in cage or not
				if (	( fabs(XPos[inear] - XPos[i])/dx < 1.5*iBody[ib].markers[m].dilation ) &&
						( fabs(YPos[jnear] - YPos[j])/dx < 1.5*iBody[ib].markers[m].dilation ) &&
						( fabs(ZPos[knear] - ZPos[k])/dx < 1.5*iBody[ib].markers[m].dilation )
					) {
						
						// Lies within support region so store information
						iBody[ib].markers[m].supp_i.push_back(i);
						iBody[ib].markers[m].supp_j.push_back(j);
						iBody[ib].markers[m].supp_k.push_back(k);

						// Store normalised area of support region (actually a volume) computed
						// from the local grid spacing in lattice units (dx = 1 / 2^level = 1 on L0)
						iBody[ib].markers[m].local_area = pow( 1 / pow(2,level) ,3) ;

						// Distance between Lagrange marker and support node in lattice units
						dist_x = (XPos[i]-iBody[ib].markers[m].position[0]) / dx;
						dist_y = (YPos[j]-iBody[ib].markers[m].position[1]) / dy;
						dist_z = (ZPos[k]-iBody[ib].markers[m].position[2]) / dz;

						// Store delta function value
						delta_x = ibm_deltakernel(dist_x, iBody[ib].markers[m].dilation);
						delta_y = ibm_deltakernel(dist_y, iBody[ib].markers[m].dilation);
						delta_z = ibm_deltakernel(dist_z, iBody[ib].markers[m].dilation);

						iBody[ib].markers[m].deltaval.push_back( delta_x * delta_y * delta_z );

				}
			}
		}
	}

#else

	// 2D check support region
	if ( (int)inear - 5 < 0 || inear + 5 >= XPos.size() ) {
		std::cout << "IB body " << std::to_string(ib) << " is too near the X boundary of the grid so support cannot be guaranteed. Exiting." << std::endl;
		system("pause");
		exit(EXIT_FAILURE);
	
	} else if ( (int)jnear - 5 < 0 || jnear + 5 >= YPos.size() ) {
		std::cout << "IB body " << std::to_string(ib) << " is too near the Y boundary of the grid so support cannot be guaranteed. Exiting." << std::endl;
		system("pause");
		exit(EXIT_FAILURE);
		
	}
	
	// 2D version to find support nodes
	for (size_t i = inear - 5; i <= inear + 5; i++) {
		for (size_t j = jnear - 5; j <= jnear + 5; j++) {
			size_t k = 0;

			// Find distance between Lagrange marker and possible support node and decide whether in cage or not
			if (	( fabs(XPos[inear] - XPos[i])/dx < 1.5*iBody[ib].markers[m].dilation ) &&
					( fabs(YPos[jnear] - YPos[j])/dx < 1.5*iBody[ib].markers[m].dilation )
				) {
						
					// Lies within support region so store information
					iBody[ib].markers[m].supp_i.push_back(i);
					iBody[ib].markers[m].supp_j.push_back(j);
					iBody[ib].markers[m].supp_k.push_back(k);

					// Store normalised area of support region = dx^2 = (1/2^level) ^ 2.
					iBody[ib].markers[m].local_area = pow( 1 / pow(2,level) ,2) ;

					//  Distance between Lagrange marker and support node in lattice units
					dist_x = (XPos[i]-iBody[ib].markers[m].position[0]) / dx;
					dist_y = (YPos[j]-iBody[ib].markers[m].position[1]) / dy;

					// Store delta function value
					delta_x = ibm_deltakernel(dist_x, iBody[ib].markers[m].dilation);
					delta_y = ibm_deltakernel(dist_y, iBody[ib].markers[m].dilation);

					iBody[ib].markers[m].deltaval.push_back( delta_x * delta_y );

			}
		}
	}
#endif
}

// ***************************************************************************************************
// Method to interpolate the velocity field onto the Lagrange markers
void GridObj::ibm_interpol(unsigned int ib) {

	// Get grid sizes
	size_t M_lim = YPos.size();
#if (dims == 3)
	size_t K_lim = ZPos.size();
#endif
	

	// For each marker
	for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
		
		// Reset the values of interpolated velocity
		std::fill(iBody[ib].markers[m].fluid_vel.begin(), iBody[ib].markers[m].fluid_vel.end(), 0.0);

		// Loop over support nodes
		for (size_t i = 0; i < iBody[ib].markers[m].deltaval.size(); i++) {
			
			// Loop over directions x y z
			for (unsigned int dir = 0; dir < dims; dir++) {

				// Read given velocity component from support node, multiply by delta function
				// for that support node and sum to get interpolated velocity.

#if (dims == 3)
			iBody[ib].markers[m].fluid_vel[dir] += u(	iBody[ib].markers[m].supp_i[i],
														iBody[ib].markers[m].supp_j[i],
														iBody[ib].markers[m].supp_k[i],
														dir,
														M_lim, K_lim, dims
														) * iBody[ib].markers[m].deltaval[i] * iBody[ib].markers[m].local_area;
#else
			iBody[ib].markers[m].fluid_vel[dir] += u(	iBody[ib].markers[m].supp_i[i],
														iBody[ib].markers[m].supp_j[i],
														dir,
														M_lim, dims
														) * iBody[ib].markers[m].deltaval[i] * iBody[ib].markers[m].local_area;
#endif
			}
		}
	}

}

// ***************************************************************************************************
// Method to compute restorative force at each marker in a body
void GridObj::ibm_computeforce(unsigned int ib) {

	// Loop over markers
	for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
		for (unsigned int dir = 0; dir < dims; dir++) {
			// Compute restorative force (in lattice units)
			iBody[ib].markers[m].force_xyz[dir] = (iBody[ib].markers[m].desired_vel[dir] - iBody[ib].markers[m].fluid_vel[dir]) / 
				1 / pow(2,level);	// Time step in lattice units dt = 1 / 2^level = dx
		}
	}

}

// ***************************************************************************************************
// Method to spread the restorative force back on to the fluid sites
void GridObj::ibm_spread(unsigned int ib) {
	
	// For each marker
	for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
		// Loop over support nodes
		for (size_t i = 0; i < iBody[ib].markers[m].deltaval.size(); i++) {

			// Get size of grid
			size_t M_lim = YPos.size();
			size_t K_lim = ZPos.size();
			
			for (size_t dir = 0; dir < dims; dir++) {
				// Add contribution of current marker force to support node Cartesian force vector using delta values computed when support was computed
				force_xyz(iBody[ib].markers[m].supp_i[i], iBody[ib].markers[m].supp_j[i], iBody[ib].markers[m].supp_k[i], dir, M_lim, K_lim, dims) +=
					iBody[ib].markers[m].deltaval[i] * iBody[ib].markers[m].force_xyz[dir] * iBody[ib].markers[m].epsilon * (iBody[ib].spacing/dx);
			
			}
		}
	}

}

// ***************************************************************************************************
// Method to find epsilon
double GridObj::ibm_findepsilon(unsigned int ib) {

	/* The Reproducing Kernel Particle Method (see Pinelli et al. 2010, JCP) requires suitable weighting
	to be computed to ensure conservation while using the interpolation functions. Epsilon is this weighting.
	We can use built-in libraries to solve the ensuing linear system in future. */

	// Declarations
	double Delta_I, Delta_J;

	///////////////////////////////////
	//	Build coefficient matrix A	//
	//		with a_ij values.		//
	///////////////////////////////////

	// Initialise 2D std vector with zeros
	std::vector< std::vector<double> > A (iBody[ib].markers.size(), std::vector<double>(iBody[ib].markers.size(), 0.0) );


	// Loop over support of marker I and integrate delta value multiplied by delta value of marker J.
	for (size_t I = 0; I < iBody[ib].markers.size(); I++) {

		// Loop over markers J
		for (size_t J = 0; J < iBody[ib].markers.size(); J++) {

			// Sum delta values evaluated for each support of I
			for (size_t s = 0; s < iBody[ib].markers[I].supp_i.size(); s++) {

				Delta_I = iBody[ib].markers[I].deltaval[s];
#if (dims == 3)
				Delta_J =	ibm_deltakernel( (iBody[ib].markers[J].position[0] - XPos[iBody[ib].markers[I].supp_i[s]])/dx, iBody[ib].markers[J].dilation ) *
							ibm_deltakernel( (iBody[ib].markers[J].position[1] - YPos[iBody[ib].markers[I].supp_j[s]])/dx, iBody[ib].markers[J].dilation ) *
							ibm_deltakernel( (iBody[ib].markers[J].position[2] - ZPos[iBody[ib].markers[I].supp_k[s]])/dx, iBody[ib].markers[J].dilation );
#else
				Delta_J =	ibm_deltakernel( (iBody[ib].markers[J].position[0] - XPos[iBody[ib].markers[I].supp_i[s]])/dx, iBody[ib].markers[J].dilation ) *
							ibm_deltakernel( (iBody[ib].markers[J].position[1] - YPos[iBody[ib].markers[I].supp_j[s]])/dx, iBody[ib].markers[J].dilation );
#endif
				// Multiply by local area (or volume in 3D)
				A[I][J] += Delta_I * Delta_J * iBody[ib].markers[I].local_area;
			}

			// Multiply by arc length between markers in lattice units
			A[I][J] = A[I][J] * (iBody[ib].spacing/dx);

		}

	}


#ifdef IBM_DEBUG
	// DEBUG -- write out A
	std::ofstream Aout;
	Aout.open("./Output/Amatrix_" + std::to_string(ib) + ".out");
	for (size_t i = 0; i < A.size(); i++) {
		Aout << "\n";
		for (size_t j = 0; j < A.size(); j++) {
			Aout << A[i][j] << "\t";
		}
	}
	Aout.close();
#endif


	// Create vectors
	std::vector<double> epsilon (iBody[ib].markers.size(), 0.0);
	std::vector<double> bVector (iBody[ib].markers.size(), 1.0);
	
	
	///////////////////
	// Solve system //
	//////////////////
		
	// Settings
    double tolerance = 1.0e-4;
	unsigned int maxiterations = 2500;
	double minimum_residual_achieved;
    
    // Biconjugate gradient stabilised method for solving asymmetric linear systems
    minimum_residual_achieved = ibm_bicgstab(A, bVector, epsilon, tolerance, maxiterations);

	// Now assign epsilon to the markers
	for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
		
		iBody[ib].markers[m].epsilon = epsilon[m];
	}

	return minimum_residual_achieved;

}

// ***************************************************************************************************
// Biconjugate gradient stabilised method of solving a linear system -- not my own code but updated from other codes
// Passes references to A, b and epsilon in system A*epsilon = b.
// Tolerance of iterator is tolerance and maximum iterations pursued is maxiterations.
// Returns the minimum residual achieved.
double GridObj::ibm_bicgstab(std::vector< std::vector<double> >& Amatrix, std::vector<double>& bVector, std::vector<double>& epsilon,
						   double tolerance, unsigned int maxiterations) {

	// Declarations //
	
	// Scalars
    double bic_alpha, bic_omega, bic_beta, res_current;
	double res_min = 100.0;		// Arbitrary selection of initial minimum residual -- deliberately big so that first loop is going to generate a epsilon vector with residual better than this.
	// Number of markers
    size_t nls = epsilon.size();
	// Vectors
    std::vector<double> bic_rho (2, 0.0); // Need both i and i-1 instances at same time so need to declare as a 2x1 vector
	std::vector<double> bic_s (nls, 0.0);
	std::vector<double> bic_t (nls, 0.0);
	std::vector<double> epsilon_best (nls, 0.0);
	std::vector<double> bic_r, bic_v, bic_p, bic_rhat;    
	
	// Step 1: Use initial guess to compute r vector
	std::vector<double> bic_Ax = matrix_multiply(Amatrix,epsilon);
	for (size_t i = 0; i < nls; i++) {
		bic_r.push_back(bVector[i] - bic_Ax[i]);
	}

	// Step 2: Choose arbitrary vector r_hat
	bic_rhat = bic_r;

	// Step 3: rho0 = alpha = omega0 = 1
	bic_alpha = 1.0;
	bic_omega = bic_alpha;
    bic_rho[0] = bic_alpha;

	// Step 4: v0 = p0 = 0
	for (size_t i = 0; i < nls; i++) {
        bic_v.push_back(0.0);
		bic_p.push_back(0.0);
    }

	// Step 5: Iterate
	for (unsigned int i = 1; i < maxiterations; i++) {

		// Step 5a: Compute new rho
		bic_rho[1] = dotprod(bic_rhat, bic_r);

		// Step 5b: Compute beta
		bic_beta = (bic_rho[1] / bic_rho[0]) * (bic_alpha / bic_omega);

		// Step 5c: Compute new p vector
		for (size_t j = 0; j < nls; j++) {
			bic_p[j] = bic_r[j] + bic_beta * ( bic_p[j] - bic_omega * bic_v[j] );
		}

		// Step 5d: Compute new v vector
		bic_v = matrix_multiply(Amatrix,bic_p);

		// Step 5e: Compute alpha
		bic_alpha = bic_rho[1] / dotprod(bic_rhat, bic_v);
		// bic_rho is not used again so copy last element back ready for next iteration
		bic_rho[0] = bic_rho[1];
		
		// Step 5f: Compute new s vector
		for (size_t j = 0; j < nls; j++) {
			bic_s[j] = bic_r[j] - (bic_alpha * bic_v[j]);
		}

		// Step 5g: Compute t vector
		bic_t = matrix_multiply(Amatrix,bic_s);

		// Step 5h: Compute new omega
		bic_omega = dotprod(bic_t, bic_s) / dotprod(bic_t, bic_t);

		// Step 5i: Update epsilon
		for (size_t j = 0; j < nls; j++) {
			epsilon[j] += bic_alpha * bic_p[j] + bic_omega * bic_s[j];
		}

		// Step 5j: Compute residual
		for (size_t j = 0; j < nls; j++) {
			bic_r[j] = bic_s[j] - bic_omega * bic_t[j];
		}

		// Step 5k: Check residual and store epsilon if best so far
		res_current = sqrt(dotprod(bic_r, bic_r));
		if ( res_current < res_min) {
			res_min = res_current;	// Note best residual
			for (size_t j = 0; j < nls; j++) {
				epsilon_best[j] = epsilon[j];
			}
			// Check tolerance of best residual
			if ( res_min <= tolerance ) {
				// Before exiting, update epsilon to best epsilon
				for (size_t j = 0; j < nls; j++) {
					epsilon[j] = epsilon_best[j];
				}
				return res_min;
			}
		}

		if (i == maxiterations-1) {
			// Warn that max iterations hit
			std::cout << "Max iterations hit -- values of epsilon might not be converged. Try adjusting the number of Lagrange markers or the grid resolution to adjust support overlap. Setting Epsilon to 2." << std::endl;
			for (size_t j = 0; j < nls; j++) {
				epsilon[j] = 2;
			}
		}
	}


	// Before exiting, update epsilon to best epsilon
	for (size_t j = 0; j < nls; j++) {
		epsilon[j] = epsilon_best[j];
	}
	return res_min;
}

// ***************************************************************************************************