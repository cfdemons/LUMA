#include "../inc/stdafx.h"
#include "../inc/GridObj.h"
#include "../inc/definitions.h"


/*
This whole file was originally lifted from the 2D IBLB code.
It needs completely rewriting to generalise.
*/



// **************************************************************************************** //
// Jacowire function for flexible cilia with one simply supported end and one free end
void GridObj::ibm_jacowire( unsigned int ib ) {

    ///////// Initialisation /////////    
    double tolerance = 1.0e-4;		// Tolerance of iterative solver
    double residual = 100;			// Set to arbitrary large number to being with
	unsigned int max_terations = iBody[ib].markers.size();	// Set maximum number of iterations
	double Froude;					// Ri g_vec / g = (Fr,0) and Fr = u^2 / gL
	size_t n = iBody[ib].markers.size() - 1;	// Number of markers - 1
	size_t i = 0;								// Iterator
	double EI = iBody[ib].flexural_rigidity;	// Flexural rigidity

	// Special variables for the Numerical Recipe code
	unsigned long *indx;
		indx = new unsigned long[3 * iBody[ib].markers.size()];
	unsigned int m1 = 3, m2 = 3;
	double **AA, **AL, *d, *res;
		AA = new double*[3 * iBody[ib].markers.size()+1];
		AL = new double*[3 * iBody[ib].markers.size()];
		d = new double[1];
	for (i = 0; i < 3 * iBody[ib].markers.size() ; ++i) {
        AA[i] = new double[m1 + m2 + 1];
        AL[i] = new double[m1 + 1];
    }
	res = new double[3 * iBody[ib].markers.size()];

		
	// Reset tension vector
	std::fill(iBody[ib].tension.begin(), iBody[ib].tension.end(), 0.0);

	// Define quantities used in the routines below
	// Length of filament in lu
	double length_lu = gUtils.vecnorm(iBody[ib].markers[iBody[ib].markers.size()-1].position[0] - iBody[ib].markers[0].position[0],
		iBody[ib].markers[iBody[ib].markers.size()-1].position[1] - iBody[ib].markers[0].position[1],
		iBody[ib].markers[iBody[ib].markers.size()-1].position[2] - iBody[ib].markers[0].position[2]) / dx;

	// Marker spacing in normalised instrinsic coordinates
	double ds_nondim = iBody[ib].spacing / (length_lu*dx);

	// Square of non-dimensional spacing
	double ds_sqrd = pow(ds_nondim,2);

#ifdef GRAVITY_ON
	Froude = pow(gUtils.vecnorm(u_0x,u_0y,u_0z),2) / grav_force * length_lu;
#else
	Froude = 0.0;
#endif

	// Beta = spacing^2 (lu) / reference time^2 (lu) = ds_nondim^2 / (dt / (length (lu) / u (lu) )^2 )
	double beta = ds_sqrd / pow( ( (1 / pow(2,level)) / (length_lu / gUtils.vecnorm(u_0x,u_0y,u_0z)) ), 2);


	// Simply supported end position in filament-normalised coordinates units and tension in between it and the next marker
	double Tension0 = 2 * iBody[ib].tension[0] - iBody[ib].tension[1];
	double x0 = 0;
	double y0 = 0;

	// Get initial positions of markers in filament-normalised coordinates [0,+/-1] for below
	std::vector<double> x, y;
	for (i = 0; i < n; i++) {
		x.push_back( (iBody[ib].markers[i+1].position[0] - iBody[ib].markers[0].position[0]) / (length_lu * dx) );
		y.push_back( (iBody[ib].markers[i+1].position[1] - iBody[ib].markers[0].position[1]) / (length_lu * dy) );
	}

	// Create body force vectors
	std::vector<double> Fx(iBody[ib].markers.size(), 0.0), Fy(iBody[ib].markers.size(), 0.0);

	// Populate force vectors with non-dimensional forces (divide by spacing/dx = marker spacing in lattice units)
	double Fref = iBody[ib].delta_rho * pow(gUtils.vecnorm(u_0x,u_0y,u_0z), 2);
	for (i = 0; i < Fx.size(); i++) {
		Fx[i] = -iBody[ib].markers[i].force_xyz[0] / (Fref / (iBody[ib].markers[i].epsilon / (iBody[ib].spacing/dx) ) );
		Fy[i] = -iBody[ib].markers[i].force_xyz[1] / (Fref / (iBody[ib].markers[i].epsilon / (iBody[ib].spacing/dx) ) );
	}
	
	// Normalised position (x,y,z)/L at t+1 for each marker are computed using extrapolation of the current (t) position 
	// and previous (t-1) position_old at neighbour node of each marker
	std::vector<double> xstar, ystar;
	for (i = 0; i < n; i++) {
		xstar.push_back( ((2 * iBody[ib].markers[i+1].position[0] - iBody[ib].markers[i+1].position_old[0]) - iBody[ib].markers[0].position[0]) / (dx * length_lu) );
		ystar.push_back( ((2 * iBody[ib].markers[i+1].position[1] - iBody[ib].markers[i+1].position_old[1]) - iBody[ib].markers[0].position[1]) / (dy * length_lu)  );
	}


	///////// Create vector of unknowns /////////
	std::vector<double> G( 3*iBody[ib].markers.size(), 0.0);

	



	/*
	The following is taken "as-is" from the 2D IB-LB code and has not been modified.
	This needs fully understanding in the future to allow different boundary conditions
	to be incorporated.
	*/
	// ********************************* Start of copy and paste ************************************** //




	// Assemble the part of rhs that stays invariant under Newton iterates
    G[1] = (Fx[0] + Froude) * x0 + Fy[0] * y0; // SIMPLY SUPPORTED constant part of Equation (15) [Resolution of the filament problem]

    int ii = -1;
    for (i = 0; i < n; ++i) {
        ii += 3;
        G[ii] = ds_sqrd;
    }
    // first two Lagrangian points: Equations (10) [Resolution of the filament problem]
    G[3] = EI * (-2 * x0 + 5 * xstar[0] - 4 * xstar[1] + xstar[2]) / ds_sqrd - ds_sqrd * (Fx[1] + Froude) - beta * xstar[0];
    G[4] = EI * (-2 * y0 + 5 * ystar[0] - 4 * ystar[1] + ystar[2]) / ds_sqrd - ds_sqrd * Fy[1] - beta * ystar[0];
    
    G[6] = EI * (x0 - 4 * xstar[0] + 6 * xstar[1] - 4 * xstar[2] + xstar[3]) / ds_sqrd - ds_sqrd * (Fx[2] + Froude) - beta * xstar[1];
    G[7] = EI * (y0 - 4 * ystar[0] + 6 * ystar[1] - 4 * ystar[2] + ystar[3]) / ds_sqrd - ds_sqrd * Fy[2] - beta * ystar[1];
    
    // bulk of Lagrangian points
    ii = 8;
    for (i = 2; i < n - 2; ++i) {
        // Unknown in Equation (3) [Resolution of the filament problem]
        G[ii + 1] = EI * (xstar[i - 2] - 4 * xstar[i - 1] + 6 * xstar[i] - 4 * xstar[i + 1] + xstar[i + 2]) / ds_sqrd - ds_sqrd * (Fx[i + 1] + Froude) - beta * xstar[i];
        // Unknown in Equation (4) [Resolution of the filament problem]
        G[ii + 2] = EI * (ystar[i - 2] - 4 * ystar[i - 1] + 6 * ystar[i] - 4 * ystar[i + 1] + ystar[i + 2]) / ds_sqrd - ds_sqrd * Fy[i + 1] - beta * ystar[i];
        ii = ii + 3;
    }
    
    // last but one Lagrangian point
    i = n - 2;
    double xsp2 = 2 * xstar[n - 1] - xstar[n - 2];
    double ysp2 = 2 * ystar[n - 1] - ystar[n - 2];
    G[ii + 1] = EI * (xstar[i - 2] - 4 * xstar[i - 1] + 6 * xstar[i] - 4 * xstar[i + 1] + xsp2) / ds_sqrd - ds_sqrd * (Fx[i + 1] + Froude) - beta * xstar[i];
    G[ii + 2] = EI * (ystar[i - 2] - 4 * ystar[i - 1] + 6 * ystar[i] - 4 * ystar[i + 1] + ysp2) / ds_sqrd - ds_sqrd * Fy[i + 1] - beta * ystar[i];
    
    // last Lagrangian point
    ii = ii + 3;
    i = n - 1;
    xsp2 = 3 * xstar[n - 1] - 2 * xstar[n - 2];
    ysp2 = 3 * ystar[n - 1] - 2 * ystar[n - 2];
    double xsp1 = 2 * xstar[n - 1] - xstar[n - 2];
    double ysp1 = 2 * ystar[n - 1] - ystar[n - 2];
    G[ii + 1] = EI * (xstar[i - 2] - 4 * xstar[i - 1] + 6 * xstar[i] - 4 * xsp1 + xsp2) / ds_sqrd - ds_sqrd * (Fx[i + 1] + Froude) - beta * xstar[i];
    G[ii + 2] = EI * (ystar[i - 2] - 4 * ystar[i - 1] + 6 * ystar[i] - 4 * ysp1 + ysp2) / ds_sqrd - ds_sqrd * Fy[i + 1] - beta * ystar[i];


#ifdef IBM_DEBUG
		// DEBUG -- write out G vector
		std::ofstream Gout;
		Gout.open("./output/Gvector_" + std::to_string(ib) + "_rank" + std::to_string(my_rank) + ".out", std::ios::app);
		Gout << "\nNEW TIME STEP" << std::endl;
		for (i = 0; i < 3*iBody[ib].markers.size(); i++) {
			Gout << G[i] << std::endl;
		}
		Gout.close();
#endif



	// ******************* //
    // *** Newton Loop *** //
	// ******************* //

    size_t iter = 0;
    while (residual > tolerance && iter < max_terations) {
        iter++;
        res[0] = 0; //  not used res is used starting from 1
        res[1] = Tension0 - iBody[ib].tension[0] + (Fx[0] + Froude) * x[0] + Fy[0] * y[0] - G[1]; // SIMPLY SUPPORTED
        res[2] = pow(x[0] - x0, 2) + pow(y[0] - y0, 2) - G[2];
        // first two Lagrangian points: Equations (11) [Resolution of the filament problem]
        res[3] = iBody[ib].tension[0] * x0 - (iBody[ib].tension[0] + iBody[ib].tension[1] + beta) * x[0] + iBody[ib].tension[1] * x[1] - G[3];
        res[4] = iBody[ib].tension[0] * y0 - (iBody[ib].tension[0] + iBody[ib].tension[1] + beta) * y[0] + iBody[ib].tension[1] * y[1] - G[4];
        
        ii = 1;
        for (i = 1; i < n - 1; i++) {
            ii = ii + 3;
            res[ii + 1] = pow(x[i] - x[i - 1], 2) + pow(y[i] - y[i - 1], 2) - G[ii + 1];
            res[ii + 2] = iBody[ib].tension[i] * x[i - 1]-(iBody[ib].tension[i] + iBody[ib].tension[i + 1] + beta) * x[i] + iBody[ib].tension[i + 1] * x[i + 1] - G[ii + 2];
            res[ii + 3] = iBody[ib].tension[i] * y[i - 1]-(iBody[ib].tension[i] + iBody[ib].tension[i + 1] + beta) * y[i] + iBody[ib].tension[i + 1] * y[i + 1] - G[ii + 3];
        }
        res[3 * n - 1] = pow(x[n - 1] - x[n - 2], 2) + pow(y[n - 1] - y[n - 2], 2) - G[3 * n - 1];
        //  FREE EXTERMITY CONDITION
        res[3 * n] = 2 * iBody[ib].tension[n - 1] * x[n - 2]-(2 * iBody[ib].tension[n - 1] + beta) * x[n - 1] - G[3 * n];
        res[3 * n + 1] = 2 * iBody[ib].tension[n - 1] * y[n - 2]-(2 * iBody[ib].tension[n - 1] + beta) * y[n - 1] - G[3 * n + 1];

        // maximum residual evaluation
        residual = 0;
        for (size_t in = 1; in <= 3 * n + 1; ++in) {
            if (abs(res[in]) > residual) {
                residual = abs(res[in]);
            }
        }
        for (i = 0; i <= 3 * n + 1; ++i) {
            indx[i] = 0;
            for (size_t j = 0; j < m1 + m2 + 2; ++j) {
                AA[i][j] = 0.0;
            }
            for (size_t j = 0; j < m1 + 1; ++j) {
                AL[i][j] = 0.0;
            }
        }

        //-- for jacobian matrix in a compact form (only non zeros entries) 

        AA[1][4] = 1.0;                 // SIMPLY SUPPORTED
        AA[1][5] = -1.0;                // SIMPLY SUPPORTED
        AA[1][6] = Fx[0] + Froude;      // SIMPLY SUPPORTED
        AA[1][7] = Fy[0];               // SIMPLY SUPPORTED
        
        // non-extensibility
        AA[2][5] = 2 * (x[0] - x0); //J[1][2]=2*(x[0]-x0);
        AA[2][6] = 2 * (y[0] - y0); //J[1][3]=2*(y[0]-y0);

        AA[3][3] = x0 - x[0]; //J[2][1]=x0-x[0];
        AA[3][4] = -(iBody[ib].tension[0] + iBody[ib].tension[1] + beta); //J[2][2]=-(Tsol[0]+Tsol[1]+beta);
        AA[3][6] = x[1] - x[0]; //J[2][4]=x[1]-x[0];
        AA[3][7] = iBody[ib].tension[1]; //J[2][5]=Tsol[1];

        AA[4][2] = y0 - y[0]; //J[3][1]=y0-y[0];
        AA[4][4] = -(iBody[ib].tension[0] + iBody[ib].tension[1] + beta); //J[3][3]=-(Tsol[0]+Tsol[1]+beta); 
        AA[4][5] = y[1] - y[0]; //J[3][4]=y[1]-y[0];
        AA[4][7] = iBody[ib].tension[1]; //J[3][6]=Tsol[1]; 
        //
        ii = 1;
        for (size_t i = 1; i < n - 1; ++i) {
            ii = ii + 3;
            // inextensibility     
            AA[ii + 1][2] = 2 * (x[i] - x[i + 1]); //J[ii][ii-2]=2*(x[i]-x[i+1]);
            AA[ii + 1][3] = 2 * (y[i] - y[i + 1]); //J[ii][ii-1]=2*(y[i]-y[i+1]); 
            AA[ii + 1][5] = 2 * (x[i + 1] - x[i]); //J[ii][ii+1]=2*(x[i+1]-x[i]); 
            AA[ii + 1][6] = 2 * (y[i + 1] - y[i]); //J[ii][ii+2]=2*(y[i+1]-y[i]);
            // x-mom                                                                                                                                            
            AA[ii + 2][1] = iBody[ib].tension[i]; //J[ii+1][ii-2]=Tsol[i]; 
            AA[ii + 2][3] = x[i - 1] - x[i]; //J[ii+1][ii]=x[i-1]-x[i];
            AA[ii + 2][4] = -(iBody[ib].tension[i] + iBody[ib].tension[i + 1] + beta); //J[ii+1][ii+1]=-(Tsol[i]+Tsol[i+1]+beta);
            AA[ii + 2][6] = x[i + 1] - x[i]; //J[ii+1][ii+3]=x[i+1]-x[i]; 
            AA[ii + 2][7] = iBody[ib].tension[i + 1]; //J[ii+1][ii+4]=Tsol[i+1];
            // y-mom                                                                                                                                            
            AA[ii + 3][1] = iBody[ib].tension[i]; //J[ii+2][ii-1]=Tsol[i]; 
            AA[ii + 3][2] = y[i - 1] - y[i]; //J[ii+2][ii]=y[i-1]-y[i]; 
            AA[ii + 3][4] = -(iBody[ib].tension[i] + iBody[ib].tension[i + 1] + beta); //J[ii+2][ii+2]=-(Tsol[i]+Tsol[i+1]+beta); 
            AA[ii + 3][5] = y[i + 1] - y[i]; //J[ii+2][ii+3]=y[i+1]-y[i];
            AA[ii + 3][7] = iBody[ib].tension[i + 1]; //J[ii+2][ii+5]=Tsol[i+1];
        }
        AA[3 * n - 1][2] = 2 * (x[n - 2] - x[n - 1]); //J[3*n-2][3*n-4]=2*(x[n-2]-x[n-1]);
        AA[3 * n - 1][3] = 2 * (y[n - 2] - y[n - 1]); //J[3*n-2][3*n-3]=2*(y[n-2]-y[n-1]);
        AA[3 * n - 1][5] = 2 * (x[n - 1] - x[n - 2]); //J[3*n-2][3*n-1]=2*(x[n-1]-x[n-2]);
        AA[3 * n - 1][6] = 2 * (y[n - 1] - y[n - 2]); //J[3*n-2][3*n]  =2*(y[n-1]-y[n-2]);
        AA[3 * n][1] = 2 * iBody[ib].tension[n - 1]; //J[3*n-1][3*n-4]=2*Tsol[n-1]; 
        AA[3 * n][3] = 2 * x[n - 2] - 2 * x[n - 1]; //J[3*n-1][3*n-2]=2*x[n-2]-2*x[n-1];
        AA[3 * n][4] = -(2 * iBody[ib].tension[n - 1] + beta); //J[3*n-1][3*n-1]=-(2*Tsol[n-1]+beta); 
        AA[3 * n + 1][1] = 2 * iBody[ib].tension[n - 1]; //J[3*n][3*n-3]=2*Tsol[n-1];
        AA[3 * n + 1][2] = 2 * y[n - 2] - 2 * y[n - 1]; //J[3*n][3*n-2]=2*y[n-2]-2*y[n-1];
        AA[3 * n + 1][4] = -(2 * iBody[ib].tension[n - 1] + beta); //J[3*n][3*n]=-(2*Tsol[n-1]+beta); 
        //---------------------------------------------------------------------------------------- 


				
        ibm_bandec(AA, 3 * n + 1, m1, m2, AL, indx, d); // LU decomposition
        ibm_banbks(AA, 3 * n + 1, m1, m2, AL, indx, res); // solve banded problem



		// Iteration update 
        Tension0 = Tension0 - res[1];
        ii = 2;
        for (size_t i = 0; i < n; ++i) {
            iBody[ib].tension[i] -= res[ii];
            x[i] -= res[ii + 1];
            y[i] -= res[ii + 2];
            ii = ii + 3;
        }



    } // end of Newton while loop


	// ********************************* End of copy and paste ************************************** //
    





	/*
	Using the output from the above we update the marker positions and desired velocities
	*/

#ifdef IBM_DEBUG
		// DEBUG -- write out res vector
		std::ofstream resout;
		resout.open("./output/res_vector_" + std::to_string(ib) + "_rank" + std::to_string(my_rank) + ".out", std::ios::app);
		resout << "\nNEW TIME STEP" << std::endl;
		for (size_t i = 0; i < 3*iBody[ib].markers.size(); i++) {
			resout << res[i] << std::endl;
		}
		resout.close();
#endif


	// Clamped boundary condition
	if (iBody[ib].BCs[0] == 2) {

		// Apply clamped BC to first flexible marker
		x[0] = x0 + ds_nondim * cos(0.0 * PI / 180.0);
		y[0] = y0 + ds_nondim * sin(0.0 * PI / 180.0);

	}

	// After the Newton iteration has completed markers are updated from the temporary vectors used in the loop
    for (i = 0; i < iBody[ib].markers.size(); i++) {
        // Initial physical position (time = t) becomes t-1 position
        iBody[ib].markers[i].position_old[0] = iBody[ib].markers[i].position[0];
        iBody[ib].markers[i].position_old[1] = iBody[ib].markers[i].position[1];
        // Current physical position comes from the newly computed positions
		if (i != 0) { // New position vectors exclude the simply supported end and start from next node in so i-1
			// Convert filament-normalised coordinates back to lu then to physical spacing then add offset of simply supported end
			iBody[ib].markers[i].position[0] = (x[i-1] * length_lu * dx) + iBody[ib].markers[0].position[0];    
			iBody[ib].markers[i].position[1] = (y[i-1] * length_lu * dy) + iBody[ib].markers[0].position[1];
		}
    }

	// Desired velocity (in lattice units) on the makers is computed using a first order estimate based on position change
    for (i = 0; i < iBody[ib].markers.size(); i++) {
        iBody[ib].markers[i].desired_vel[0] = ( (iBody[ib].markers[i].position[0] - iBody[ib].markers[i].position_old[0]) / dx) / (1 / pow(2,level));
        iBody[ib].markers[i].desired_vel[1] = ( (iBody[ib].markers[i].position[1] - iBody[ib].markers[i].position_old[1]) / dy) / (1 / pow(2,level));
    }


}


// **************************************************************************************** //
/*
* Given an n by n band diagonal matrix A with 
* m1 subdiagonal rows and m2 superdiagonal rows,
* compactly stored in the array A[1..n][1..m1+m2+1] , this routine
* constructs an LU decomposition of a rowwise permutation of A. The upper
* triangular matrix replaces A, while the lower triangular matrix is
* returned in AL[1..n][1..m1]. INDX[1..n] is an output vector
* which records the row permutation effected by the partial
* pivoting; D is output as +/-1 depending on whether the number
* of row interchanges was even or odd, respectively.
* This routine is used in combination with ibm_banbks() to solve
* band-diagonal sets of equations. Once the matrix A has been 
* decomposed, any number of right-hand sides can be solved in
* turn by repeated calls to ibm_banbks()
* A		= array of subdiagonal and superdiagonals rows
* n		= size of the square matrix \f$ A \f$
* m1	= number of subdiagonal rows
* m2	= number of superdiagonal rows
* AL	= lower triangular matrix
* INDX	= row permutation vector
* D		= odd or even number of row interchages
*/
#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}
#define TINY 1.0e-20
void GridObj::ibm_bandec(double **a, unsigned long n, unsigned int m1, unsigned int m2, double **al,
	unsigned long indx[], double *d)
{
	unsigned long i,j,k,l;
	unsigned int mm;
	double dum;
    
	mm=m1+m2+1;
	l=m1;
	for (i=1;i<=m1;i++) {
		for (j=m1+2-i;j<=mm;j++) a[i][j-l]=a[i][j];
		l--;
		for (j=mm-l;j<=mm;j++) a[i][j]=0.0;
	}
	*d=1.0;
	l=m1;
	for (k=1;k<=n;k++) {
		dum=a[k][1];
		i=k;
		if (l < n) l++;
		for (j=k+1;j<=l;j++) {
			if (fabs(a[j][1]) > fabs(dum)) {
				dum=a[j][1];
				i=j;
			}
		}
		indx[k]=i;
		if (dum == 0.0) a[k][1]=TINY;
		if (i != k) {
			*d = -(*d);
			for (j=1;j<=mm;j++) SWAP(a[k][j],a[i][j])
		}
		for (i=k+1;i<=l;i++) {
			dum=a[i][1]/a[k][1];
			al[k][i-k]=dum;
			for (j=2;j<=mm;j++) a[i][j-1]=a[i][j]-dum*a[k][j];
			a[i][mm]=0.0;
		}
	}
   
        //gettimeofday(&t2, NULL);
        //elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
        //elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
        //printf("bandec time %e ms\n",elapsedTime);  
 
}
#undef SWAP
#undef TINY
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */


// **************************************************************************************** //
/*
* Given the arrays A, AL, and INDX as returned from ibm_bandec(), 
* and given a right-hand side vector B[1..n], solves the band diagonal
* linear equations AX = B. The solution vector X overwrites
* B. The other input arrays are not modified, and can be left
* in place for successive calls with different right-hand sides.
* A		= array of subdiagonal and superdiagonals rows
* n		= size of the square matrix A
* m1	= number of subdiagonal rows
* m2	= number of superdiagonal rows
* AL	= lower triangular matrix
* INDX	= row permutation vector
* B		= right hand side vector
*/
#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}
void GridObj::ibm_banbks(double **a, unsigned long n, unsigned int m1, unsigned int m2, double **al,
	unsigned long indx[], double b[])
{
	unsigned long i,k,l;
	unsigned int mm;
	double dum;
        //timeval t1, t2;
        //double elapsedTime;

        //gettimeofday(&t1, NULL);
    
	mm=m1+m2+1;
	l=m1;
	for (k=1;k<=n;k++) {
		i=indx[k];
		if (i != k) SWAP(b[k],b[i])
		if (l < n) l++;
		for (i=k+1;i<=l;i++) b[i] -= al[k][i-k]*b[k];
	}
	l=1;
	for (i=n;i>=1;i--) {
		dum=b[i];
		for (k=2;k<=l;k++) dum -= a[i][k]*b[k+i-1];
		b[i]=dum/a[i][1];
		if (l < mm) l++;
	}

}
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software ?421.1-9. */

// **************************************************************************************** //
// Routine to update the position of deformable bodies
void GridObj::ibm_position_update( unsigned int ib ) {

	// If a flexible body then launch structural solver to find new positions
	if (iBody[ib].flex_rigid) {

		// Do jacowire calculation
		ibm_jacowire(ib);

	} else {	// Body is deformable but not flexible so positional update comes from 
				// external forcing

		
		// Call some routine for external forcing here...


	}

	// Recompute support for new marker positions
	for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
					
		// Erase support vectors
		// Note:	Fixing the number of support sites to 9 negates 
		//			the need to erase the vectors and will improve speed.
		iBody[ib].markers[m].supp_i.clear();
		iBody[ib].markers[m].supp_j.clear();
		iBody[ib].markers[m].supp_k.clear();
		iBody[ib].markers[m].deltaval.clear();
					
		// Recompute support
		ibm_findsupport(ib, m);
	}


}
// **************************************************************************************** //

// Overloaded function to update position of body in a group using the position of the flexible
// body in the group. Must be called after all previous positional update routines have been called.
// group_update flag is only there to ensure the overload is defined.
void GridObj::ibm_position_update_grp( unsigned int group ) {

	// Find flexible body in group and store index
	unsigned int ib_flex;
	for (size_t i = 0; i < iBody.size(); i++) {

		if (iBody[i].flex_rigid && iBody[i].groupID == group) {
			ib_flex = i;
			break;
		}

	}

	// Loop over bodies in group
	for (size_t ib = 0; ib < iBody.size(); ib++) {

		// If body is deformable but not flexible give it a positional update from flexible
		if (iBody[ib].groupID == group && 
			iBody[ib].deformable && 
			iBody[ib].flex_rigid == false) {

			// Copy the position vectors of flexible markers in x and y directions
			for (size_t m = 0; m < iBody[ib].markers.size(); m++) {

				for (size_t d = 0; d < 2; d++) {
					iBody[ib].markers[m].position_old[d] = iBody[ib].markers[m].position[d];
					iBody[ib].markers[m].position[d] = iBody[ib_flex].markers[m].position[d];
				}

			}

			// Recompute support for new marker positions
			for (size_t m = 0; m < iBody[ib].markers.size(); m++) {
					
				// Erase support vectors
				// Note:	Fixing the number of support sites to 9 negates 
				//			the need to erase the vectors and will improve speed.
				iBody[ib].markers[m].supp_i.clear();
				iBody[ib].markers[m].supp_j.clear();
				iBody[ib].markers[m].supp_k.clear();
				iBody[ib].markers[m].deltaval.clear();
					
				// Recompute support
				ibm_findsupport(ib, m);
			}
		}
	}


}
// **************************************************************************************** //
