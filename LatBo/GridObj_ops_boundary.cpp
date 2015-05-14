/* This file holds all the code for applying standard LBM boundary conditions such as bounce-back, inlet and outlet. */

#include "stdafx.h"
#include "GridObj.h"
#include <vector>

// ***************************************************************************************************
/*
	Boundary condition application routine
	Supply an integer r indicating from which level the algorithm is to be executed plus a flag to specify
	which type of condition should be applied.
	0 == apply all simultaneously
	1 == apply solid site
	2 == apply inlet and outlet

	Boundary label types are:
	0 == solid site (no-slip)
	7 == inlet site
	8 == outlet site
*/
void GridObj::LBM_boundary (int bc_type_flag) {
/*
	// Get grid sizes
	size_t N_lim = Grids[r].XPos.size();
	size_t M_lim = Grids[r].YPos.size();
	size_t K_lim = Grids[r].ZPos.size();

	// Loop over grid, identify BC required & apply BC
	for (size_t i = 0; i < N_lim; i++) {
		for (size_t j = 0; j < M_lim; j++) {
			for (size_t k = 0; k < K_lim; k++) {
				
				// Declare vector to store missing populations
				std::vector<size_t> missing_pops;

				// Get index for current site
				size_t idx = idxmap(i,j,k,M_lim,K_lim);
*/
				/*	******************************************************************************************
					*************************************** ZOU - HE *****************************************
					***********************************  for solid sites *************************************
					******************************************************************************************
				*/
    /*
				if (Grids[r].LatTyp[idx] == 0 && (bc_type_flag == 1 || bc_type_flag == 0)) {
					
					// Look at each direction (excluding rest particle)
					for (size_t v = 0; v < nVels-1; v++) {

						size_t idx_v = idxmap(i,j,k,v,M_lim,N_lim,nVels);

						// Compute destination coordinates for stream (periodicity applied)
						size_t dest_x = (i+c[0][v] + N_lim) % N_lim;
						size_t dest_y = (j+c[1][v] + M_lim) % M_lim;
						size_t dest_z = (k+c[2][v] + K_lim) % K_lim;
						size_t idx_dest = idxmap(dest_x,dest_y,dest_z,M_lim,K_lim);
						
						// If destination not another solid site, population needs correcting
						if (Grids[r].LatTyp[idx_dest] != 0) {

							// Log that this direction needs adjustment
							missing_pops.push_back(v);

						}

					}

					// Based on population availability, decide where wall is and apply suitable Zou-He
					// calculation to fill in or correct missing populations
					applyZouHe(missing_pops, Grids[r].LatTyp[idx]);

*/
				/*	******************************************************************************************
					*************************************** ZOU - HE *****************************************
					*********************************** for inlet sites **************************************
					******************************************************************************************
				*/
    /*
				} else if (Grids[r].LatTyp[idx] == 7 && (bc_type_flag == 2 || bc_type_flag == 0)) {


					// DO SOMETHING
    */
				/*	******************************************************************************************
					************************************* EXTRAPOLATION **************************************
					************************************ for outlet sites ************************************
					******************************************************************************************
				*/
    /*
				} else if (Grids[r].LatTyp[idx] == 8 && (bc_type_flag == 2 || bc_type_flag == 0)) {

					// DO SOMETHING

				}
    */
}

// ***************************************************************************************************

// Routine to apply Zou-He boundary conditions
void applyZouHe( std::vector<size_t> missing_pops, int label ) {

	/* Zou-He velocity boundary condition computed from the following equations
	rho = sum ( fi )
	rho*ux = sum( fi * cxi )
	rho*uy = sum( fi * cyi )
	rho*uz = sum( fi * czi )
	(fi - feq)_in = (fi - feq)_out ------ normal to wall
	
	.... more for 3D

	3 populations (2D) or 5 populations (3D) will be unknown for the boundary site
	*/

/*
	// Get wall identification
	unsigned wall_type;
	wall_type = whichwall(missing_pops);

	// Set velocity vector based on label
	if (label == 0) {
		double u_wall[3] = {0.0, 0.0, 0.0};
	} else if (label == 7) {
		double u_wall[3] = {u_0x, u_0y, u_0z};
	}
	
	// Apply appropriate BCs -- switch might be slow depending on compiler options
	switch (wall_type)
	{
	case 1: // Left Wall





		
		break;

	case 2: // Right Wall

		
		break;

	case 3: // Top Wall

		
		break;

	case 4: // Bottom Wall

		
		break;

	case 5: // Front Wall

		
		break;

	case 6: // Back Wall

		
		break;

	// What do we do about corners?

	}
*/
}

// ***************************************************************************************************

// Routine to compute the opposite direction fo the one supplied based on D2Q9 or D3Q19 numbering
int getOpposite(int direction) {
    /*

	int direction_opposite;

	// If rest particle then opposite is simply itself
	if (direction == nVels-1) {
		
		direction_opposite = direction;

	} else if (dims == 3) {
		
		// Using D3Q19
		direction_opposite = direction + 1;

	} else {

		// Using D2Q9
		direction_opposite = (direction + 4 % 8);
	}

	return direction_opposite;
    */
    return 0;

}

// ***************************************************************************************************

// Routine to identify which wall the boundary site sits in based on missing populations
unsigned whichwall( std::vector<size_t> pops ) {

	// Declare unsigned integer type
	unsigned wall_type = 0;

	/*
	Key:
	1 == Left
	2 == Right
	3 == Top
	4 == Bottom
	5 == Front
	6 == Back

	Produces an error when called if the boundary condition has not been implemented yet
	*/

	// Expecting 3 missing populations (2D) or 5 (3D)
	// Knowledge of normal population will tell us which boundary it is.


	return wall_type;
}

/*





            
            
                    
                    // Incoming direction from source
                    kopp = 1 + mod(k+3,8);
                    
            
        % INLET CONDITIONS
        elseif Grid(r+1).LatTyp(j,i) == 7 && type_flag == 2
            
            % Implement using 4 equations Zou-He
            % rho_in = sum( fi )
            % rho_in * ux = (f1 + f2 + f8) - (f4 + f5 + f6)
            % rho_in * uy = (f2 + f3 + f4) - (f6 + f7 + f8)
            % f1 - feq1 = f5 - feq5 (equilibrium normal to boundary)
            
            % Compute desired u at the inlet site
            [ux,uy] = InletVel(Grid(r+1).YPos(j),u_in,nu,rho_in,r);
         
            % Find density on wall corresponding to given velocity
            rho_w = (1 / (1 - ux)) * (...
                Grid(r+1).f(j,i,9) + ...
                Grid(r+1).f(j,i,3) + ...
                Grid(r+1).f(j,i,7) + ...
                2 * (...
                Grid(r+1).f(j,i,4) + ...
                Grid(r+1).f(j,i,5) + ...
                Grid(r+1).f(j,i,6)...
            ) );
            
            % Find f1 using equations above
            Grid(r+1).f(j,i,1) = Grid(r+1).f(j,i,5) + (2/3) * rho_w * ux;
            
            % Find f2 using equations above
            Grid(r+1).f(j,i,2) = 0.5 * ( (rho_w * ux) - ...
                (Grid(r+1).f(j,i,1) + Grid(r+1).f(j,i,3)) + ...
                Grid(r+1).f(j,i,5) + ...
                2*Grid(r+1).f(j,i,6) + ...
                Grid(r+1).f(j,i,7) );
            
            % Find f8 using equations above
            Grid(r+1).f(j,i,8) = 0.5 * ( (rho_w * ux) - ...
                (Grid(r+1).f(j,i,1) + Grid(r+1).f(j,i,7)) + ...
                Grid(r+1).f(j,i,3) + ...
                2*Grid(r+1).f(j,i,4) + ...
                Grid(r+1).f(j,i,5) );
            
            
        % OUTLET CONDITIONS
        elseif Grid(r+1).LatTyp(j,i) == 8 && type_flag == 2
            
%             % OPTION1: Use second order polynomial extrapolation
%             for k =  [4 5 6]
%                 
%                 Grid(r+1).f(j,i,k) = 2 * Grid(r+1).f(j,i-1,k) - Grid(r+1).f(j,i-2,k);
%                 
%             end
            
            
            % OPTION 2: Use linear extrapolation (more stable than OPT 1)
            for k =  [4 5 6]
                
                y2 = Grid(r+1).f(j,i-1,k);
                y1 = Grid(r+1).f(j,i-2,k);
                x1 = 0;
                x2 = Grid(r+1).dx;
                x3 = 2*x2;
                
                lin_m = (y2 - y1) / (x2 - x1);
                lin_c = y1;
                
                Grid(r+1).f(j,i,k) = lin_m * x3 + lin_c;
                
            end            
            
        end
        
    end
end

end

	

    */

// ***************************************************************************************************