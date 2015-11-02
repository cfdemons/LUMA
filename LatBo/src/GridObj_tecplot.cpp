// Routine for doing a TecPlot write out to file

#include "../inc/stdafx.h"
#include "../inc/GridObj.h"


// ************************************************************** //
void GridObj::io_tecplot(double tval) {

	std::ofstream tecfile;
	// Get grid sizes
	unsigned int N_lim, M_lim, K_lim, min, minz;
	N_lim = XInd.size();
	M_lim = YInd.size();
	K_lim = ZInd.size();
	min = 0; minz = min;

	// If building for MPI then correct grid sizes to avoid writing out outer overlap
	if (level == 0) {

#ifdef BUILD_FOR_MPI
		min = 1; minz = min;
		N_lim = XInd.size()-1;
		M_lim = YInd.size()-1;
#if (dims == 3)
		K_lim = ZInd.size()-1;
#else
		K_lim = 1;
		minz = 0;
#endif
#endif

	}

	// If the first rank or a subgrid then create file anew and add header details
	if (my_rank == 0 || level != 0) {
		tecfile.open("./output/tecplotout.Lev" + std::to_string(level) + ".Reg" + std::to_string(region_number)
			+ "." + std::to_string((int)tval) + ".dat", std::ios::out);

		// Add header
		tecfile << "TITLE = All grid quantities" << std::endl;
		tecfile << "FILETYPE = FULL" << std::endl;
		tecfile << "VARIABLES = \"X\" \"Y\" \"Z\" \"RHO\" \"UX\" \"UY\" \"UZ\" \"TA_RHO\" \"TA_UX\" \"TA_UY\" \"TA_UZ\" "<<
			"\"TA_UXUX\" \"TA_UXUY\" \"TA_UXUZ\" \"TA_UYUY\" \"TA_UYUZ\" \"TA_UZUZ\"" << std::endl;
		tecfile << "ZONE" << std::endl;

		tecfile << "I = " << std::to_string(N)
				<< ", J = " << std::to_string(M)
				<< ", K = " <<
#if (dims == 3)
				std::to_string(K)
#else
				std::to_string(1)
#endif				
				<< std::endl;

		tecfile << "ZONETYPE = Ordered, DATAPACKING = POINT" << std::endl;
		tecfile << "SOLUTIONTIME = " << std::to_string(tval) << std::endl;



	} else {
		// If not the first rank, append data
		tecfile.open("./output/tecplotout.Lev" + std::to_string(level) + ".Reg" + std::to_string(region_number)
			+ "." + std::to_string((int)tval) + ".dat", std::ios::out|std::ios::app);
	}


	// Counters
	unsigned int i,j,k,v;
		
	// Write out values
	for (k = minz; k < K_lim; k++) {
		for (j = min; j < M_lim; j++) {
			for (i = min; i < N_lim; i++) {
				
				// Write out X, Y, Z
				tecfile << XPos[i] << "\t" << YPos[j] << "\t" << ZPos[k] << "\t";

				// Write out rho and u
				tecfile << rho(i,j,k,YInd.size(),ZInd.size()) << "\t";
				for (v = 0; v < dims; v++) {
					tecfile << u(i,j,k,v,YInd.size(),ZInd.size(),dims) << "\t";
				}
#if (dims != 3)
				tecfile << 0 << "\t";
#endif
				
				// Write out time averaged rho and u
				tecfile << rho_timeav(i,j,k,YInd.size(),ZInd.size()) << "\t";
				for (v = 0; v < dims; v++) {
					tecfile << ui_timeav(i,j,k,v,YInd.size(),ZInd.size(),dims) << "\t";
				}
#if (dims != 3)
				tecfile << 0 << "\t";
#endif

				// Write out time averaged u products
				tecfile << uiuj_timeav(i,j,k,0,YInd.size(),ZInd.size(),(3*dims-3)) << "\t";
				tecfile << uiuj_timeav(i,j,k,1,YInd.size(),ZInd.size(),(3*dims-3)) << "\t";
#if (dims == 3)
				tecfile << uiuj_timeav(i,j,k,2,YInd.size(),ZInd.size(),(3*dims-3)) << "\t";
#else
				tecfile << 0 << "\t";
#endif
#if (dims == 3)
				tecfile << uiuj_timeav(i,j,k,3,YInd.size(),ZInd.size(),(3*dims-3)) << "\t";
				tecfile << uiuj_timeav(i,j,k,4,YInd.size(),ZInd.size(),(3*dims-3)) << "\t";
				tecfile << uiuj_timeav(i,j,k,5,YInd.size(),ZInd.size(),(3*dims-3)) << "\t";
#else
				tecfile << uiuj_timeav(i,j,k,2,YInd.size(),ZInd.size(),(3*dims-3)) << "\t";
				tecfile << 0 << "\t" << 0 << "\t";
#endif

				// New line
				tecfile << std::endl;


			}
		}
	}


}