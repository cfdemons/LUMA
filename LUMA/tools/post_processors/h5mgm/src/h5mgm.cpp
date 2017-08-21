/*
* --------------------------------------------------------------
*
* ------ Lattice Boltzmann @ The University of Manchester ------
*
* -------------------------- L-U-M-A ---------------------------
*
*  Copyright (C) The University of Manchester 2017
*  E-mail contact: info@luma.manchester.ac.uk
*
* This software is for academic use only and not available for
* further distribution commericially or otherwise without written consent.
*
*/

#include "VelocitySorter.h"
#include "ThirdParty/vtkCleanUnstructuredGrid.h"

// Method to construct a 3D polyhedron out of faces and add to an unstructured grid
size_t addCell(vtkSmartPointer<vtkPoints> global_pts,
	vtkSmartPointer<vtkUnstructuredGrid> grid,
	std::vector<vtkIdType>& global_ids,
	std::vector< std::vector<double> >& cell_pts,
	int *point_count, int dimensions_p)
{

	// Local declarations
	vtkSmartPointer<vtkCellArray> face_array = vtkSmartPointer<vtkCellArray>::New();
	std::vector<vtkIdType> face_ids;
	size_t status = 0;

	// Add new points to global points and assign new id
	for (size_t c = 0; c < cell_pts.size(); c++) {
		global_pts->InsertNextPoint(cell_pts[c][0], cell_pts[c][1], cell_pts[c][2]);
		global_ids.push_back(*point_count + c);
	}

	if (dimensions_p == 3)
	{

		// Specify points that make up voxel in specific order
		face_ids.push_back(*point_count + 0);
		face_ids.push_back(*point_count + 4);
		face_ids.push_back(*point_count + 2);
		face_ids.push_back(*point_count + 6);
		face_ids.push_back(*point_count + 1);
		face_ids.push_back(*point_count + 5);
		face_ids.push_back(*point_count + 3);
		face_ids.push_back(*point_count + 7);
		

		// Update offset
		*point_count += 8;
	}
	else
	{

		// Specify points that make up quad in specific order
		face_ids.push_back(*point_count + 0);
		face_ids.push_back(*point_count + 2);
		face_ids.push_back(*point_count + 1);
		face_ids.push_back(*point_count + 3);

		// Update offset
		*point_count += 4;

	}

	// Update grid with new point list and add cell
	grid->SetPoints(global_pts);

	if (dimensions_p == 3) {
		grid->InsertNextCell(VTK_VOXEL, 8, &face_ids[0]);
	}
	else {
		grid->InsertNextCell(VTK_PIXEL, 4, &face_ids[0]);

	}

	// Free face array (forces destructor to be called)
	face_array = NULL;

	return status;

}

/* H5 Multi-Grid Merge Tool for post-processing HDF5 files written by LUMA */
int main(int argc, char* argv[])
{

	// Parse arguments and handle
	std::string case_num("000");
	for (int a = 1; a < argc; ++a)
	{
		std::string arg_str = std::string(argv[a]);

		if (arg_str == "version")
		{
			std::cout << "H5MultiGridMerge (h5mgm) Version " << H5MGM_VERSION << std::endl;
			return 0;
		}
		else if (arg_str == "quiet")
		{
			bQuiet = true;
		}
		else if (arg_str == "loud")
		{
			bLoud = true;
		}
		else if (arg_str == "cut")
		{
			bCutSolid = true;
		}
		else if (arg_str == "legacy")
		{
			bLegacy = true;
		}
		else if (arg_str == "sorter")
		{
			bSorter = true;
		}
		else
		{
			case_num = std::string(argv[a]);
		}
	}

	// Print out to screen
	std::cout << "H5MultiGridMerge (h5mgm) Version " << H5MGM_VERSION << ". Running..." << std::endl;

	// Path for output
	std::string path_str(H5MGM_OUTPUT_PATH);

	// Create directory
	std::string command = "mkdir -p " + path_str;
#ifdef _WIN32   // Running on Windows
	CreateDirectoryA((LPCSTR)path_str.c_str(), NULL);
#else   // Running on Unix system
	system(command.c_str());
#endif // _WIN32

	// Open log file if not set to quiet
	if (bQuiet == false)
	{
		std::string logpath = H5MGM_OUTPUT_PATH;
		logpath += "/h5mgm.log";
		logfile.open(logpath, std::ios::out | std::ios::app);
		logfile << "---------------------------------------" << std::endl;
	}

	// Start process
	std::cout << "Reconstructing HDF data..." << std::endl;

	// Turn auto error printing off
	H5Eset_auto(H5E_DEFAULT, NULL, NULL);


	// If using the sorter use the standalone class then return
	if (bSorter)
	{
		VelocitySorter<double> *vs = new VelocitySorter<double>();
		vs->readAndSort();
		delete vs;
		return 0;
	}



	// Construct L0 filename
	std::string IN_FILE_NAME("./hdf_R0N0.h5");

	// Set T = 0 time string
	std::string TIME_STRING = "/Time_0";

	// Declarations
	herr_t status = 0;
	hid_t output_fid = NULL;
	hid_t input_fid = NULL;
	hid_t input_aid = NULL;
	int dimensions_p, levels, regions, timesteps, out_every, mpi_flag;
	double dx;
	std::string VAR;
	int* LatTyp = nullptr;
	double *dummy_d = nullptr;
	int *dummy_i = nullptr;

	// Open L0 input file
	input_fid = H5Fopen(IN_FILE_NAME.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (input_fid == NULL) writeInfo("Cannot open input file!", eHDF);

	// Read in key attributes from the file
	input_aid = H5Aopen(input_fid, "Dimensions", H5P_DEFAULT);
	if (input_aid == NULL) writeInfo("Cannot open attribute!", eHDF);
	status = H5Aread(input_aid, H5T_NATIVE_INT, &dimensions_p);
	if (status != 0) writeInfo("Cannot read attribute!", eHDF);
	status = H5Aclose(input_aid);
	if (status != 0) writeInfo("Cannot close attribute!", eHDF);

	input_aid = H5Aopen(input_fid, "NumberOfGrids", H5P_DEFAULT);
	if (input_aid <= 0) writeInfo("Cannot open attribute!", eHDF);
	status = H5Aread(input_aid, H5T_NATIVE_INT, &levels);
	if (status != 0) writeInfo("Cannot read attribute!", eHDF);
	status = H5Aclose(input_aid);
	if (status != 0) writeInfo("Cannot close attribute!", eHDF);

	input_aid = H5Aopen(input_fid, "NumberOfRegions", H5P_DEFAULT);
	if (input_aid <= 0) writeInfo("Cannot open attribute!", eHDF);
	status = H5Aread(input_aid, H5T_NATIVE_INT, &regions);
	if (status != 0) writeInfo("Cannot read attribute!", eHDF);
	status = H5Aclose(input_aid);
	if (status != 0) writeInfo("Cannot close attribute!", eHDF);

	input_aid = H5Aopen(input_fid, "Timesteps", H5P_DEFAULT);
	if (input_aid <= 0) writeInfo("Cannot open attribute!", eHDF);
	status = H5Aread(input_aid, H5T_NATIVE_INT, &timesteps);
	if (status != 0) writeInfo("Cannot read attribute!", eHDF);
	status = H5Aclose(input_aid);
	if (status != 0) writeInfo("Cannot close attribute!", eHDF);

	input_aid = H5Aopen(input_fid, "OutputFrequency", H5P_DEFAULT);
	if (input_aid <= 0) writeInfo("Cannot open attribute!", eHDF);
	status = H5Aread(input_aid, H5T_NATIVE_INT, &out_every);
	if (status != 0) writeInfo("Cannot read attribute!", eHDF);
	status = H5Aclose(input_aid);
	if (status != 0) writeInfo("Cannot close attribute!", eHDF);

	input_aid = H5Aopen(input_fid, "Mpi", H5P_DEFAULT);
	if (input_aid == NULL) writeInfo("Cannot open attribute!", eHDF);
	status = H5Aread(input_aid, H5T_NATIVE_INT, &mpi_flag);
	if (status != 0) writeInfo("Cannot read attribute!", eHDF);
	status = H5Aclose(input_aid);
	if (status != 0) writeInfo("Cannot close attribute!", eHDF);

	// Store basic data
	int *gridsize = (int*)malloc(3 * sizeof(int));	// Space for grid size
	gridsize[2] = 1;	// Set 3D dimension to 1, will get overwritten if actually 3D
	int num_grids = ((levels - 1) * regions + 1);	// Total number of grids

	// Close file
	status = H5Fclose(input_fid);
	if (status != 0) writeInfo("Cannot close file!", eHDF);

	// Create cell node position list for VTK
	std::vector< std::vector<double> > cell_points;
	if (dimensions_p == 3) {
		cell_points.resize(8, std::vector<double>(3));
	}
	else {
		cell_points.resize(4, std::vector<double>(3));
	}


		

	// Create VTK grid
	vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
		vtkSmartPointer<vtkUnstructuredGrid>::New();

	// Create VTK grid points
	vtkSmartPointer<vtkPoints> points =
		vtkSmartPointer<vtkPoints>::New();

	// Create point ID list
	std::vector<vtkIdType> *pointIds = 
		new std::vector<vtkIdType>();

	// Total point count
	int point_count = 0;
		
	// Status indicator
	size_t res = 0;

	std::cout << "Building mesh..." << std::endl;

	// Loop 1 to get mesh details
	for (int lev = 0; lev < levels; lev++) {
		for (int reg = 0; reg < regions; reg++) {

			// L0 doesn't have different regions
			if (lev == 0 && reg != 0) continue;

			// Debug
			int local_point_count = 0;
			std::cout << "Adding cells from L" << lev << " R" << reg << "..." << std::endl;			

			// Construct input file name
			std::string IN_FILE_NAME("./hdf_R" + std::to_string(reg) + "N" + std::to_string(lev) + ".h5");
			input_fid = H5Fopen(IN_FILE_NAME.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
			if (input_fid <= 0) writeInfo("Cannot open input file!", eHDF);

			// Get local grid size
			input_aid = H5Aopen(input_fid, "GridSize", H5P_DEFAULT);
			if (input_aid <= 0) writeInfo("Cannot open attribute!", eHDF);
			status = H5Aread(input_aid, H5T_NATIVE_INT, gridsize);
			if (status != 0) writeInfo("Cannot read attribute!", eHDF);
			status = H5Aclose(input_aid);
			if (status != 0) writeInfo("Cannot close attribute!", eHDF);

			// Get local dx
			input_aid = H5Aopen(input_fid, "Dx", H5P_DEFAULT);
			if (input_aid == NULL) writeInfo("Cannot open attribute!", eHDF);
			status = H5Aread(input_aid, H5T_NATIVE_DOUBLE, &dx);
			if (status != 0) writeInfo("Cannot read attribute!", eHDF);
			status = H5Aclose(input_aid);
			if (status != 0) writeInfo("Cannot close attribute!", eHDF);

			// Allocate space for LatTyp data
			int *Type = (int*)malloc((gridsize[0] * gridsize[1] * gridsize[2]) * sizeof(int));
			if (Type == NULL) {
				writeInfo("Not enough Memory!!!!", eFatal);
				exit(EXIT_FAILURE);
			}

			// Allocate space for X, Y, Z coordinates
			double *X = (double*)malloc((gridsize[0] * gridsize[1] * gridsize[2]) * sizeof(double));
			if (X == NULL) {
				writeInfo("Not enough Memory!!!!", eFatal);
				exit(EXIT_FAILURE);
			}
			double *Y = (double*)malloc((gridsize[0] * gridsize[1] * gridsize[2]) * sizeof(double));
			if (Y == NULL) {
				writeInfo("Not enough Memory!!!!", eFatal);
				exit(EXIT_FAILURE);
			}
			double *Z = (double*)malloc((gridsize[0] * gridsize[1] * gridsize[2]) * sizeof(double));
			if (Z == NULL) {
				writeInfo("Not enough Memory!!!!", eFatal);
				exit(EXIT_FAILURE);
			}

			// Create input dataspace
			hsize_t dims_input[1];
			dims_input[0] = gridsize[0] * gridsize[1] * gridsize[2];
			hid_t input_sid = H5Screate_simple(1, dims_input, NULL);
			if (input_sid <= 0) writeInfo("Cannot create input dataspace!", eHDF);

			// Open, read into buffers information required for mesh definition and close input datasets
			status = readDataset("/LatTyp", TIME_STRING, input_fid, input_sid, H5T_NATIVE_INT, Type);
			if (status != 0)
			{
				writeInfo("Typing matrix read failed -- exiting early.", eFatal);
				exit(EARLY_EXIT);
			}
			status = readDataset("/XPos", TIME_STRING, input_fid, input_sid, H5T_NATIVE_DOUBLE, X);
			if (status != 0)
			{
				writeInfo("X position vector read failed -- exiting early.", eFatal);
				exit(EARLY_EXIT);
			}
			status = readDataset("/YPos", TIME_STRING, input_fid, input_sid, H5T_NATIVE_DOUBLE, Y);
			if (status != 0)
			{
				writeInfo("Y position vector read failed -- exiting early.", eFatal);
				exit(EARLY_EXIT);
			}
			if (dimensions_p == 3) {
				status = readDataset("/ZPos", TIME_STRING, input_fid, input_sid, H5T_NATIVE_DOUBLE, Z);
				if (status != 0)
				{
					writeInfo("Z position vector read failed -- exiting early.", eFatal);
					exit(EARLY_EXIT);
				}
			}

			// Loop over grid sites
			for (int c = 0; c < gridsize[0] * gridsize[1] * gridsize[2]; c++) {

				if (bLoud || c % 100000 == 0) {
					std::cout << "\r" << "Examining cell " << std::to_string(c) << "/" <<
						std::to_string(gridsize[0] * gridsize[1] * gridsize[2]) << std::flush;
				}

				// Ignore list
				if (isOnIgnoreList(static_cast<eType>(Type[c]))) continue;

				local_point_count++;

				if (dimensions_p == 3) {

					// Find positions of corner points and store
					for (int p = 0; p < 8; ++p) {
						cell_points[p][0] = X[c] + e[0][p] * (dx / 2);
						cell_points[p][1] = Y[c] + e[1][p] * (dx / 2);
						cell_points[p][2] = Z[c] + e[2][p] * (dx / 2);
					}

				}
				else {

					// Find positions of corner points and store
					for (int p = 0; p < 4; ++p) {
						cell_points[p][0] = X[c] + e2[0][p] * (dx / 2);
						cell_points[p][1] = Y[c] + e2[1][p] * (dx / 2);
						cell_points[p][2] = 0.0;
					}

				}

				// Build cell from points and add to mesh
				res = addCell(points, unstructuredGrid, *pointIds, cell_points, &point_count, dimensions_p);

			}
			std::cout << "\r" << "Examining cell " << std::to_string(gridsize[0] * gridsize[1] * gridsize[2]) << "/" <<
				std::to_string(gridsize[0] * gridsize[1] * gridsize[2]) << std::flush;
			std::cout << std::endl;

			// Close input dataspace
			status = H5Sclose(input_sid);
			if (status != 0) writeInfo("Cannot close input dataspace!", eHDF);

			// Close input file
			status = H5Fclose(input_fid);
			if (status != 0) writeInfo("Cannot close input file!", eHDF);

			// Free memory
			free(Type);
			free(X);
			free(Y);
			free(Z);

			// Debug
			std::cout << "Valid Point Count  = " << local_point_count << std::endl;

		}
	}

	// Delete pointIds as we are finished with them
	delete pointIds;

	// Debug
	std::cout << "Total number of cells retained for merged mesh = " << unstructuredGrid->GetNumberOfCells() << std::endl;
	std::cout << "Adding data for each time step to mesh..." << std::endl;

	// Clean mesh to remove duplicate points (Paraview not VTK)
	vtkSmartPointer<vtkCleanUnstructuredGrid> cleaner =
		vtkSmartPointer<vtkCleanUnstructuredGrid>::New();
	cleaner->SetInputData(unstructuredGrid);
	cleaner->Update();

	// Replace the unstructured grid with the cleaned version
	unstructuredGrid->Reset();
	unstructuredGrid = cleaner->GetOutput();

	// Time loop and data addition
	int count = 0;
	for (size_t t = 0; t <= (size_t)timesteps; t += out_every)
	{

		// Create filename
		std::string vtkFilename = path_str + "/luma_" + case_num + "." + std::to_string(t);
		if (bLegacy)
			vtkFilename += ".vtk";
		else
			vtkFilename += ".vtu";

		// Add data
		vtkSmartPointer<vtkIntArray> LatTyp = vtkSmartPointer<vtkIntArray>::New();
		LatTyp->Allocate(unstructuredGrid->GetNumberOfCells());
		LatTyp->SetName("LatTyp");
		status = addDataToGrid("/LatTyp", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_i, H5T_NATIVE_INT, LatTyp);
		
		// If no typing matrix then assume the time step is not available and exit
		if (status == DATASET_READ_FAIL)
		{
			writeInfo("Couldn't find time step " + std::to_string(t) + ". Read failed -- exiting early.", eFatal);
			exit(EARLY_EXIT);
		}

		// MPI block data always read from Time_0
		if (mpi_flag)
		{
			vtkSmartPointer<vtkIntArray> Block = vtkSmartPointer<vtkIntArray>::New();
			Block->Allocate(unstructuredGrid->GetNumberOfCells());
			Block->SetName("MpiBlockNumber");
			status = addDataToGrid("/MpiBlock", "/Time_0", levels, regions, gridsize, unstructuredGrid, dummy_i, H5T_NATIVE_INT, Block);
		}

		vtkSmartPointer<vtkDoubleArray> RhoArray = vtkSmartPointer<vtkDoubleArray>::New();
		RhoArray->Allocate(unstructuredGrid->GetNumberOfCells());
		RhoArray->SetName("Rho");
		status = addDataToGrid("/Rho", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, RhoArray);

		vtkSmartPointer<vtkDoubleArray> Rho_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
		Rho_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
		Rho_TimeAv->SetName("Rho_TimeAv");
		status = addDataToGrid("/Rho_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, Rho_TimeAv);

		vtkSmartPointer<vtkDoubleArray> Ux = vtkSmartPointer<vtkDoubleArray>::New();
		Ux->Allocate(unstructuredGrid->GetNumberOfCells());
		Ux->SetName("Ux");
		status = addDataToGrid("/Ux", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, Ux);

		vtkSmartPointer<vtkDoubleArray> Uy = vtkSmartPointer<vtkDoubleArray>::New();
		Uy->Allocate(unstructuredGrid->GetNumberOfCells());
		Uy->SetName("Uy");
		status = addDataToGrid("/Uy", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, Uy);

		vtkSmartPointer<vtkDoubleArray> Ux_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
		Ux_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
		Ux_TimeAv->SetName("Ux_TimeAv");
		status = addDataToGrid("/Ux_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, Ux_TimeAv);

		vtkSmartPointer<vtkDoubleArray> Uy_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
		Uy_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
		Uy_TimeAv->SetName("Uy_TimeAv");
		status = addDataToGrid("/Uy_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, Uy_TimeAv);

		vtkSmartPointer<vtkDoubleArray> UxUx_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
		UxUx_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
		UxUx_TimeAv->SetName("UxUx_TimeAv");
		status = addDataToGrid("/UxUx_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, UxUx_TimeAv);

		vtkSmartPointer<vtkDoubleArray> UxUy_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
		UxUy_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
		UxUy_TimeAv->SetName("UxUy_TimeAv");
		status = addDataToGrid("/UxUy_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, UxUy_TimeAv);

		vtkSmartPointer<vtkDoubleArray> UyUy_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
		UyUy_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
		UyUy_TimeAv->SetName("UyUy_TimeAv");
		status = addDataToGrid("/UyUy_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, UyUy_TimeAv);

		if (dimensions_p == 3)
		{

			vtkSmartPointer<vtkDoubleArray> Uz = vtkSmartPointer<vtkDoubleArray>::New();
			Uz->Allocate(unstructuredGrid->GetNumberOfCells());
			Uz->SetName("Uz");
			status = addDataToGrid("/Uz", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, Uz);

			vtkSmartPointer<vtkDoubleArray> Uz_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
			Uz_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
			Uz_TimeAv->SetName("Uz_TimeAv");
			status = addDataToGrid("/Uz_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, Uz_TimeAv);

			vtkSmartPointer<vtkDoubleArray> UxUz_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
			UxUz_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
			UxUz_TimeAv->SetName("UxUz_TimeAv");
			status = addDataToGrid("/UxUz_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, UxUz_TimeAv);

			vtkSmartPointer<vtkDoubleArray> UyUz_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
			UyUz_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
			UyUz_TimeAv->SetName("UyUz_TimeAv");
			status = addDataToGrid("/UyUz_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, UyUz_TimeAv);

			vtkSmartPointer<vtkDoubleArray> UzUz_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
			UzUz_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
			UzUz_TimeAv->SetName("UzUz_TimeAv");
			status = addDataToGrid("/UzUz_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, UzUz_TimeAv);
		}

		// Otherwise, must be simply a missing dataset so continue to next time step
		TIME_STRING = "/Time_" + std::to_string(t + out_every);

		// Print progress to screen
		std::cout << "\r" << std::to_string((int)(((float)(t + out_every) / 
			(float)(timesteps + out_every)) * 100.0f)) << "% complete." << std::flush;

		// Write grid to file
		void *writer = nullptr;
		if (bLegacy)
		{
			vtkSmartPointer<vtkUnstructuredGridWriter> writer =
				vtkSmartPointer<vtkUnstructuredGridWriter>::New();
			writer->SetFileName(vtkFilename.c_str());
			writer->SetFileTypeToBinary();
			writer->SetFileName(vtkFilename.c_str());
#if VTK_MAJOR_VERSION <= 5
			writer->SetInput(unstructuredGrid);
#else
			writer->SetInputData(unstructuredGrid);
#endif
			writer->Write();
		}
		else
		{
			vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
				vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
			writer->SetFileName(vtkFilename.c_str());
#if VTK_MAJOR_VERSION <= 5
			writer->SetInput(unstructuredGrid);
#else
			writer->SetInputData(unstructuredGrid);
#endif
			writer->Write();
		}

		// Free the memory used (forces destructor call)
		LatTyp = NULL;
		RhoArray = NULL;
		Rho_TimeAv = NULL;
		Ux = NULL;
		Uy = NULL;
		Ux_TimeAv = NULL;
		Uy_TimeAv = NULL;
		UxUx_TimeAv = NULL;
		UxUy_TimeAv = NULL;
		UyUy_TimeAv = NULL;
	}

	return 0;
}

