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

#include "h5mgm.h"

// Method to construct a 3D polyhedron out of faces and add to an unstructured grid
size_t addCell(vtkSmartPointer<vtkPoints> global_pts,
	vtkSmartPointer<vtkUnstructuredGrid> grid,
	std::vector<vtkIdType>& global_ids,
	std::vector< std::vector<double> >& cell_pts,
	int *point_count, int dimensions_p)
{

	// Local declarations
	vtkSmartPointer<vtkCellArray> face_array =
		vtkSmartPointer<vtkCellArray>::New();
	std::vector<vtkIdType> face_ids;
	size_t status = 0;

	// Add new points to global points and assign new id
	for (size_t c = 0; c < cell_pts.size(); c++) {
		global_pts->InsertNextPoint(cell_pts[c][0], cell_pts[c][1], cell_pts[c][2]);
		global_ids.push_back(*point_count + c);
	}

	/* In 3D we construct a polyhedron out of faces.
	 * In 2D we simply construct a single polygon. */

	if (dimensions_p == 3) {

		// Construct one face at a time out of vertices
		for (int c = 0; c < 27; c++) {

			// If it is a face
			if (abs(e[0][c]) + abs(e[1][c]) + abs(e[2][c]) == 1) {

				// Get face centre vector
				int face_centre[3];
				face_centre[0] = e[0][c];
				face_centre[1] = e[1][c];
				face_centre[2] = e[2][c];

				// Clear temporary storage
				face_ids.clear();

				// Left face
				if (face_centre[0] == -1) {

					face_ids.push_back(*point_count + 0);
					face_ids.push_back(*point_count + 2);
					face_ids.push_back(*point_count + 8);
					face_ids.push_back(*point_count + 6);

				}
				// Right face
				else if (face_centre[0] == 1) {

					face_ids.push_back(*point_count + 18);
					face_ids.push_back(*point_count + 20);
					face_ids.push_back(*point_count + 26);
					face_ids.push_back(*point_count + 24);


				}
				// Bottom face
				else if (face_centre[1] == -1) {

					face_ids.push_back(*point_count + 0);
					face_ids.push_back(*point_count + 2);
					face_ids.push_back(*point_count + 20);
					face_ids.push_back(*point_count + 18);

				}
				// Top face
				else if (face_centre[1] == 1) {

					face_ids.push_back(*point_count + 6);
					face_ids.push_back(*point_count + 8);
					face_ids.push_back(*point_count + 26);
					face_ids.push_back(*point_count + 24);

				}
				// Front face
				else if (face_centre[2] == -1) {

					face_ids.push_back(*point_count + 0);
					face_ids.push_back(*point_count + 6);
					face_ids.push_back(*point_count + 24);
					face_ids.push_back(*point_count + 18);

				}
				// Back face
				else if (face_centre[2] == 1) {


					face_ids.push_back(*point_count + 2);
					face_ids.push_back(*point_count + 8);
					face_ids.push_back(*point_count + 26);
					face_ids.push_back(*point_count + 20);

				}

				// Add face to array
				face_array->InsertNextCell(face_ids.size(), &face_ids[0]);
			}
		}

		// Update offset
		*point_count += 27;
	}
	else {

		// Specify points that make up face
		face_ids.push_back(*point_count + 0);
		face_ids.push_back(*point_count + 2);
		face_ids.push_back(*point_count + 8);
		face_ids.push_back(*point_count + 6);

		// Update offset
		*point_count += 9;

	}

	// Update grid with new point list and add cell
	grid->SetPoints(global_pts);

	if (dimensions_p == 3) {
		grid->InsertNextCell(VTK_POLYHEDRON, 27, &global_ids[0],
			6, face_array->GetPointer());
	}
	else {
		grid->InsertNextCell(VTK_POLYGON, 4, &face_ids[0]);

	}

	return status;

}


/* H5 Multi-Grid Merge Tool for post-processing HDF5 files written by LUMA */
int main(int argc, char* argv[])
{

	// Parse arguments and handle
	std::string case_num("000");
	if (argc == 2)
	{
		std::string arg_str = std::string(argv[1]);

		if (arg_str == "version") {
			std::cout << "H5MultiGridMerge (h5mgm) Version " << H5MGM_VERSION << std::endl;
			return 0;
		}
		else {
			case_num = std::string(argv[1]);
		}
	}
	
	// Print out to screen
	std::cout << "H5MultiGridMerge (h5mgm) Version " << H5MGM_VERSION << std::endl;
	std::cout << "Reconstructing HDF data..." << std::endl;

	// Turn auto error printing off
	H5Eset_auto(H5E_DEFAULT, NULL, NULL);

	// Construct L0 filename
	std::string IN_FILE_NAME("./hdf_R0N0.h5");

	// Set T = 0 time string
	std::string TIME_STRING = "/Time_0";

	// Path for output
	std::string path_str("./vtk_output");

	// Declarations
	herr_t status = 0;
	hid_t output_fid = NULL;
	hid_t input_fid = NULL;
	hid_t input_aid = NULL;
	int dimensions_p, levels, regions, timesteps, out_every;
	double dx;
	std::string VAR;
	int* LatTyp = nullptr;
	double *dummy_d = nullptr;
	int *dummy_i = nullptr;

	// Open L0 input file
	input_fid = H5Fopen(IN_FILE_NAME.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (input_fid == NULL) std::cout << "HDF5 ERROR: Cannot open input file!" << std::endl;

	// Read in key attributes from the file
	input_aid = H5Aopen(input_fid, "Dimensions", H5P_DEFAULT);
	if (input_aid == NULL) std::cout << "HDF5 ERROR: Cannot open attribute!" << std::endl;
	status = H5Aread(input_aid, H5T_NATIVE_INT, &dimensions_p);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot read attribute!" << std::endl;
	status = H5Aclose(input_aid);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot close attribute!" << std::endl;

	input_aid = H5Aopen(input_fid, "NumberOfGrids", H5P_DEFAULT);
	if (input_aid <= 0) std::cout << "HDF5 ERROR: Cannot open attribute!" << std::endl;
	status = H5Aread(input_aid, H5T_NATIVE_INT, &levels);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot read attribute!" << std::endl;
	status = H5Aclose(input_aid);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot close attribute!" << std::endl;

	input_aid = H5Aopen(input_fid, "NumberOfRegions", H5P_DEFAULT);
	if (input_aid <= 0) std::cout << "HDF5 ERROR: Cannot open attribute!" << std::endl;
	status = H5Aread(input_aid, H5T_NATIVE_INT, &regions);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot read attribute!" << std::endl;
	status = H5Aclose(input_aid);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot close attribute!" << std::endl;

	input_aid = H5Aopen(input_fid, "Timesteps", H5P_DEFAULT);
	if (input_aid <= 0) std::cout << "HDF5 ERROR: Cannot open attribute!" << std::endl;
	status = H5Aread(input_aid, H5T_NATIVE_INT, &timesteps);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot read attribute!" << std::endl;
	status = H5Aclose(input_aid);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot close attribute!" << std::endl;

	input_aid = H5Aopen(input_fid, "OutputFrequency", H5P_DEFAULT);
	if (input_aid <= 0) std::cout << "HDF5 ERROR: Cannot open attribute!" << std::endl;
	status = H5Aread(input_aid, H5T_NATIVE_INT, &out_every);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot read attribute!" << std::endl;
	status = H5Aclose(input_aid);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot close attribute!" << std::endl;

	// Store basic data
	int *gridsize = (int*)malloc(3 * sizeof(int));	// Space for grid size
	gridsize[2] = 1;	// Set 3D dimension to 1, will get overwritten if actually 3D
	int num_grids = ((levels - 1) * regions + 1);	// Total number of grids

	// Close file
	status = H5Fclose(input_fid);
	if (status != 0) std::cout << "HDF5 ERROR: Cannot close file!" << std::endl;

	// Create cell node position list for VTK
	std::vector< std::vector<double> > cell_points;
	if (dimensions_p == 3) {
		cell_points.resize(27, std::vector<double>(3));
	}
	else {
		cell_points.resize(9, std::vector<double>(3));
	}

	// Create directory
	std::string command = "mkdir -p " + path_str;
#ifdef _WIN32   // Running on Windows
	CreateDirectoryA((LPCSTR)path_str.c_str(), NULL);
#else   // Running on Unix system
	system(command.c_str());
#endif // _WIN32

	
	// Time loop
	int count = 0; size_t res = 0;
	for (size_t t = 0; t <= (size_t)timesteps; t += out_every) {

		// Create VTK filename
		std::string vtkFilename = path_str + "/luma_" + case_num + "." + std::to_string(t) + ".vtk";

		// Create VTK grid
		vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
			vtkSmartPointer<vtkUnstructuredGrid>::New();

		// Create VTK grid points
		vtkSmartPointer<vtkPoints> points =
			vtkSmartPointer<vtkPoints>::New();

		// Create point ID list
		std::vector<vtkIdType> pointIds;

		// Total point count
		int point_count = 0;

		// Loop 1 to get mesh details
		for (int lev = 0; lev < levels; lev++) {
			for (int reg = 0; reg < regions; reg++) {

				// Debug
				int local_point_count = 0;
				if (t == 0) std::cout << "Interrogating L" << lev << " R" << reg << "...";

				// L0 doesn't have different regions
				if (lev == 0 && reg != 0) continue;

				// Construct input file name
				std::string IN_FILE_NAME("./hdf_R" + std::to_string(reg) + "N" + std::to_string(lev) + ".h5");
				input_fid = H5Fopen(IN_FILE_NAME.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
				if (input_fid <= 0) std::cout << "HDF5 ERROR: Cannot open input file!" << std::endl;

				// Get local grid size
				input_aid = H5Aopen(input_fid, "GridSize", H5P_DEFAULT);
				if (input_aid <= 0) std::cout << "HDF5 ERROR: Cannot open attribute!" << std::endl;
				status = H5Aread(input_aid, H5T_NATIVE_INT, gridsize);
				if (status != 0) std::cout << "HDF5 ERROR: Cannot read attribute!" << std::endl;
				status = H5Aclose(input_aid);
				if (status != 0) std::cout << "HDF5 ERROR: Cannot close attribute!" << std::endl;

				// Get local dx
				input_aid = H5Aopen(input_fid, "Dx", H5P_DEFAULT);
				if (input_aid == NULL) std::cout << "HDF5 ERROR: Cannot open attribute!" << std::endl;
				status = H5Aread(input_aid, H5T_NATIVE_DOUBLE, &dx);
				if (status != 0) std::cout << "HDF5 ERROR: Cannot read attribute!" << std::endl;
				status = H5Aclose(input_aid);
				if (status != 0) std::cout << "HDF5 ERROR: Cannot close attribute!" << std::endl;

				// Allocate space for LatTyp data
				int *Type = (int*)malloc((gridsize[0] * gridsize[1] * gridsize[2]) * sizeof(int));
				if (Type == NULL) {
					std::cout << "Not enough Memory!!!!" << std::endl;
					exit(EXIT_FAILURE);
				}

				// Allocate space for X, Y, Z coordinates
				double *X = (double*)malloc((gridsize[0] * gridsize[1] * gridsize[2]) * sizeof(double));
				if (X == NULL) {
					std::cout << "Not enough Memory!!!!" << std::endl;
					exit(EXIT_FAILURE);
				}
				double *Y = (double*)malloc((gridsize[0] * gridsize[1] * gridsize[2]) * sizeof(double));
				if (Y == NULL) {
					std::cout << "Not enough Memory!!!!" << std::endl;
					exit(EXIT_FAILURE);
				}
				double *Z = (double*)malloc((gridsize[0] * gridsize[1] * gridsize[2]) * sizeof(double));
				if (Z == NULL) {
					std::cout << "Not enough Memory!!!!" << std::endl;
					exit(EXIT_FAILURE);
				}

				// Create input dataspace
				hsize_t dims_input[1];
				dims_input[0] = gridsize[0] * gridsize[1] * gridsize[2];
				hid_t input_sid = H5Screate_simple(1, dims_input, NULL);
				if (input_sid <= 0) std::cout << "HDF5 ERROR: Cannot create input dataspace!" << std::endl;

				// Open, read into buffers and close input datasets
				status = readDataset("/LatTyp", TIME_STRING, input_fid, input_sid, H5T_NATIVE_INT, Type);
				if (status != 0) goto readFailed;
				status = readDataset("/XPos", TIME_STRING, input_fid, input_sid, H5T_NATIVE_DOUBLE, X);
				if (status != 0) goto readFailed;
				status = readDataset("/YPos", TIME_STRING, input_fid, input_sid, H5T_NATIVE_DOUBLE, Y);
				if (status != 0) goto readFailed;
				if (dimensions_p == 3) {
					status = readDataset("/ZPos", TIME_STRING, input_fid, input_sid, H5T_NATIVE_DOUBLE, Z);
					if (status != 0) goto readFailed;
				}

				// Loop over grid sites
				for (int c = 0; c < gridsize[0] * gridsize[1] * gridsize[2]; c++) {

					// Ignore list
					if (
						Type[c] == eRefined ||
						Type[c] == eRefinedInlet ||
						Type[c] == eRefinedSolid ||
						Type[c] == eRefinedSymmetry ||
						Type[c] == eTransitionToCoarser
						) continue;

					local_point_count++;

					if (dimensions_p == 3) {

						// Find positions of all 27 possible points and store
						for (int p = 0; p < 27; p++) {
							cell_points[p][0] = X[c] + e[0][p] * (dx / 2);
							cell_points[p][1] = Y[c] + e[1][p] * (dx / 2);
							cell_points[p][2] = Z[c] + e[2][p] * (dx / 2);
						}

					}
					else {

						// Find positions of all 9 possible points and store
						int vert = 0;
						for (int p = 1; p < 27; p += 3) {
							cell_points[vert][0] = X[c] + e[0][p] * (dx / 2);
							cell_points[vert][1] = Y[c] + e[1][p] * (dx / 2);
							cell_points[vert][2] = 0.0;
							vert++;
						}

					}

					// Build cell from points and add to mesh
					res = addCell(points, unstructuredGrid, pointIds, cell_points, &point_count, dimensions_p);

				}

				// Close input dataspace
				status = H5Sclose(input_sid);
				if (status != 0) std::cout << "HDF5 ERROR: Cannot close input dataspace!" << std::endl;

				// Close input file
				status = H5Fclose(input_fid);
				if (status != 0) std::cout << "HDF5 ERROR: Cannot close input file!" << std::endl;

				// Free memory
				free(Type);
				free(X);
				free(Y);
				free(Z);

				// Debug
				if (t == 0) std::cout << "Valid Point Count  = " << local_point_count << std::endl;

			}
		}

		// Debug
		if (t == 0)	std::cout << "Number of cells retained = " << unstructuredGrid->GetNumberOfCells() << std::endl;

		// Add data
		vtkSmartPointer<vtkIntArray> LatTyp = vtkSmartPointer<vtkIntArray>::New();
		LatTyp->Allocate(unstructuredGrid->GetNumberOfCells());
		LatTyp->SetName("LatTyp");
		addDataToGrid("/LatTyp", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_i, H5T_NATIVE_INT, LatTyp);

		vtkSmartPointer<vtkDoubleArray> RhoArray = vtkSmartPointer<vtkDoubleArray>::New();
		RhoArray->Allocate(unstructuredGrid->GetNumberOfCells());
		RhoArray->SetName("Rho");
		addDataToGrid("/Rho", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, RhoArray);

		vtkSmartPointer<vtkDoubleArray> Rho_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
		Rho_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
		Rho_TimeAv->SetName("Rho_TimeAv");
		addDataToGrid("/Rho_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, Rho_TimeAv);

		vtkSmartPointer<vtkDoubleArray> Ux = vtkSmartPointer<vtkDoubleArray>::New();
		Ux->Allocate(unstructuredGrid->GetNumberOfCells());
		Ux->SetName("Ux");
		addDataToGrid("/Ux", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, Ux);

		vtkSmartPointer<vtkDoubleArray> Uy = vtkSmartPointer<vtkDoubleArray>::New();
		Uy->Allocate(unstructuredGrid->GetNumberOfCells());
		Uy->SetName("Uy");
		addDataToGrid("/Uy", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, Uy);

		vtkSmartPointer<vtkDoubleArray> Ux_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
		Ux_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
		Ux_TimeAv->SetName("Ux_TimeAv");
		addDataToGrid("/Ux_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, Ux_TimeAv);

		vtkSmartPointer<vtkDoubleArray> Uy_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
		Uy_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
		Uy_TimeAv->SetName("Uy_TimeAv");
		addDataToGrid("/Uy_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, Uy_TimeAv);

		vtkSmartPointer<vtkDoubleArray> UxUx_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
		UxUx_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
		UxUx_TimeAv->SetName("UxUx_TimeAv");
		addDataToGrid("/UxUx_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, UxUx_TimeAv);

		vtkSmartPointer<vtkDoubleArray> UxUy_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
		UxUy_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
		UxUy_TimeAv->SetName("UxUy_TimeAv");
		addDataToGrid("/UxUy_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, UxUy_TimeAv);

		vtkSmartPointer<vtkDoubleArray> UyUy_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
		UyUy_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
		UyUy_TimeAv->SetName("UyUy_TimeAv");
		addDataToGrid("/UyUy_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, UyUy_TimeAv);

		if (dimensions_p == 3) {

			vtkSmartPointer<vtkDoubleArray> Uz_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
			Uz_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
			Uz_TimeAv->SetName("Uz_TimeAv");
			addDataToGrid("/Uz_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, Uz_TimeAv);

			vtkSmartPointer<vtkDoubleArray> UxUz_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
			UxUz_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
			UxUz_TimeAv->SetName("UxUz_TimeAv");
			addDataToGrid("/UxUz_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, UxUz_TimeAv);

			vtkSmartPointer<vtkDoubleArray> UyUz_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
			UyUz_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
			UyUz_TimeAv->SetName("UyUz_TimeAv");
			addDataToGrid("/UyUz_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, UyUz_TimeAv);

			vtkSmartPointer<vtkDoubleArray> UzUz_TimeAv = vtkSmartPointer<vtkDoubleArray>::New();
			UzUz_TimeAv->Allocate(unstructuredGrid->GetNumberOfCells());
			UzUz_TimeAv->SetName("UzUz_TimeAv");
			addDataToGrid("/UzUz_TimeAv", TIME_STRING, levels, regions, gridsize, unstructuredGrid, dummy_d, H5T_NATIVE_DOUBLE, UzUz_TimeAv);
		}

		TIME_STRING = "/Time_" + std::to_string(t + out_every);

		// Print progress to screen
		std::cout << std::to_string((int)(((float)(t + out_every) / (float)(timesteps + out_every)) * 100.0f)) << "% complete." << std::endl;

		// Clean mesh to remove duplicates (Paraview not VTK)
		/*vtkSmartPointer<vtkCleanUnstructuredGrid> cleaner =
		vtkSmartPointer<vtkCleanUnstructuredGrid>::New();
		cleaner->SetInputData(unstructuredGrid);
		cleaner->Update();
		cleaner->GetOutput();*/

		// Write grid to file
		vtkSmartPointer<vtkUnstructuredGridWriter> writer =
			vtkSmartPointer<vtkUnstructuredGridWriter>::New();
		writer->SetFileName(vtkFilename.c_str());
		writer->SetFileTypeToBinary();
		writer->SetInputData(unstructuredGrid);
		writer->Write();
	}

	return 0;


readFailed:
	std::cout << "Read failed -- exiting early.";
	return 999;
}

