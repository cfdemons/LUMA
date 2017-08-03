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

#include "../inc/stdafx.h"
#include "../inc/FEMBody.h"
#include "../inc/IBMarker.h"

// *****************************************************************************
/// \brief	Default constructor for FEM body
FEMBody::FEMBody () {

	// Set members to default values
	iBodyPtr = NULL;
	DOFsPerNode = 0;
	DOFsPerElement = 0;
	systemDOFs = 0;
	it = 0;
	res = 0.0;
	timeav_FEMIterations = 0.0;
	timeav_FEMResidual = 0.0;
	BC_DOFs = 0;
}

// *****************************************************************************
/// \brief	Default destructor for FEM body.
FEMBody::~FEMBody () {
}

// *****************************************************************************
///	\brief	Custom constructor for building FEM body from inputs
///
///	\param	iBody			pointer to owning IBBody
///	\param	start_position	start position of FEMBody
///	\param	length			length of FEMBody
///	\param	height			height of FEMBody
///	\param	depth			depth of FEMBody (set to dh in 2D cases)
///	\param	angles			angles of FEMBody
///	\param	nElements		number of elements in FEMBody
///	\param	clamped			boundary condition
///	\param	density			material density
///	\param	E				Young's Modulus
FEMBody::FEMBody (IBBody *iBody, std::vector<double> &start_position, double length,
		double height, double depth, std::vector<double> &angles, int nElements, bool clamped, double density, double E) {

	// Set members to default values
	iBodyPtr = iBody;
	DOFsPerNode = 3;
	DOFsPerElement = 6;
	systemDOFs = (nElements + 1) * DOFsPerNode;
	it = 0;
	res = 0.0;
	timeav_FEMIterations = 0.0;
	timeav_FEMResidual = 0.0;

	// Set number of DOFs to remove in BC
	if (clamped == true)
		BC_DOFs = 3;
	else
		BC_DOFs = 2;

	// Get horizontal and vertical angles
	double body_angle_v = angles[0];
#if (L_DIM == 3)
	double body_angle_h = angles[1];
#else
	double body_angle_h = 0.0;
#endif

	// Compute spacing
	double spacing = length / nElements;							// Physical spacing between markers
	double spacing_h = spacing * cos(body_angle_v * L_PI / 180);	// Local spacing projected onto the horizontal plane

	// Add all FEM nodes
	double x, y, z;
	for (int i = 0; i < nElements + 1; i++) {
		x = start_position[eXDirection] + i * spacing_h * cos(body_angle_h * L_PI / 180.0);
		y = start_position[eYDirection] + i * spacing * sin(body_angle_v * L_PI / 180.0);
		z = start_position[eZDirection] + i * spacing_h * sin(body_angle_h * L_PI / 180.0);
		nodes.emplace_back(i, x, y, z, angles);
	}

// If 2D then set depth to 1 lattice unit
#if (L_DIMS == 2)
	depth = iBodyPtr->_Owner->dh;
#endif

	// Add all FEM elements
	for (int i = 0; i < nElements; i++) {
		elements.emplace_back(i, DOFsPerElement, spacing, height, depth, angles, density, E);
	}

	// Get number of IBM nodes
	int nIBMNodes = floor(length / iBodyPtr->_Owner->dh) + 1;
	int nFEMNodes = nodes.size();

	// Compute IBM-FEM conforming parameters
	computeNodeMapping(nIBMNodes, nFEMNodes);

	// Resize the matrices and set to zero
	M.resize(systemDOFs, std::vector<double>(systemDOFs, 0.0));
	K_L.resize(systemDOFs, std::vector<double>(systemDOFs, 0.0));
	K_NL.resize(systemDOFs, std::vector<double>(systemDOFs, 0.0));
	R.resize(systemDOFs, 0.0);
	F.resize(systemDOFs, 0.0);
	U.resize(systemDOFs, 0.0);
	U_n.resize(systemDOFs, 0.0);
	delU.resize(systemDOFs, 0.0);
	Udot.resize(systemDOFs, 0.0);
	Udot_n.resize(systemDOFs, 0.0);
	Udotdot.resize(systemDOFs, 0.0);
	Udotdot_n.resize(systemDOFs, 0.0);
}


// *****************************************************************************
///	\brief	Main outer routine for FEM solver
void FEMBody::dynamicFEM () {

	// Set values to start of timestep
	U = U_n;
	Udot = Udot_n;
	Udotdot = Udotdot_n;

	// Set elemental values
	for (size_t el = 0; el < elements.size(); el++) {
		elements[el].length = elements[el].length_n;
		elements[el].angles = elements[el].angles_n;
		elements[el].T = elements[el].T_n;
	}

	// Construct load vector as invariant during Newton-Raphson iterations
	constructRVector();

	// While loop parameters
	double TOL = 1e-10;
	double MAXIT = 20;

	// Set while counter to zero
	it = 0;

	// While loop for FEM solver
	do {

		// Solve and iterate over the system
		newtonRaphsonIterator();

		// Check residual
		res = checkNRConvergence();

		// Increment counter
		it++;

	} while (res > TOL && it < MAXIT);

	// Calculate velocities and accelerations
	updateVelocityAndAcceleration();

	// Update IBM markers
	updateIBMarkers();
}

// *****************************************************************************
///	\brief	Newton-Raphson routine for solving non-linear FEM
void FEMBody::newtonRaphsonIterator () {

	// Construct mass matrix
	constructMassMat();

	// Construct stiffness matrix
	constructStiffMat();

	// Construct the internal force vector
	constructFVector();

	// Construct the full non-linear stiffness matrix
	constructNLStiffMat();

	// Declare reduced matrices with BCs applied
	std::vector<std::vector<double>> M_hat(systemDOFs - BC_DOFs, std::vector<double>(systemDOFs - BC_DOFs, 0.0));
	std::vector<std::vector<double>> K_hat(systemDOFs - BC_DOFs, std::vector<double>(systemDOFs - BC_DOFs, 0.0));
	std::vector<double> RmF_hat(systemDOFs - BC_DOFs, 0.0);
	std::vector<double> delU_hat(systemDOFs - BC_DOFs, 0.0);

	// Apply the boundary conditions
	bcFEM(M_hat, K_hat, RmF_hat);

	// Apply Newmark scheme (using Newmark coefficients)
	setNewmark(M_hat, K_hat, RmF_hat);

	// Solve linear system using LAPACK library
	delU_hat = GridUtils::solveLinearSystem(K_hat, RmF_hat);

	// Assign displacement to delU
	for (int i = 0; i < systemDOFs - BC_DOFs; i++) {
		delU[i+BC_DOFs] = delU_hat[i];
	}

	// Add deltaU to U
	for (int i = 0; i < systemDOFs; i++) {
		U[i] = U[i] + delU[i];
	}

	// Update FEM positions
	updateFEMNodes();
}

// *****************************************************************************
///	\brief	Check convergence of the Newton-Raphson scheme
///
/// \return	residual
double FEMBody::checkNRConvergence () {

	// Get original length of body
	double L = elements[0].length0 * elements.size();

	// Normalise the displacement
	std::vector<double> delUStar(delU.size(), 0.0);
	for (size_t i = 0; i < delUStar.size(); i++)
		delUStar[i] = delU[i] / L;

	// Get the norm but divide by length so it is resolution independent and return
	return sqrt(GridUtils::dotprod(delUStar, delUStar) / delUStar.size());
}


// *****************************************************************************
///	\brief	Newmark-Beta scheme for getting FEM velocities and accelerations
void FEMBody::updateVelocityAndAcceleration () {

	// Get velocity and acceleration using Newmark-Beta scheme
	double Dt = iBodyPtr->_Owner->dt;
	double UdotdotStar;

	// Newmark coefficients
	double a6 = 1.0 / (L_NB_ALPHA * SQ(Dt));
	double a7 = -1.0 / (L_NB_ALPHA * Dt);
	double a8 = -(1.0 / (2.0 * L_NB_ALPHA) - 1.0);
	double a9 = Dt * (1.0 - L_NB_DELTA);
	double a10 = L_NB_DELTA * Dt;

	// Update velocities and accelerations
	for (int i = 0; i < systemDOFs; i++) {
		UdotdotStar = a6 * (U[i] - U_n[i]) + a7 * Udot[i] + a8 * Udotdot[i];
		Udot[i] = Udot[i] + a9 * Udotdot[i] + a10 * UdotdotStar;
		Udotdot[i] = UdotdotStar;
	}
}


// *****************************************************************************
///	\brief	Construct load vector
void FEMBody::constructRVector() {

	// Initialise arrays and matrices for calculating load vector
	std::vector<std::vector<double>> T(L_DIMS, std::vector<double>(L_DIMS, 0.0));
	std::vector<double> Rlocal(DOFsPerElement, 0.0);
	std::vector<double> RGlobal(DOFsPerElement, 0.0);
	std::vector<double> F(L_DIMS, 0.0);

	// Required parameters
	double forceScale, length, a, b;
	int IBnode;

	// Set R vector to zero
	fill(R.begin(), R.end(), 0.0);

	// Get force scaling parameter
	forceScale = iBodyPtr->_Owner->dm / SQ(iBodyPtr->_Owner->dt);

	// Loop through all FEM elements
	for (size_t el = 0; el < elements.size(); el++) {

		// Get element parameters
		length = elements[el].length;

		// Now loop through all child IB nodes this element has
		for (size_t node = 0; node < elements[el].IBChildNodes.size(); node++) {

			// Get IB node and integration ranges
			IBnode = elements[el].IBChildNodes[node].nodeID;
			a = elements[el].IBChildNodes[node].zeta1;
			b = elements[el].IBChildNodes[node].zeta2;

			// Get subset of transpose matrix
			T = {{elements[el].T[0][0], elements[el].T[0][1]},
				 {elements[el].T[1][0], elements[el].T[1][1]}};

			// Convert force to local coordinates
			F = GridUtils::matrix_multiply(T, GridUtils::vecmultiply(iBodyPtr->markers[IBnode].epsilon * 1.0 * forceScale, iBodyPtr->markers[IBnode].force_xyz));

			// Get the nodal values by integrating over range of IB point
			Rlocal[0] = F[0] * 0.5 * length * (0.5 * b - 0.5 * a + 0.25 * SQ(a) - 0.25 * SQ(b));
			Rlocal[1] = F[1] * 0.5 * length * (0.5 * b - 0.5 * a - SQ(a) * SQ(a) / 16.0 + SQ(b) * SQ(b) / 16.0 + 3.0 * SQ(a) / 8.0 - 3.0 * SQ(b) / 8.0);
			Rlocal[2] = F[1] * 0.5 * length * (length * (-SQ(a) * SQ(a) + SQ(b) * SQ(b)) / 32.0 - length * (-TH(a) + TH(b)) / 24.0 - length * (-SQ(a) + SQ(b)) / 16.0 + length * (b - a) / 8.0);
			Rlocal[3] = F[0] * 0.5 * length * (-0.25 * SQ(a) + 0.25 * SQ(b) + 0.5 * b - 0.5 * a);
			Rlocal[4] = F[1] * 0.5 * length * (0.5 * b - 0.5 * a + SQ(a) * SQ(a) / 16.0 - SQ(b) * SQ(b) / 16.0 - 3.0 * SQ(a) / 8.0 + 3.0 * SQ(b) / 8.0);
			Rlocal[5] = F[1] * 0.5 * length * (length * (-SQ(a) * SQ(a) + SQ(b) * SQ(b)) / 32.0 + length * (-TH(a) + TH(b)) / 24.0 - length * (-SQ(a) + SQ(b)) / 16.0 - length * (b - a) / 8.0);

			// Get element internal forces
			RGlobal = GridUtils::matrix_multiply(GridUtils::matrix_transpose(elements[el].T), Rlocal);

			// Add to global vector
			GridUtils::assembleGlobalVec(el, DOFsPerNode, RGlobal, R);
		}
	}
}

// *****************************************************************************
///	\brief	Construct mass matrix
void FEMBody::constructMassMat () {

	// Initialise arrays and matrices for calculating mass matrix
	std::vector<std::vector<double>> Mglobal(DOFsPerElement, std::vector<double>(DOFsPerElement, 0.0));
	std::vector<std::vector<double>> Mlocal(DOFsPerElement, std::vector<double>(DOFsPerElement, 0.0));

	// Reset mass matrix to zero
	for (size_t i = 0; i < M.size(); i++) {
		fill(M[i].begin(), M[i].end(), 0.0);
	}

	// Coefficients
	double A, rho, L0, C1;

	// Loop through each element and create stiffness matrix
	for (size_t el = 0; el < elements.size(); el++) {

		// Coefficients
		A = elements[el].area;
		rho = elements[el].density;
		L0 = elements[el].length0;
		C1 = rho * A * L0 / 420.0;

		// Construct part 1 of mass matrix (axial and transverse)
		Mlocal[0][0] = C1 * 140.0;
		Mlocal[0][3] = C1 * 70.0;
		Mlocal[1][1] = C1 * 156.0;
		Mlocal[1][2] = C1 * 22.0 * L0;
		Mlocal[1][4] = C1 * 54;
		Mlocal[1][5] = C1 * (-13.0 * L0);
		Mlocal[2][2] = C1 * 4.0 * SQ(L0);
		Mlocal[2][4] = C1 * 13.0 * L0;
		Mlocal[2][5] = C1 * (-3.0 * SQ(L0));
		Mlocal[3][3] = C1 * 140.0;
		Mlocal[4][4] = C1 * 156.0;
		Mlocal[4][5] = C1 * (-22.0 * L0);
		Mlocal[5][5] = C1 * 4.0 * SQ(L0);

		// Copy to the lower half (symmetrical matrix)
		for (int i = 1; i < DOFsPerElement; i++) {
			for (int j = 0; j < i; j++) {
				Mlocal[i][j] = Mlocal[j][i];
			}
		}

		// Multiply by transformation matrices to get global matrix for single element
		Mglobal = GridUtils::matrix_multiply(GridUtils::matrix_multiply(GridUtils::matrix_transpose(elements[el].T), Mlocal), elements[el].T);

		// Add to global matrix
		GridUtils::assembleGlobalMat(el, DOFsPerNode, Mglobal, M);
	}
}

// *****************************************************************************
///	\brief	Construct linear stiffness matrix
void FEMBody::constructStiffMat () {

	// Initialise arrays and matrices for calculating linear stiffness matrix
	std::vector<std::vector<double>> Kglobal(DOFsPerElement, std::vector<double>(DOFsPerElement, 0.0));
	std::vector<std::vector<double>> Klocal(DOFsPerElement, std::vector<double>(DOFsPerElement, 0.0));


	// Reset linear stiffness matrix to zero
	for (size_t i = 0; i < K_L.size(); i++) {
		fill(K_L[i].begin(), K_L[i].end(), 0.0);
	}

	// Coefficients
	double A, E, I, L0;

	// Loop through each element and create stiffness matrix
	for (size_t el = 0; el < elements.size(); el++) {

		// Coefficients
		A = elements[el].area;
		E = elements[el].E;
		I = elements[el].I;
		L0 = elements[el].length0;

		// Construct upper half of local stiffness matrix for single element
		Klocal[0][0] = E * A / L0;
		Klocal[0][3] = -E * A / L0;
		Klocal[1][1] = 12.0 * E * I / TH(L0);
		Klocal[1][2] = 6.0 * E * I / SQ(L0);
		Klocal[1][4] = -12.0 * E * I / TH(L0);
		Klocal[1][5] = 6.0 * E * I / SQ(L0);
		Klocal[2][2] = 4.0 * E * I / L0;
		Klocal[2][4] = -6.0 * E * I / SQ(L0);
		Klocal[2][5] = 2.0 * E * I / L0;
		Klocal[3][3] = E * A / L0;
		Klocal[4][4] = 12.0 * E * I / TH(L0);
		Klocal[4][5] = -6.0 * E * I / SQ(L0);
		Klocal[5][5] = 4.0 * E * I / L0;

		// Copy to the lower half (symmetrical matrix)
		for (int i = 1; i < DOFsPerElement; i++) {
			for (int j = 0; j < i; j++) {
				Klocal[i][j] = Klocal[j][i];
			}
		}

		// Multiply by transformation matrices to get global matrix for single element
		Kglobal = GridUtils::matrix_multiply(GridUtils::matrix_multiply(GridUtils::matrix_transpose(elements[el].T), Klocal), elements[el].T);

		// Add to global matrix
		GridUtils::assembleGlobalMat(el, DOFsPerNode, Kglobal, K_L);
	}
}

// *****************************************************************************
///	\brief	Construct internal force vector
void FEMBody::constructFVector () {

	// Element internal forces in global coordinates
	std::vector<double> FGlobal(DOFsPerElement, 0.0);

	// Declare values
	double E, I, A, L0, L;
	double angleElement, angleNode1, angleNode2;
	double u, theta1, theta2;
	double F0, M1, M2;

	// Reset force vector
	fill(F.begin(), F.end(), 0.0);

	// Loop through each element
	for (size_t el = 0; el < elements.size(); el++) {

		// Set various values for element
		E = elements[el].E;
		I = elements[el].I;
		A = elements[el].area;
		L0 = elements[el].length0;
		L = elements[el].length;
		angleElement = elements[el].angles;
		angleNode1 = nodes[el].angles;
		angleNode2 = nodes[el+1].angles;

		// Calculate the local displacements
		u = (SQ(L) - SQ(L0)) / (L + L0);
		theta1 = atan2(cos(angleElement) * sin(angleNode1) - sin(angleElement) * cos(angleNode1), cos(angleElement) * cos(angleNode1) + sin(angleElement) * sin(angleNode1));
		theta2 = atan2(cos(angleElement) * sin(angleNode2) - sin(angleElement) * cos(angleNode2), cos(angleElement) * cos(angleNode2) + sin(angleElement) * sin(angleNode2));

		// Calculate internal forces for beam
		F0 = (E * A / L0) * u;
		M1 = (2 * E * I / L0) * (2.0 * theta1 + theta2);
		M2 = (2 * E * I / L0) * (theta1 + 2.0 * theta2);

		// Set the internal local nodal forces
		elements[el].F[0] = -F0;
		elements[el].F[1] = (1.0 / L0) * (M1 + M2);
		elements[el].F[2] = M1;
		elements[el].F[3] = F0;
		elements[el].F[4] = -(1.0 / L0) * (M1 + M2);
		elements[el].F[5] = M2;

		// Get element internal forces
		FGlobal = GridUtils::matrix_multiply(GridUtils::matrix_transpose(elements[el].T), elements[el].F);

		// Add to global vector
		GridUtils::assembleGlobalVec(el, DOFsPerNode, FGlobal, F);
	}
}

// *****************************************************************************
///	\brief	Construct nonlinear stiffness matrix
void FEMBody::constructNLStiffMat () {

	// Initialise arrays and matrices for calculating non-linear stiffness matrix
	std::vector<std::vector<double>> Kglobal(DOFsPerElement, std::vector<double>(DOFsPerElement, 0.0));
	std::vector<std::vector<double>> Klocal(DOFsPerElement, std::vector<double>(DOFsPerElement, 0.0));

	// Reset non-linear stiffness matrix to zero
	for (size_t i = 0; i < K_NL.size(); i++) {
		fill(K_NL[i].begin(), K_NL[i].end(), 0.0);
	}

	// Coefficients
	double L0, F0, V0;

	// Loop through each element and create stiffness matrix
	for (size_t el = 0; el < elements.size(); el++) {

		// Calculate length of element
		L0 = elements[el].length0;

		// Internal forces
		F0 = -elements[el].F[0];
		V0 = elements[el].F[4];

		// Construct upper half of local stiffness matrix for single element
		Klocal[0][1] = -V0 / L0;
		Klocal[0][4] = -(F0 / L0) + V0 / L0;
		Klocal[1][0] = -V0 / L0;
		Klocal[1][1] = F0 / L0;
		Klocal[1][3] = V0 / L0;
		Klocal[3][1] = V0 / L0;
		Klocal[3][4] = -V0 / L0;
		Klocal[4][0] = V0 / L0;
		Klocal[4][1] = -F0 / L0;
		Klocal[4][3] = -V0 / L0;
		Klocal[4][4] = F0 / L0;

		// Multiply by transformation matrices to get global matrix for single element
		Kglobal = GridUtils::matrix_multiply(GridUtils::matrix_multiply(GridUtils::matrix_transpose(elements[el].T), Klocal), elements[el].T);

		// Add to global matrix
		GridUtils::assembleGlobalMat(el, DOFsPerNode, Kglobal, K_NL);
	}
}


// *****************************************************************************
///	\brief	Apply BCs by removing elements in global matrices
///
///	\param	M_hat		mass matrix with BCs about to be applied
///	\param	K_hat		stiffness matrix with BCs about to be applied
///	\param	RmF_hat		balanced load vector with BCs about to be applied
void FEMBody::bcFEM (std::vector<std::vector<double>> &M_hat, std::vector<std::vector<double>> &K_hat, std::vector<double> &RmF_hat) {

	// Loop through reduced size
	for (int i = 0; i < systemDOFs - BC_DOFs; i++) {

		// Unbalanced load vector
		RmF_hat[i] = R[i+BC_DOFs] - F[i+BC_DOFs];

		// Now loop through mass and stiffness matrices
		for (int j = 0; j < systemDOFs - BC_DOFs; j++) {
			M_hat[i][j] = M[i+BC_DOFs][j+BC_DOFs];
			K_hat[i][j] = K_L[i+BC_DOFs][j+BC_DOFs] + K_NL[i+BC_DOFs][j+BC_DOFs];
		}
	}
}

// *****************************************************************************
///	\brief	First step in Newmark-Beta time integration
///
///	\param	M_hat		mass matrix with BCs applied
///	\param	K_hat		stiffness matrix with BCs applied
///	\param	RmF_hat		balanced load vector with BCs applied
void FEMBody::setNewmark (std::vector<std::vector<double>> &M_hat, std::vector<std::vector<double>> &K_hat, std::vector<double> &RmF_hat) {

	// Newmark-beta method for time integration
	double Dt = iBodyPtr->_Owner->dt;
	double a0, a2, a3;
	a0 = 1.0 / (L_NB_ALPHA * SQ(Dt));
	a2 = 1.0 / (L_NB_ALPHA * Dt);
	a3 = 1.0 / (2.0 * L_NB_ALPHA) - 1.0;

	// Calculate effective load vector
	std::vector<double> Meff_hat(systemDOFs - BC_DOFs, 0.0);
	for (int i = 0; i < systemDOFs - BC_DOFs; i++) {
		Meff_hat[i] = a0 * (U_n[i+BC_DOFs] - U[i+BC_DOFs]) + a2 * Udot[i+BC_DOFs] + a3 * Udotdot[i+BC_DOFs];
	}

	// Multiply with mass matrix to get inertia forces
	std::vector<double> MF_hat(systemDOFs - BC_DOFs, 0.0);
	MF_hat = GridUtils::matrix_multiply(M_hat, Meff_hat);

	// Calculate effective load vector and stiffness matrix
	for (int i = 0; i < systemDOFs - BC_DOFs; i++) {

		// Effective load
		RmF_hat[i] += MF_hat[i];

		// Effective stiffness
		for (int j = 0; j < systemDOFs - BC_DOFs; j++) {
			K_hat[i][j] += a0 * M_hat[i][j];
		}
	}
}

// *****************************************************************************
///	\brief	Update the new FEM node data using the displacements
void FEMBody::updateFEMNodes () {

	// Set the new positions in the particle_struct
	for (size_t n = 0; n < nodes.size(); n++) {

		// Positions
		for (int d = 0; d < L_DIMS; d++)
			nodes[n].position[d] = nodes[n].position0[d] + U[n*DOFsPerNode+d];

		// Set angle
		nodes[n].angles = nodes[n].angles0 + U[n*DOFsPerNode+L_DIMS];
	}

	// Set the new angles and lengths of the elements
	std::vector<double> elVector;
	for (size_t el = 0; el < elements.size(); el++) {

		// Set the new angles and lengths of the elements
		elVector = GridUtils::subtract(nodes[el+1].position, nodes[el].position);
		elements[el].angles = atan2(elVector[eYDirection], elVector[eXDirection]);
		elements[el].length = GridUtils::vecnorm(elVector);

		// Set new transformation matrix
		elements[el].T[0][0] = elements[el].T[1][1] =  elements[el].T[3][3] = elements[el].T[4][4] = cos(elements[el].angles);
		elements[el].T[0][1] = elements[el].T[3][4] = sin(elements[el].angles);
		elements[el].T[1][0] = elements[el].T[4][3] = -sin(elements[el].angles);
		elements[el].T[2][2] = elements[el].T[5][5] =  1.0;
	}
}

// *****************************************************************************
///	\brief	Update the new IBMaker data using the FEM data
void FEMBody::updateIBMarkers() {

	// Parameters
	int el;
	double zeta, length;
	std::vector<double> ULocal, UnodeLocal, UnodeGlobal;
	std::vector<double> UDotLocal, UDotNodeLocal, UDotNodeGlobal;
	std::vector<double> UGlobal(DOFsPerElement, 0.0);
	std::vector<double> UDotGlobal(DOFsPerElement, 0.0);
	std::vector<std::vector<double>> T(L_DIMS, std::vector<double>(L_DIMS, 0.0));

	// Loop through all IBM nodes
	for (size_t node = 0; node < IBNodeParents.size(); node++) {

		// Get element this IB node exists within
		el = IBNodeParents[node].elementID;
		zeta = IBNodeParents[node].zeta;
		length = elements[el].length;

		// Get the element values
		GridUtils::disassembleGlobalVec(el, DOFsPerNode, U, UGlobal);
		GridUtils::disassembleGlobalVec(el, DOFsPerNode, Udot, UDotGlobal);

		// Get element values in local coordinates
		ULocal = GridUtils::matrix_multiply(elements[el].T, UGlobal);
		UDotLocal = GridUtils::matrix_multiply(elements[el].T, UDotGlobal);

		// Get the local displacement of the IB node
		UnodeLocal = shapeFunctions(ULocal, zeta, length);
		UDotNodeLocal = shapeFunctions(UDotLocal, zeta, length);

		// Get subset of transformation matrix
		T = {{elements[el].T[0][0], elements[el].T[0][1]},
			 {elements[el].T[1][0], elements[el].T[1][1]}};

		// Get the IB node displacements in global coordinates
		UnodeGlobal = GridUtils::matrix_multiply(GridUtils::matrix_transpose(T), UnodeLocal);
		UDotNodeGlobal = GridUtils::vecmultiply(iBodyPtr->_Owner->dt / iBodyPtr->_Owner->dh, GridUtils::matrix_multiply(GridUtils::matrix_transpose(T), UDotNodeLocal));

		// Set the IBM node
		for (int d = 0; d < L_DIMS; d++) {
			iBodyPtr->markers[node].position[d] = iBodyPtr->markers[node].position0[d] + UnodeGlobal[d];
			iBodyPtr->markers[node].markerVel_km1[d] = iBodyPtr->markers[node].markerVel[d];
			iBodyPtr->markers[node].markerVel[d] = L_RELAX * UDotNodeGlobal[d] + (1.0 - L_RELAX) * iBodyPtr->markers[node].markerVel_km1[d];
		}
	}
}

// *****************************************************************************
///	\brief	Get shape functions given natural value
///
///	\param DOFVec	vector containing the values in natural coordinates
///	\param zeta		the natural coordinate
/// \param length	length of element
/// \return value at the coordinate
std::vector<double> FEMBody::shapeFunctions (std::vector<double> &DOFVec, double zeta, double length) {

	// Results vector
	std::vector<double> resVec(L_DIMS, 0.0);

	// Use shape functions to calculate values
	double N0 = 1.0 - (zeta + 1.0) / 2.0;
	double N1 = 1.0 - 3.0 * SQ((zeta + 1.0) / 2.0) + 2.0 * TH((zeta + 1.0) / 2.0);
	double N2 = ((zeta + 1.0) / 2.0 - 2.0 * SQ((zeta + 1.0) / 2.0) + TH((zeta + 1.0) / 2.0)) * length;
	double N3 = (zeta + 1.0) / 2.0;
	double N4 = 3.0 * SQ((zeta + 1.0) / 2.0) - 2.0 * TH((zeta + 1.0) / 2.0);
	double N5 = (-SQ((zeta + 1.0) / 2.0) + TH((zeta + 1.0) / 2.0)) * length;

	// Calculate values using shape functions
	resVec[eXDirection] = DOFVec[0] * N0 + DOFVec[3] * N3;
	resVec[eYDirection] = DOFVec[1] * N1 + DOFVec[2] * N2 + DOFVec[4] * N4 + DOFVec[5] * N5;

	// Return
	return resVec;
}


// *****************************************************************************
///	\brief	Compute the mapping parameters for FEM-IBM nodes
///
///	\param nIBMNodes	number of IBM nodes in body
///	\param nFEMNodes	number of FEM nodes in body
void FEMBody::computeNodeMapping (int nIBMNodes, int nFEMNodes) {

	// Set the IBM parent elements
	IBNodeParents.resize(nIBMNodes);

	// Set first node first
	IBNodeParents[0].elementID = 0;
	IBNodeParents[0].zeta = -1.0;

	// Now loop through and set values
	for (int i = 1; i < nIBMNodes; i++) {

		// Check if remainder is zero
		if (fmod(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)), 1.0) == 0.0) {
			IBNodeParents[i].elementID = floor(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0))) - 1;
			IBNodeParents[i].zeta = 1.0;
		}
		else {
			IBNodeParents[i].elementID = floor(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)));
			IBNodeParents[i].zeta = fmod(i * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)), 1.0) * 2.0 - 1.0;
		}
	}

	// Loop through elements
	double node1, node2;
	for (size_t el = 0; el < elements.size(); el++) {

		// Loop through all nodes and scale range to local coordinates for element
		for (int node = 0; node < nIBMNodes; node++) {
			node1 = -1 + 2.0 * ((node - 0.5) * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)) - el) / 1.0;
			node2 = -1 + 2.0 * ((node + 0.5) * ((nFEMNodes - 1.0) / (nIBMNodes - 1.0)) - el) / 1.0;

			// Check if any points lie within element coordinate range
			if ((node1 > -1.0 && node1 < 1.0) || (node2 > -1.0 && node2 < 1.0)) {

				// Sort out end nodes where one point lie outside element
				if (node1 < -1.0)
					node1 = -1.0;
				if (node2 > 1.0)
					node2 = 1.0;

				// Call constructor for chile IB point
				elements[el].IBChildNodes.emplace_back(node, node1, node2);
			}
		}
	}
}

// *****************************************************************************
///	\brief	Default constructor for parent element class
FEMBody::IBMParentElements::IBMParentElements() {

	// Set default values
	elementID = 0;
	zeta = 0.0;
}

// *****************************************************************************
///	\brief	Default destructor for parent element class
FEMBody::IBMParentElements::~IBMParentElements() {
}

