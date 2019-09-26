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

#ifndef FEMBODY_H
#define FEMBODY_H

// Forward declarations
class IBBody;
#include "FEMNode.h"
#include "FEMElement.h"

/// \brief	Finite element body
///
///			Class for the finite elements body which is contained within the IBBody.
class FEMBody {

	/************** Friends **************/
	friend class IBBody;
	friend class ObjectManager;
	friend class FEMElement;


	/************** Nested classes **************/
	/// \brief	Nested parent element class for each IBM node
	///
	///			Stores which element and where in the element
	///			it each IBM node sits.
	class IBMParentElements {

		// Set friend class
		friend class FEMBody;
		friend class ObjectManager;

		// Constructor and destructor
	public:
		IBMParentElements();
		~IBMParentElements();

		// Members
	private:
		int elementID;
		double zeta;
	};


	/************** Constructors **************/
	FEMBody();
	~FEMBody();

	// Custom constructor for building FEM body from input parameters
	FEMBody(IBBody *iBody, std::vector<double> &start_position, double length,
			double height, double depth, std::vector<double> &angles, int nElements, bool clamped, double density, double E);


	/************** Member Data **************/

	// System values
	IBBody *iBodyPtr;				///< Pointer to owning iBody
	int DOFsPerNode;				///< DOFs per node
	int DOFsPerElement;				///< DOFs per element
	int systemDOFs;					///< DOFs for whole system
	int BC_DOFs;					///< Number of DOFs removed when applying BCs
	int it;							///< Number of iterations for Newton-Raphson solver
	double res;						///< Residual Newton-Raphson solver reached
	double timeav_FEMIterations;	///< Number of iterations for Newton-Raphson solver (time-averaged)
	double timeav_FEMResidual;		///< Residual Newton-Raphson solver reached (time-averaged)

	// Nodes and elements
	std::vector<FEMNode> nodes;				///< Vector of FEM nodes
	std::vector<FEMElement> elements;		///< Vector of FEM elements

	// System matrices
	std::vector<std::vector<double>> M;			///< Mass matrix
	std::vector<std::vector<double>> K;			///< Linear stiffness matrix
	std::vector<double> R;						///< Load vector
	std::vector<double> F;						///< Vector of internal forces
	std::vector<double> U;						///< Vector of displacements
	std::vector<double> U_n;					///< Vector of displacements at start of current time step
	std::vector<double> delU;					///< Vector of incremental displacements
	std::vector<double> Udot;					///< Vector velocities
	std::vector<double> Udot_n;					///< Vector velocities at start of current time step
	std::vector<double> Udotdot;				///< Vector of accelerations
	std::vector<double> Udotdot_n;				///< Vector of accelerations at start of current time step

	// Vector of parent elements for each IBM node
	std::vector<IBMParentElements> IBNodeParents;


	/************** Member Methods **************/

	// Main FEM solver methods
	void dynamicFEM();											// Main outer routine for solving FEM
	void newtonRaphsonIterator();								// Newton-Raphson routine for solve non-linear FEM
	void setNewmark();											// First step in Newmar-Beta time integration
	void finishNewmark();										// Newmark-Beta scheme for getting FEM velocities and accelerations
	void updateFEMValues();										// Update the FEM node data using the new displacements
	void updateIBMarkers();										// Update the IBM markers using new FEM node vales

	// Helper methods
	double checkNRConvergence();								// Check convergence of the Newton-Raphson scheme
	void computeNodeMapping(int nIBMNodes, int nFEMNodes);		// Compute mapping between FEM and IBM nodes
};

#endif
