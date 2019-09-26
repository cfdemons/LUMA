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

#ifndef FEMELEMENT_H
#define FEMELEMENT_H

// Forward declarations
class FEMBody;
class IBBody;
class IBMarker;

/// \brief	Finite element class
///
///			Class for the finite elements which are contained within the FEMBody.
class FEMElement {

	/************** Friends **************/
	friend class FEMBody;
	friend class ObjectManager;

	/************** Nested classes **************/
	/// \brief	Nested child IBM node class for FEM elements
	///
	///			Stores which nodes and where along itself each IBM
	///			node sits.
	class FEMChildNodes {

		// Set friend class
		friend class FEMElement;
		friend class ObjectManager;

		// Constructor and destructor
	public:
		FEMChildNodes();
		~FEMChildNodes();

		// Custom constructor for initialising vector element
		FEMChildNodes(int node, double zetaA, double zetaB);

		// Members
	private:
		int nodeID;
		double zeta1;
		double zeta2;
	};

	/************** Constructors **************/
public:
	FEMElement();
	~FEMElement();

	// Custom constructor for building FEM body from input parameters
	FEMElement(FEMBody *fPtr, int i, int DOFs, double spacing, double height, double depth, std::vector<double> &angles, double density, double E);


	/************** Member Data **************/
private:

	// Pointer to fBody
	FEMBody *fPtr;

	// Element ID
	int ID;											///< Element ID in FEM body

	// Geometry properties
	double length0;									///< Initial length of element
	double length;									///< Current length of element
	double angles;									///< Current orientation of element
	double area;									///< Cross-sectional area of element
	double I;										///< Second moment areas

	// Structural properties
	double E;										///< Youngs modulus
	double density;									///< Material density

	// Transformation matrix
	std::vector<std::vector<double>> T;				///< Local transformation matrix

	// Internal forces
	std::vector<double> F;							///< Vector of internal forces

	// Global DOFs
	std::vector<int> DOFs;							///< Global DOFs for this element

	// Vector of child IBM nodes which exist along this element
	std::vector<FEMChildNodes> IBChildNodes;		///< Vector of child IBM nodes which exist along this element


	/************** Member Methods **************/

	// Elemental shape functions & matrices
	std::vector<double> shapeFuns(const std::vector<double> &vec, double zeta);		// Elemental shape functions
	void loadVector();																// Construct elemental load vector
	void massMatrix();																// Construct elemental mass matrix
	void stiffMatrix();																// Construct elemental stiffness matrix
	void forceVector();																// Construct elemental internal force vector

	// Assembly methods
	void assembleGlobalMat(const std::vector<double> &localVec, std::vector<double> &globalVec);							// Assemble into global vector
	void assembleGlobalMat(const std::vector<std::vector<double>> &localMat, std::vector<std::vector<double>> &globalMat);	// Assemble into global matrix
	std::vector<double> disassembleGlobalMat(const std::vector<double> &globalVec);											// Disassemble global vector

};

#endif
