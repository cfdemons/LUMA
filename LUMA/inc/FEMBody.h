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
#ifndef FEMBODY_H
#define FEMBODY_H

// Forward declarations
#include "IBBody.h"
#include "FEMNode.h"
#include "FEMElement.h"

class IBBody;

/// \brief	Finite element body
///
///			Class for the finite elements body which is contained within the IBBody.
class FEMBody {

	/************** Friends **************/
	friend class IBBody;
	friend class ObjectManager;


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
	std::vector<std::vector<double>> K_L;		///< Linear stiffness matrix
	std::vector<std::vector<double>> K_NL;		///< Non-linear stiffness matrix
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
	void updateVelocityAndAcceleration();						// Newmark-Beta scheme for getting FEM velocities and accelerations
	void constructRVector();									// Construct load vector
	void constructMassMat();									// Construct mass matrix
	void constructStiffMat();									// Construct linear stiffness matrix
	void constructFVector();									// Construct internal force vector
	void constructNLStiffMat();									// Construct non-linear stiffness matrix
	void updateFEMNodes();										// Update the FEM node data using the new displacements
	void updateIBMarkers(double relax = L_RELAX);				// Update the IBM markers using new FEM node vales
	std::vector<double> shapeFunctions(std::vector<double> &vec, double zeta, double length);											// Sum the shape functions to get displacement/velocity
	void bcFEM(std::vector<std::vector<double>> &M_hat, std::vector<std::vector<double>> &K_hat, std::vector<double> &RmF_hat);			// Apply BCs by removing elements in global matrices
	void setNewmark(std::vector<std::vector<double>> &M_hat, std::vector<std::vector<double>> &K_hat, std::vector<double> &RmF_hat);	// First step in Newmar-Beta time integration

	// Helper methods
	double checkNRConvergence();								// Check convergence of the Newton-Raphson scheme
	void computeNodeMapping(int nIBMNodes, int nFEMNodes);		// Compute mapping between FEM and IBM nodes

};

#endif
