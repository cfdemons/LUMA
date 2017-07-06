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

/// \brief	Finite element body.
class FEMBody {

	/************** Friends **************/
	friend class IBBody;
	friend class ObjectManager;


	/************** Constructors **************/
	FEMBody();
	~FEMBody();

	// Custom constructor for building FEM body from input parameters
	FEMBody(IBBody *iBody, std::vector<double> &start_position, double length,
			double height, double depth, std::vector<double> &angles, int nElements, bool clamped, double density, double E);


	/************** Member Data **************/

	// System values
	IBBody *iBodyPtr;			///< Pointer to owning iBody
	int DOFsPerNode;			///< DOFs per node
	int DOFsPerElement;			///< DOFs per element
	int systemDOFs;				///< DOFs for whole system
	int BC_DOFs;				///< Number of DOFs removed when applying BCs
	int NRIterations;			///< Number of iterations for Newton-Raphson solver
	double NRResidual;			///< Residual Newton-Raphson solver reached

	// Nodes and elements
	std::vector<FEMNode> node;				///< Vector of FEM nodes
	std::vector<FEMElement> element;		///< Vector of FEM elements

	// System matrices
	std::vector<std::vector<double>> M;			///< Mass matrix
	std::vector<std::vector<double>> K_L;		///< Linear stiffness matrix
	std::vector<std::vector<double>> K_NL;		///< Non-linear stiffness matrix
	std::vector<double> R;						///< Load vector
	std::vector<double> F;						///< Vector of internal forces
	std::vector<double> U;						///< Vector of displacements
	std::vector<double> delU;					///< Vector of incremental displacements
	std::vector<double> Udot;					///< Vector velocities
	std::vector<double> Udotdot;				///< Vector of accelerations


	/************** Member Methods **************/

};

#endif
