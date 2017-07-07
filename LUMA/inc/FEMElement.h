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
#ifndef FEMELEMENT_H
#define FEMELEMENT_H

// Forward declarations

/// \brief	Finite element body.
class FEMElement {

	/************** Friends **************/
	friend class FEMBody;


	/************** Nested classes **************/
	/// \brief	Nested child IBM node class for FEM elements.
	///
	///			Stores which nodes and where along itself each IBM
	///			node sits.
	class FEMChildNodes {

		// Set friend class
		friend class FEMBody;

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
	FEMElement(int i, int DOFs, double spacing, double height, double depth, std::vector<double> &angles, double density, double E);


	/************** Member Data **************/
private:

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

	// Vector of child IBM nodes which exist along this element
	std::vector<FEMChildNodes> IBChildNodes;		///< Vector of child IBM nodes which exist along this element


	/************** Member Methods **************/

};

#endif
