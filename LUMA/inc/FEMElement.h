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


	/************** Constructors **************/
public:
	FEMElement();
	~FEMElement();

	// Custom constructor for building FEM body from input parameters
	FEMElement(int i, int DOFs, double spacing, double height, double depth, std::vector<double> &angles, double density, double E);


	/************** Member Data **************/
private:

	// Element ID
	int ID;								///< Element ID in FEM body

	// Geometry properties
	double length0;						///< Initial length of element
	double length;						///< Current length of element
	std::vector<double> angles;			///< Current orientation of element
	double area;						///< Cross-sectional area of element
	std::vector<double> I;				///< Second moment areas

	// Structural properties
	double E;							///< Youngs modulus
	double density;						///< Material density

	// Internal forces
	std::vector<double> F;				///< Vector of internal forces


	/************** Member Methods **************/

};

#endif
