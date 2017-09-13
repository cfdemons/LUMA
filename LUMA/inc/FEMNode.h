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
#ifndef FEMNODE_H
#define FEMNODE_H

// Forward declarations

/// \brief	Finite element node class
///
///			Class for the finite element nodes which are contained within the FEMBody.
class FEMNode {

	/************** Friends **************/
	friend class FEMBody;
	friend class ObjectManager;

	/************** Constructors **************/
public:
	FEMNode();
	~FEMNode();

	// Custom constructor for building FEM node from input parameters
	FEMNode(int idx, double x, double y, double z, std::vector<double> &angles);


	/************** Member Data **************/

	// Node values
private:
	int ID;								///< Node ID in FEM body
	std::vector<double> position0;		///< Initial position of FEM node
	std::vector<double> position;		///< Current position of FEM node
	double angles0;						///< Initial angles of FEM node
	double angles;						///< Current angles of FEM node


	/************** Member Methods **************/

};

#endif
