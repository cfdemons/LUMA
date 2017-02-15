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
#ifndef BFLBODY_H
#define BFLBODY_H

#include "stdafx.h"
#include "Body.h"
#include "BFLMarker.h"
class PCpts;

/// \brief	BFL body.
///
///			A BFL body is made up of a collection of BFLMarkers.
class BFLBody :
	public Body<BFLMarker>
{
	friend class ObjectManager;
	friend class GridObj;

public:
	// Default constructor and destructor
	BFLBody(void);
	~BFLBody(void);
	// Custom constructor which takes pointer to point cloud data and a pointer to the grid owner for the labelling
	BFLBody(GridObj *g, int bodyID, PCpts *_PCpts);
	BFLBody(GridObj* g, int bodyID, int lev, int reg, std::vector<double> &centre_point, double radius);
	BFLBody(GridObj* g, int bodyID, int lev, int reg, std::vector<double> &centre_point,
			std::vector<double> &dimensions, std::vector<double> &angles);
	BFLBody(GridObj* g, int bodyID, int lev, int reg, std::vector<double> &start_position,
		double length, std::vector<double> &angles);

protected:

	/************** Member Data **************/

	/// \brief	Distance between adjacent lattice site and the surface of the body.
	///
	///			There are two stores of values. Store 1 is the distance on one 
	///			side of the wall and store 2 the distance on the other side. 
	///			One store is appended to the other in this structure.
	std::vector< std::vector<double> > Q;


	/************** Member Methods **************/

	// Compute Q routine + overload
	void computeQ(int i, int j, int k, GridObj* g);
	void computeQ(int i, int j, GridObj* g);

};

#endif
