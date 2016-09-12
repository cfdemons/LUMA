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

#pragma once
#include "BFLMarker.h"
#include "Body.h"
#include "PCpts.h"
#include "ObjectManager.h"

/// \brief	BFL body.
///
///			A BFL body is made up of a collection of BFLMarkers.
class BFLBody :
	public Body<BFLMarker>
{

	friend class GridObj;

public:
	// Default constructor and destructor
	BFLBody(void);
	~BFLBody(void);
	// Custom constructor which takes pointer to point cloud data and a pointer to the grid hierarchy for the labelling
	BFLBody(PCpts *_PCpts, GridObj *g_hierarchy);

protected:

	/*
	***************************************************************************************************************
	********************************************* Member Data *****************************************************
	***************************************************************************************************************
	*/

	/// \brief	Distance between adjacent lattice site and the surface of the body.
	///
	///			There are two stores of values. Store 1 is the distance on one 
	///			side of the wall and store 2 the distance on the other side. 
	///			One store is appended to the other in this structure.
	std::vector< std::vector<double> > Q;


	/*
	***************************************************************************************************************
	********************************************* Member Methods **************************************************
	***************************************************************************************************************
	*/

	// Compute Q routine + overload
	void computeQ(int i, int j, int k, int N_lim, int M_lim, int K_lim, GridObj* g);
	void computeQ(int i, int j, int N_lim, int M_lim, GridObj* g);

};