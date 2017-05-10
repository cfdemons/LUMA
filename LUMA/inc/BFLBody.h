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

	// Custom constructor for building prefab circle or sphere
	BFLBody(GridObj* g, int bodyID, std::vector<double> &centre_point, double radius);

	// Custom constructor for building prefab filament
	BFLBody(GridObj* g, int bodyID, std::vector<double> &centre_point, std::vector<double> &dimensions, std::vector<double> &angles);

	// Custom constructor for building prefab square or cuboid
	BFLBody(GridObj* g, int bodyID, std::vector<double> &start_position, double length, std::vector<double> &angles);

protected:

	/************** Member Data **************/

	/// \brief	Distance between adjacent lattice site and the surface of the body.
	///
	///			Q values are stored in a flattened 2D array with the fatest 
	///			changing index the velocity direction and the second the marker
	///			ID.
	std::vector<double>	Q;


	/************** Member Methods **************/
private :

	// Initialiser (wrapper for labeller and Q computation)
	void initialise();

	// Compute Q routine + overload
	void computeQ(int i, int j, int k, GridObj* g);
	void computeQ(int i, int j, GridObj* g);

	// Surface closure
	void enforceSurfaceClosure();

};

#endif
