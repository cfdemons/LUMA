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
#ifndef IBMARKER_H
#define IBMARKER_H

#include "stdafx.h"
#include "Marker.h"

/// \brief	Immersed boundary marker.
///
///			This class declaration is for an immersed boundary Lagrange point.
///			A collection of these points form an immersed boundary body.
class IBMarker : public Marker {

	// Make ObjectManager a friend class so it can access the protected data of IBMarker objects
	friend class ObjectManager;

	// Make MPIManager a friend so it can access body data
	friend class MpiManager;

	// Same for IBBody
	friend class IBBody;
	friend class IBInfo;
	friend class FEMBody;

public:

	/// Default constructor
	IBMarker(void)
	{

		// TODO: Initialise the rest of the IB Marker properties here. 

	};

	/// Default destructor
	~IBMarker(void)
	{
	};

	// Custom constructor to add support etc.
	IBMarker(double xPos, double yPos, double zPos, int markerID, GridObj const * const body_owner);

protected:

	/************** Member Data **************/
	
	// General vectors
	std::vector<double> interpMom;		///< Fluid velocity interpolated from lattice nodes
	double interpRho;					///< Fluid density interpolate from lattice nodes
	std::vector<double> markerVel;		///< Desired velocity at marker
	std::vector<double> markerVel_km1;	///< Marker vel at last FSI sub-iteration
	std::vector<double> force_xyz;		///< Restorative force vector on marker
	std::vector<double> position0;	///< Vector containing the physical coordinates (x,y,z) of the marker at t-1. Used for moving bodies.

	// Support quantities
	std::vector<double> deltaval;		///< Value of delta function for a given support node

	// Scalars
	double epsilon;			///< Scaling parameter
	double local_area;		///< Area associated with support node in lattice units (same for all points if from same grid and regularly spaced like LBM)
	double dilation;		///< Dilation parameter in lattice units (same in all directions for uniform Eulerian grid)
	double ds;				///< Spacing of IBM marker in lattice units

	int owningRank;			///< Rank which owns the region where this marker exists
};

#endif
