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

#ifndef IBMARKER_H
#define IBMARKER_H

#include "stdafx.h"
#include "Marker.h"

/// \brief	Immersed boundary marker
///
///			This class declaration is for an immersed boundary Lagrange point.
///			A collection of these points form an immersed boundary body.
class IBMarker : public Marker {

	/************** Friends **************/
	friend class ObjectManager;
	friend class MpiManager;
	friend class IBBody;
	friend class IBInfo;
	friend class FEMBody;
	friend class FEMElement;

public:

	/************** Constructors **************/
	IBMarker();
	~IBMarker();

	// Custom constructor to add support etc.
	IBMarker(double xPos, double yPos, double zPos, int markerID, GridObj const * const body_owner);

protected:

	/************** Member Data **************/
	
	// General vectors
	std::vector<double> interpMom;		///< Fluid velocity interpolated from lattice nodes
	double interpRho;					///< Fluid density interpolated from lattice nodes
	std::vector<double> markerVel;		///< Desired velocity at marker
	std::vector<double> markerVel_km1;	///< Marker velocity at last FSI sub-iteration
	std::vector<double> force_xyz;		///< Restorative force vector on marker
	std::vector<double> position0;		///< Vector containing the initial physical coordinates (x,y,z)

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
