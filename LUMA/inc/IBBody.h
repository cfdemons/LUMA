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

#ifndef IBBODY_H
#define IBBODY_H

#include "stdafx.h"

// Forward declarations
#include "Body.h"		// This is a templated class so include the whole file
#include "FEMBody.h"
class IBMarker;
class PCpts;
class GridObj;
class FEMBody;

/// \brief	Immersed boundary body class
///
///			Class for immersed boundary bodies (inherits from Body class with IB markers).
class IBBody : public Body<IBMarker> {

	/************** Friends **************/
	friend class ObjectManager;
	friend class IBInfo;
	friend class FEMBody;
	friend class FEMElement;
	friend class MpiManager;

public:

	/************** Constructors **************/
	IBBody(void);
	~IBBody(void);

	// Custom constructor which takes pointer to point cloud data and a pointer to the grid owner for the labelling
	IBBody(GridObj* g, int bodyID, PCpts* _PCpts, eMoveableType moveProperty);

	// Custom constructor for building prefab circle or sphere
	IBBody(GridObj* g, int bodyID, std::vector<double> &centre_point, double radius, eMoveableType moveProperty);

	// Custom constructor for building prefab square or cuboid
	IBBody(GridObj* g, int bodyID, std::vector<double> &centre_point,
		std::vector<double> &dimensions, std::vector<double> &angles, eMoveableType moveProperty);

	// Custom constructor for building prefab plate
	IBBody(GridObj* g, int bodyID, std::vector<double> &centre_point,
		double length, double width, std::vector<double> &angles, eMoveableType moveProperty);

	// Custom constructor for building prefab filament
	IBBody(GridObj* g, int bodyID, std::vector<double> &start_position,
		double length, double height, double depth, std::vector<double> &angles, eMoveableType moveProperty,
		int nElement, bool clamped, double density, double E);



protected:

	/************** Member Data **************/

	bool isFlexible;					///< Flag to indicate flexibility: false == rigid body; true == flexible filament
	bool isMovable;						///< Flag to indicate if body is movable or not.

	double dh;							///< Local grid spacing for this body

	FEMBody *fBody;						///< Pointer to FEM body object


	/************** Member Methods **************/

private:

	void initialise(eMoveableType moveProperty);		// Initialisation wrapper for setting flags
	void getValidMarkers();								// Get valid markers in iBody (only relevant for owning rank)
	void sortPtCloudMarkers();							// Sort pt cloud markers and IDs

};

#endif
