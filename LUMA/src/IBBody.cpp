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

/* This file defines the constructors and methods for the immersed boundary body object.
*/

#include "../inc/stdafx.h"
#include "../inc/IBBody.h"
#include "../inc/IBMarker.h"
#include "../inc/PCpts.h"
#include "../inc/ObjectManager.h"

// *****************************************************************************
///	\brief	Default constructor for immersed boundary body
IBBody::IBBody()
{
	this->_Owner = nullptr;
	this->id = 0;
	fBody = NULL;
	isFlexible = false;
	isMovable = false;
}


// *****************************************************************************
///	\brief	Default destructor for immersed boundary body
IBBody::~IBBody(void)
{
}

// *****************************************************************************
///	\brief	Initialiser for a general IB body
///
///			This initialiser performs debugging IO and sets the property flags. 
///			It is called immediately after the appropriate constructor by all 
///			derived constructors.
///
///	\param	moveProperty	movable body flag passed on by the constructor
void IBBody::initialise(eMoveableType moveProperty)
{
	// Set local grid spacing
	dh = _Owner->dh;

	// Set FEM body pointer to NULL (if required it will be set properly later)
	fBody = NULL;

	// Set movable and flexible parameters
	if (moveProperty == eFlexible) {
		this->isFlexible = true;
		this->isMovable = true;
	}
	else if (moveProperty == eMovable) {
		this->isFlexible = false;
		this->isMovable = true;
	}
	else if (moveProperty == eRigid) {
		this->isFlexible = false;
		this->isMovable = false;
	}

	// Get rank
	int rank = GridUtils::safeGetRank();
	
	// Delete markers which exist off rank
	if (rank != owningRank) {
		*GridUtils::logfile << "Deleting IB markers which exist on receiver layer..." << std::endl;
		deleteRecvLayerMarkers();
	}

	// Get indices of valid markers
	getValidMarkers();
}


// *****************************************************************************
///	\brief	Get the indexes of the valid markers (only relevant for owning ranks in parallel)
void IBBody::getValidMarkers() {

	// Clear the vector
	validMarkers.clear();

	// Sort out the markers which do actually exist on this rank
#ifdef L_BUILD_FOR_MPI

	// Get rank
	int rank = GridUtils::safeGetRank();

	// If owning rank then not all markers are valid
	if (rank == owningRank) {

		// Declare values
		double x, y, z;
		eLocationOnRank loc;

		// Loop through all markers
		for (size_t m = 0; m < markers.size(); m++) {

			// Get positions
			x = markers[m].position[eXDirection];
			y = markers[m].position[eYDirection];
			z = markers[m].position[eZDirection];
			loc = eNone;

			// Check if on rank and not a receiver layer
			if (GridUtils::isOnThisRank(x, y, z, &loc, _Owner) && !GridUtils::isOnRecvLayer(x, y, z))
				validMarkers.push_back(static_cast<int>(m));
		}
	}
	else {

		// If not owning rank then all markers on this rank are valid
		validMarkers = GridUtils::onespace(0, static_cast<int>(markers.size()) - 1);
	}

#else
	// For serial all markers are valid
	validMarkers = GridUtils::onespace(0, static_cast<int>(markers.size()) - 1);
#endif
}


// *****************************************************************************
///	\brief	Point cloud marker IDs are not built in order - this step is required to sort them
void IBBody::sortPtCloudMarkers() {

	// If serial build then sorting is trivial
#ifndef L_BUILD_FOR_MPI

	// Loop through markers and assign value
	for (int i = 0; i < markers.size(); i++)
		markers[i].id = i;

#else

	// Parallel build is more tricky as owning rank needs to gather all markers, sort then scatter to other ranks
	MpiManager *mpim = MpiManager::getInstance();

	// Declare vectors for MPI comms
	std::vector<double> recvPositionBuffer;
	std::vector<int> recvIDBuffer;
	std::vector<int> recvSizeBuffer;
	std::vector<int> recvDisps;

	// Gather in marker IDs and positions
	mpim->mpi_ptCloudMarkerGather(this, recvPositionBuffer, recvIDBuffer, recvSizeBuffer, recvDisps);

	// Now do the sorting
	std::vector<int> sendSortedIDBuffer;
	if (mpim->my_rank == owningRank) {

		// Create vector for each marker which contains rank ID and an index array
		std::vector<int> rankIDs(recvIDBuffer.size(), 0);
		std::vector<int> indexIDs(recvIDBuffer.size(), 0);
		int count = 0;
		for (size_t i = 0; i < recvSizeBuffer.size(); i++) {
			for (int j = 0; j < recvSizeBuffer[i]; j++) {
				rankIDs[count] = static_cast<int>(i);
				indexIDs[count] = count;
				count++;
			}
		}

		// Sort the IDs according to the current ID values and return the indices
		std::sort(indexIDs.begin(), indexIDs.end(), [&](int a, int b) {return recvIDBuffer[a] < recvIDBuffer[b];});

		// Pack into send buffer
		sendSortedIDBuffer.resize(indexIDs.size(), 0);
		for (size_t i = 0; i < indexIDs.size(); i++) {
			sendSortedIDBuffer[indexIDs[i]] = static_cast<int>(i);
		}

		// First clear the markers
		markers.clear();

		// Recreate the markers
		double x, y, z;
		for (size_t i = 0; i < indexIDs.size(); i++) {

			// Get position
			x = recvPositionBuffer[indexIDs[i] * 3];
			y = recvPositionBuffer[indexIDs[i] * 3 + 1];
			z = recvPositionBuffer[indexIDs[i] * 3 + 2];

			// Add marker
			addMarker(x, y, z, static_cast<int>(i));
		}
	}

	// Gather in marker IDs and positions
	mpim->mpi_ptCloudMarkerScatter(this, sendSortedIDBuffer, recvSizeBuffer, recvDisps);
#endif
}


// *****************************************************************************
///	\brief	Custom constructor to populate body from array of points
///
///	\param 	g				hierarchy pointer to grid hierarchy
///	\param 	bodyID			global ID of body in array of bodies
///	\param 	_PCpts			pointer to point cloud data
///	\param 	moveProperty	determines if body is moveable, flexible or rigid
IBBody::IBBody(GridObj* g, int bodyID, PCpts* _PCpts, eMoveableType moveProperty)
	: Body(g, bodyID, _PCpts)
{

	// IBM-specific initialisation
	initialise(moveProperty);

	// Get rank
	int rank = GridUtils::safeGetRank();

	// Delete markers which exist off rank
	if (rank == owningRank) {
		*GridUtils::logfile << "Deleting IB markers which exist on receiver layer..." << std::endl;
		deleteRecvLayerMarkers();
	}

	// Sort the marker IDs as the point cloud reader does not build them consecutively
	sortPtCloudMarkers();

	// Get indices of valid markers
	getValidMarkers();
}


// *****************************************************************************
///	\brief	Custom constructor for building prefab filament
///
///	\param 	g					hierarchy pointer to grid hierarchy
///	\param 	bodyID				global ID of body in array of bodies
///	\param 	start_position		base of filament
/// \param	length				length of filament
/// \param	height				height of filament
/// \param	depth				depth of filament
/// \param	angles				angle of filament
///	\param 	moveProperty		determines if body is moveable, flexible or rigid
///	\param 	nElements			number of FEM elements
///	\param 	clamped				boundary condition for structural solver
///	\param 	density				material density
///	\param 	E					Young's modulus
IBBody::IBBody(GridObj* g, int bodyID, std::vector<double> &start_position,
		double length, double height, double depth, std::vector<double> &angles, eMoveableType moveProperty, int nElements, bool clamped, double density, double E)
		: Body(g, bodyID, start_position, length, angles)
{

	// IBM-specific initialisation
	initialise(moveProperty);

	// Get current rank
	int rank = GridUtils::safeGetRank();

	// If body is flexible then create an FEM body
	if (isFlexible == true && owningRank == rank)
		fBody = new FEMBody(this, start_position, length, height, depth, angles, nElements, clamped, density, E);
}


// *****************************************************************************
///	\brief	Custom constructor for building prefab circle/sphere
///
///	\param 	g					hierarchy pointer to grid hierarchy
///	\param 	bodyID				global ID of body in array of bodies
///	\param 	centre_point		centre of body
/// \param	radius				radius of body
///	\param 	moveProperty		determines if body is moveable, flexible or rigid
IBBody::IBBody(GridObj* g, int bodyID, std::vector<double> &centre_point,
		double radius, eMoveableType moveProperty)
		: Body(g, bodyID, centre_point, radius)
{
	// IBM-specific initialisation
	initialise(moveProperty);
}


// *****************************************************************************
///	\brief	Custom constructor for building square/cuboid
///
///	\param 	g					hierarchy pointer to grid hierarchy
///	\param 	bodyID				global ID of body in array of bodies
///	\param 	centre_point		centre of body
/// \param	dimensions			dimensions of body
/// \param	angles				angles of body
///	\param 	moveProperty		determines if body is moveable, flexible or rigid
IBBody::IBBody(GridObj* g, int bodyID, std::vector<double> &centre_point,
		std::vector<double> &dimensions, std::vector<double> &angles, eMoveableType moveProperty)
		: Body(g, bodyID, centre_point, dimensions, angles)
{
	// IBM-specific initialisation
	initialise(moveProperty);
}


// *****************************************************************************
///	\brief	Custom constructor for building plate
///
///	\param 	g					hierarchy pointer to grid hierarchy
///	\param 	bodyID				global ID of body in array of bodies
///	\param 	centre_point		centre of body
/// \param	length				length of body
/// \param	width				width of body
/// \param	angles				angles of body
///	\param 	moveProperty		determines if body is moveable, flexible or rigid
IBBody::IBBody(GridObj* g, int bodyID, std::vector<double> &centre_point,
		double length, double width, std::vector<double> &angles, eMoveableType moveProperty)
		: Body(g, bodyID, centre_point, length, width, angles)
{
	// IBM-specific initialisation
	initialise(moveProperty);
}
