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

/* This file defines the constructors and methods for the immersed boundary body object.
*/

#include "../inc/stdafx.h"
#include "../inc/IBBody.h"
#include "../inc/IBMarker.h"
#include "../inc/PCpts.h"
#include "../inc/ObjectManager.h"

// ***************************************************************************************************
/// \brief Constructor which sets group ID to zero by default.
IBBody::IBBody()
{
	this->_Owner = nullptr;
	this->id = 0;
}

/// Default destructor
IBBody::~IBBody(void)
{
}

/******************************************************************************/
/// \brief	Initialiser for a general IB body.
///
///			This initialiser performs debugging IO and sets the property flags. 
///			It is called immediately after the appropriate constructor by all 
///			derived constructors.
///
///	\param	moveProperty	movable body flag passed on by the constructor.
void IBBody::initialise(eMoveableType moveProperty)
{

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


/*********************************************/
/// \brief 		Get the indexes of the valid markers (only relevant for owning ranks in parallel)
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
		for (int m = 0; m < markers.size(); m++) {

			// Get positions
			x = markers[m].position[eXDirection];
			y = markers[m].position[eYDirection];
			z = markers[m].position[eZDirection];
			loc = eNone;

			// Check if on rank and not a receiver layer
			if (GridUtils::isOnThisRank(x, y, z, &loc, _Owner) && !GridUtils::isOnRecvLayer(x, y, z))
				validMarkers.push_back(m);
		}
	}
	else {

		// If not owning rank then all markers on this rank are valid
		validMarkers = GridUtils::onespace(0, markers.size()-1);
	}

#else
	// For serial all markers are valid
	validMarkers = GridUtils::onespace(0, markers.size()-1);
#endif
}


/*********************************************/
/// \brief Custom constructor to populate body from array of points.
/// \param g				hierarchy pointer to grid hierarchy
/// \param bodyID			ID of body in array of bodies.
/// \param _PCpts			pointer to point cloud data
/// \param moveProperty		determines if body is moveable, flexible or rigid
/// \param clamped			boundary condition for fixed end (only relevant if body is flexible)
IBBody::IBBody(GridObj* g, int bodyID, PCpts* _PCpts, eMoveableType moveProperty, bool clamped)
	: Body(g, bodyID, _PCpts)
{

	// IBM-specific initialisation
	initialise(moveProperty);

	// Sort the marker IDs as the point cloud reader does not build them consecutively
	sortMarkerIDs();
}


/*********************************************/
/// \brief 	Custom constructor for building prefab filament
/// \param g				hierarchy pointer to grid hierarchy
/// \param bodyID			ID of body in array of bodies.
/// \param start_position	start position of base of filament
/// \param length			length of filament
/// \param angles			angle of filament
/// \param moveProperty		determines if body is moveable, flexible or rigid
/// \param clamped			boundary condition for fixed end (only relevant if body is flexible)
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


/*********************************************/
/// \brief 	Custom constructor for building prefab circle/sphere
/// \param g				hierarchy pointer to grid hierarchy
/// \param bodyID			ID of body in array of bodies.
/// \param centre_point		centre point of circle
/// \param radius			radius of circle
/// \param moveProperty		determines if body is moveable, flexible or rigid
IBBody::IBBody(GridObj* g, int bodyID, std::vector<double> &centre_point,
		double radius, eMoveableType moveProperty)
		: Body(g, bodyID, centre_point, radius)
{
	// IBM-specific initialisation
	initialise(moveProperty);
}

/*********************************************/
/// \brief 	Custom constructor for building square/cuboid
/// \param g				hierarchy pointer to grid hierarchy
/// \param bodyID			ID of body in array of bodies.
/// \param centre_point		centre point of square
/// \param dimensions		dimensions of square
/// \param angles			angle of square
/// \param moveProperty		determines if body is moveable, flexible or rigid
IBBody::IBBody(GridObj* g, int bodyID, std::vector<double> &centre_point,
		std::vector<double> &dimensions, std::vector<double> &angles, eMoveableType moveProperty)
		: Body(g, bodyID, centre_point, dimensions, angles)
{
	// IBM-specific initialisation
	initialise(moveProperty);
}


/*********************************************/
/// \brief 		Custom constructor for building dummy iBody for epsilon calculation
/// \param iBody				original iBody that this rank owns
/// \param recvBuffer			all additional information required for epsilon calculation
IBBody::IBBody(IBBody &iBody, std::vector<std::vector<double>> &recvBuffer) {

	// Get MPI manager instance
	MpiManager *mpim = MpiManager::getInstance();

	// Set grid pointer
	_Owner = iBody._Owner;

	// Set body ID
	id = iBody.id;

	// Set spacing
	spacing = iBody.spacing;

	// First size markers from original body
	int nMarkers = iBody.markers.size();

	// Now get extra markers that have been received from other ranks
	for (int i = 0; i < mpim->markerCommOwnerSide.size(); i++) {
		if (mpim->markerCommOwnerSide[i].bodyID == iBody.id)
			nMarkers++;
	}

	// Now resize the markers
	markers.resize(nMarkers);

	// Now assign values from the original body first
	int markerID;
	for (int m = 0; m < iBody.markers.size(); m++) {

		// Get markerID
		markerID = iBody.markers[m].id;

		// Assign the position of the marker and first support site
		markers[markerID].position[eXDirection] = iBody.markers[m].position[eXDirection];
		markers[markerID].position[eYDirection] = iBody.markers[m].position[eYDirection];
		markers[markerID].supp_x[0] = iBody.markers[m].supp_x[0];
		markers[markerID].supp_y[0] = iBody.markers[m].supp_y[0];
		markers[markerID].deltaval.push_back(iBody.markers[m].deltaval[0]);
#if (L_DIMS == 3)
		markers[markerID].position[eZDirection] = iBody.markers[m].position[eZDirection];
		markers[markerID].supp_z[0] = iBody.markers[m].supp_z[0];
#endif
		markers[markerID].local_area = iBody.markers[m].local_area;
		markers[markerID].dilation = iBody.markers[m].dilation;

		// Assign rest of the support sites
		for (int s = 1; s < iBody.markers[m].deltaval.size(); s++) {
			markers[markerID].supp_x.push_back(iBody.markers[m].supp_x[s]);
			markers[markerID].supp_y.push_back(iBody.markers[m].supp_y[s]);
#if (L_DIMS == 3)
			markers[markerID].supp_z.push_back(iBody.markers[m].supp_z[s]);
#endif
			markers[markerID].deltaval.push_back(iBody.markers[m].deltaval[s]);
		}
	}

	// Index vector for looping through recvBuffer
	std::vector<int> idx(mpim->num_ranks, 0);

	// Get starting indices
	for (int i = 0; i < mpim->markerCommOwnerSide.size(); i++) {
		if (mpim->markerCommOwnerSide[i].bodyID < iBody.id)
			idx[mpim->markerCommOwnerSide[i].rankComm] += L_DIMS + 2 + mpim->markerCommOwnerSide[i].nSupportSites * (L_DIMS + 1);
	}

	// Rank to collect data from
	int fromRank;

	// Now fill in the rest of the values from other ranks
	for (int i = 0; i < mpim->markerCommOwnerSide.size(); i++) {
		if (mpim->markerCommOwnerSide[i].bodyID == iBody.id) {

			// Get ID info
			fromRank = mpim->markerCommOwnerSide[i].rankComm;
			markerID = mpim->markerCommOwnerSide[i].markerID;

			// Assign the position of the marker
			markers[markerID].position[eXDirection] = recvBuffer[fromRank][idx[fromRank]]; idx[fromRank]++;
			markers[markerID].position[eYDirection] = recvBuffer[fromRank][idx[fromRank]]; idx[fromRank]++;
#if (L_DIMS == 3)
			markers[markerID].position[eZDirection] = recvBuffer[fromRank][idx[fromRank]]; idx[fromRank]++;
#endif
			markers[markerID].local_area = recvBuffer[fromRank][idx[fromRank]]; idx[fromRank]++;
			markers[markerID].dilation = recvBuffer[fromRank][idx[fromRank]]; idx[fromRank]++;

			// Assign first support site
			markers[markerID].supp_x[0] = recvBuffer[fromRank][idx[fromRank]]; idx[fromRank]++;
			markers[markerID].supp_y[0] = recvBuffer[fromRank][idx[fromRank]]; idx[fromRank]++;
#if (L_DIMS == 3)
			markers[markerID].supp_z[0] = recvBuffer[fromRank][idx[fromRank]]; idx[fromRank]++;
#endif
			markers[markerID].deltaval.push_back(recvBuffer[fromRank][idx[fromRank]]); idx[fromRank]++;

			// Loop through rest of support sites
			for (int s = 1; s < mpim->markerCommOwnerSide[i].nSupportSites; s++) {
				markers[markerID].supp_x.push_back(recvBuffer[fromRank][idx[fromRank]]); idx[fromRank]++;
				markers[markerID].supp_y.push_back(recvBuffer[fromRank][idx[fromRank]]); idx[fromRank]++;
	#if (L_DIMS == 3)
				markers[markerID].supp_z.push_back(recvBuffer[fromRank][idx[fromRank]]); idx[fromRank]++;
	#endif
				markers[markerID].deltaval.push_back(recvBuffer[fromRank][idx[fromRank]]); idx[fromRank]++;
			}
		}
	}
}
