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

	// Delete any markers which are on the receiver layer as not needed for IBM
#ifdef L_BUILD_FOR_MPI
	*GridUtils::logfile << "Deleting IB markers which exist on receiver layer..." << std::endl;
	deleteRecvLayerMarkers();
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
	initialise(moveProperty);	
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
		double length, std::vector<double> &angles, eMoveableType moveProperty, bool clamped)
		: Body(g, bodyID, start_position, length, angles)
{
	initialise(moveProperty);
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
	initialise(moveProperty);
}


/*********************************************/
/// \brief 		Custom constructor for building dummy iBody for epsilon calculation
/// \param iBody				original iBody that this rank owns
/// \param recvBuffer			all additional information required for epsilon calculation
IBBody::IBBody(IBBody &iBody, std::vector<std::vector<double>> &recvBuffer) {

	// Get MPI manager instance
	MpiManager *mpim = MpiManager::getInstance();

	// Set body ID
	id = iBody.id;

	// First size markers from original body
	int nMarkers = iBody.markers.size();

	// Now get extra markers that have been received from other ranks
	for (int i = 0; i < mpim->epsCommOwnerSide.size(); i++) {
		if (mpim->epsCommOwnerSide[i].bodyID == iBody.id)
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
	for (int i = 0; i < mpim->epsCommOwnerSide.size(); i++) {
		if (mpim->epsCommOwnerSide[i].bodyID < iBody.id)
			idx[mpim->epsCommOwnerSide[i].rankComm] += L_DIMS + 2 + mpim->epsCommOwnerSide[i].nSupportSites * (L_DIMS + 1);
	}

	// Rank to collect data from
	int fromRank;

	// Now fill in the rest of the values from other ranks
	for (int i = 0; i < mpim->epsCommOwnerSide.size(); i++) {
		if (mpim->epsCommOwnerSide[i].bodyID == iBody.id) {

			// Get ID info
			fromRank = mpim->epsCommOwnerSide[i].rankComm;
			markerID = mpim->epsCommOwnerSide[i].markerID;

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
			for (int s = 1; s < mpim->epsCommOwnerSide[i].nSupportSites; s++) {
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
