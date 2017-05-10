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
#ifndef MARKER_H
#define MARKER_H

#include "stdafx.h"

/// \brief	Generic marker class.
class Marker
{

	// Members //
public:
	std::vector<double> position;		///< Position vector of marker location in physical units

	/* Vector of indices for lattice sites considered to be in support of the marker:
	* In IBM these are all the lattice support sites for spreading and interpolating;
	* In BFL these are the voxel indices in which the BFL marker resides;
	*/
	std::vector<int> supp_i;	///< X-indices of lattice sites in support of this marker
	std::vector<int> supp_j;	///< Y-indices of lattice sites in support of this marker
	std::vector<int> supp_k;	///< Z-indices of lattice sites in support of this marker

	/// Array of indices indicating on which rank the given support point resides
	std::vector<int> support_rank;

	/// ID of marker within its owning body
	int id;

	/// Spacing between this marker and neighbours
	double ds;


public:
	/// Default constructor
	Marker(void)
	{
		// Set marker position to zero by default
		position.push_back(0.0);
		position.push_back(0.0);
		position.push_back(0.0);

		// Set the i, j, and k for the closest voxel to zero
		supp_i.push_back(0);
		supp_j.push_back(0);
		supp_k.push_back(0);

		// Set rank of first support marker and marker ID in body
		support_rank.push_back(GridUtils::safeGetRank());
		id = 0;
		ds = 0.0;
	};

	/// Default destructor
	~Marker(void)
	{
	};

	/// \brief Custom constructor which locates marker.
	///
	///			In order to properley initialise during construction, a
	///			grid should be passed in on which primary support point
	///			can be found.
	///
	/// \param	x			X-position of marker
	/// \param	y			Y-position of marker
	/// \param	z			Z-position of marker
	/// \param markerID		ID of marker within body
	///	\param	body_owner	Grid on which primary support point is to be found.
	Marker(double x, double y, double z, int markerID, GridObj const * const body_owner)
	{

		// Set marker position
		position.push_back(x);
		position.push_back(y);
		position.push_back(z);

		// Get the i, j, and k for the closest voxel to markers
		std::vector<int> ijk;
		GridUtils::getEnclosingVoxel(x, y, z, body_owner, &ijk);
		supp_i.push_back(ijk[0]);
		supp_j.push_back(ijk[1]);
		supp_k.push_back(ijk[2]);

		// Set rank of first support marker and marker ID in body
		support_rank.push_back(GridUtils::safeGetRank());
		id = markerID;
		ds = 0.0;
	}
};

#endif
