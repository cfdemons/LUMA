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


public:
	/// Default constructor
	Marker(void)
	{
		position.push_back(0.0);
		position.push_back(0.0);
		position.push_back(0.0);

		supp_i.push_back(0);
		supp_j.push_back(0);
		supp_k.push_back(0);

#ifdef L_BUILD_FOR_MPI
		support_rank.push_back(MpiManager::getInstance()->my_rank);
#else
		support_rank.push_back(0);
#endif

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
	///	\param	body_owner	Grid on which primary support point is to be found.
	Marker(double x, double y, double z, GridObj const * const body_owner)
	{
		position.push_back(x);
		position.push_back(y);
		position.push_back(z);

		std::vector<int> ijk;
		GridUtils::getEnclosingVoxel(x, y, z, body_owner, &ijk);
		supp_i.push_back(ijk[0]);
		supp_j.push_back(ijk[1]);
		supp_k.push_back(ijk[2]);

#ifdef L_BUILD_FOR_MPI
		support_rank.push_back(MpiManager::getInstance()->my_rank);
#else
		support_rank.push_back(0);
#endif
	}

};

#endif