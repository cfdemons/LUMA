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
#ifndef MARDAT_H
#define MARDAT_H

#include "stdafx.h"

/// \brief	Container class to hold marker information.
class MarkerData {

public:

	/// \brief Constructor.
	/// \param i i-index of primary support site
	/// \param j j-index of primary support site
	/// \param k k-index of primary support site
	/// \param x x-position of marker
	/// \param y y-position of marker
	/// \param z z-position of marker
	/// \param ID marker number in a given body
	MarkerData(int i, int j, int k, double x, double y, double z, int ID)
		: i(i), j(j), k(k), x(x), y(y), z(z), ID(ID)
	{ };

	/// \brief	Default Constructor.
	///
	///			Initialise with invalid marker indicator which 
	///			is to set the x position to NaN.
	MarkerData(void) {

		// Give invalid ID number at construction time
		this->ID = -1;

	};

	/// Default destructor
	~MarkerData(void) {};

	/// Method to check if marker data structure is valid
	bool isValid()
	{
		if (ID == -1) return false;
		return true;
	}

	// Voxel indices
	int i;	///< i-index of primary support site
	int j;	///< j-index of primary support site
	int k;	///< k-index of primary support site

	/// Marker ID (position in array of markers)
	int ID;

	// Marker position
	double x;	///< x-position of marker
	double y;	///< y-position of marker
	double z;	///< z-position of marker

};

#endif // MARDAT_H