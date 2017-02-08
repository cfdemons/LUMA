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

#include "../inc/stdafx.h"
#include "../inc/ObjectManager.h"

/// \brief	Prefab body building routine.
///
///			Not implemented in this version.
///
/// \param	body_type	type of prefab body to be built.
void ObjectManager::bfl_buildBody(int body_type) {

	// This builder is like the IB body builder where we can use it to create prefab objects


}

/// \brief	Wrapper for building BFL body from point cloud.
/// \param	_PCpts	pointer to point cloud data.
void ObjectManager::bfl_buildBody(PCpts* _PCpts) {	

	// Call body constructor and pass on pointer to hierarchy
	pBody.emplace_back(_Grids, pBody.size(), _PCpts);

}