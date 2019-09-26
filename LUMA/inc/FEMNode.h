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

#ifndef FEMNODE_H
#define FEMNODE_H

/// \brief	Finite element node class
///
///			Class for the finite element nodes which are contained within the FEMBody.
class FEMNode {

	/************** Friends **************/
	friend class FEMBody;
	friend class FEMElement;
	friend class ObjectManager;

	/************** Constructors **************/
public:
	FEMNode();
	~FEMNode();

	// Custom constructor for building FEM node from input parameters
	FEMNode(int idx, double x, double y, double z, std::vector<double> &angles);


	/************** Member Data **************/

	// Node values
private:
	int ID;								///< Node ID in FEM body
	std::vector<double> position0;		///< Initial position of FEM node
	std::vector<double> position;		///< Current position of FEM node
	double angles0;						///< Initial angles of FEM node
	double angles;						///< Current angles of FEM node


	/************** Member Methods **************/

};

#endif
