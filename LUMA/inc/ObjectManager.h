/*
* --------------------------------------------------------------
*
* ------ Lattice Boltzmann @ The University of Manchester ------
*
* -------------------------- L-U-M-A ---------------------------
*
* Copyright 2018 The University of Manchester
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

#ifndef OBJMAN_H
#define OBJMAN_H

#include "stdafx.h"
#include "IVector.h"
#include "IBInfo.h"
#include "IBMarker.h"
#include "IBBody.h"
#include "BFLBody.h"

class PCpts;
class GridObj;

/// \brief	Object Manager class.
///
///			Class to manage all objects in the domain from creation through 
///			manipulation to destruction.
class ObjectManager
{

	// Make Grid a friend so boundary conditions can access the body data
	friend class GridObj;
	friend class MpiManager;

	/// \brief	Nested Geometry data structure class.
	///
	///			Stores information about a body in a coherent manner but not 
	///			intended to be instantiated outside the Object Manager.
	class GeomPacked
	{
	public:

		GeomPacked();
		GeomPacked(
			eObjectType objtype, int bodyID, std::string fileName, 
			int onGridLev, int onGridReg,
			bool isCentreX, double refX,
			bool isCentreY, double refY, 
			bool isCentreZ, double refZ, 
			double bodyLength, eCartesianDirection scaleDirection, 
			eMoveableType moveProperty, bool isClamped
			);
		~GeomPacked();

		static bool interpretRef(std::string refType);
		
	
		// Members
		eObjectType objtype;
		int bodyID;
		std::string fileName;
		int onGridLev;
		int onGridReg;
		bool isRefXCentre;
		bool isRefYCentre;
		bool isRefZCentre;
		double bodyRefX;
		double bodyRefY;
		double bodyRefZ;
		double bodyLength;
		eCartesianDirection scaleDirection;
		eMoveableType moveProperty;
		bool isClamped;

	};

	/* Members */

private:

	// Private file stream for debugging momentum exchange
	std::ofstream debugstream;

	// Bounce-back object fields
	double bbbForceOnObjectX = 0.0;			///< Instantaneous X-direction force on BB bodies in domain
	double bbbForceOnObjectY = 0.0;			///< Instantaneous Y-direction force on BB bodies in domain
	double bbbForceOnObjectZ = 0.0;			///< Instantaneous Z-direction force on BB bodies in domain
	int bbbOnGridLevel = -1;				///< Grid level on which the BB body resides
	int bbbOnGridReg = -1;					///< Grid region on which the BB body resides

	// Objects (could be stored in a single Body array if we use pointers)
	std::vector<IBBody> iBody;				///< Array of immersed boundary bodies
	std::vector<BFLBody> pBody;				///< Array of BFL bodies

	/// Pointer to grid hierarchy
	GridObj* _Grids;

	/// Pointer to self
	static ObjectManager* me;

	// Flag for if there are any flexible bodies in the simulation
	std::vector<bool> hasIBMBodies;
	std::vector<bool> hasFlexibleBodies;

	// Map global body ID to an index in the iBody vector
	std::vector<int> bodyIDToIdx;

	// Vector of indices for iBody vector for which this rank owns and is flexible
	std::vector<int> idxFEM;

	// Subiteration loop parameters
	double timeav_subResidual;
	double timeav_subIterations;

	/* Methods */

private:
	ObjectManager(void);		///< Private constructor
	~ObjectManager(void);		///< Private destructor
	/// Overloaded constructor to set pointer to grid hierarchy
	ObjectManager(GridObj* _Grids);		

public:
	// Singleton design
	static ObjectManager *getInstance();			///< Get instance method
	static void destroyInstance();					///< Destroy instance method
	static ObjectManager *getInstance(GridObj* g);	///< Overloaded get instance passing in pointer to grid hierarchy

	// IBM methods //
	void ibm_apply(GridObj *g, bool doSubIterate);									// Apply interpolate, compute and spread operations for all bodies.
	void ibm_initialise();															// Initialise a built immersed body with support.
	double ibm_deltaKernel(double rad, double dilation);							// Evaluate kernel (delta function approximation).
	void ibm_interpolate(int level);												// Interpolation of velocity field onto markers of ib-th body.
	void ibm_spread(int level);														// Spreading of restoring force from ib-th body.
	void ibm_updateMacroscopic(int level);											// Update the macroscopic values with the IBM force
	void ibm_findSupport(int ib);													// Populates support information for the m-th marker of ib-th body.
	void ibm_initialiseSupport(int ib, int m, std::vector<double> &estimated_position);	// Initialises data associated with the support points.
	void ibm_computeForce(int level);												// Compute restorative force at each marker in ib-th body.
	void ibm_findEpsilon(int level);												// Method to find epsilon weighting parameter for ib-th body.
	void ibm_computeDs(int level);
	void ibm_moveBodies(int level);													// Update all IBBody positions and support.
	void ibm_finaliseReadIn(int iBodyID);											// Do some house-keeping after geometry read in
	void ibm_universalEpsilonGather(int level, IBBody &iBodyTmp);					// Gather all the markers into the temporary iBody vector
	void ibm_universalEpsilonScatter(int level, IBBody &iBodyTmp);					// Gather all the markers into the temporary iBody vector
	void ibm_subIterate(GridObj *g);												// Subiterate to enforce correct kinematic conditions at interface
	double ibm_checkVelDiff(int level);												// Check residual from sub-iteration step

	// IBM Debug methods //
	void ibm_debug_epsilon(int ib);
	void ibm_debug_interpVel(int ib);
	void ibm_debug_markerForce(int ib);
	void ibm_debug_markerPosition(int ib);
	void ibm_debug_supportInfo(int ib);
	void ibm_debug_supportVel(int ib);
	void ibm_debug_supportForce(int ib);

	// IBM-MPI methods
	void ibm_updateMPIComms(int level);
	void ibm_interpolateOffRankVels(int level);
	void ibm_spreadOffRankForces(int level);
	void ibm_updateMarkers(int level);

	// Bounceback Body Methods
	void addBouncebackObject(GeomPacked *geom, PCpts *_PCpts);				// Override method to add BBB from cloud reader.
	void addBouncebackObject(GridObj *g, GeomPacked *geom, PCpts *_PCpts);	// Method to add a BBB from the cloud reader.
	void computeLiftDrag(int i, int j, int k, GridObj *g);			// Compute force using Momentum Exchange for BBB on supplied grid.
	void computeLiftDrag(int v, int id, GridObj *g, int markerID);	// Compute force using Momentum Exchange for BFL on supplied grid.
	void resetMomexBodyForces(GridObj * grid);						// Reset the force stores for Momentum Exchange

	// IO methods //
	void io_vtkBodyWriter(int tval);						// VTK body writer wrapper
	void io_vtkFEMWriter(int tval);							// VTK FEM writer
	void io_writeBodyPosition(int timestep);				// Write out IBBody positions at specified timestep to text files
	void io_writeLiftDrag();								// Write out IBBody lift and drag at specified timestep
	void io_restart(eIOFlag IO_flag, int level);			// Restart read and write for IBBodies given grid level
	void io_readInCloud(PCpts*& _PCpts, GeomPacked *geom);	// Method to read in Point Cloud data
	void io_writeForcesOnObjects(double tval);				// Method to write object forces to a csv file
	void io_readInGeomConfig();								// Read in geometry configuration file
	void io_writeTipPositions(int t);						// Write out tip positions of flexible filaments

	// Debug
	void toggleDebugStream(GridObj *g);		// Method to open/close a debugging file
};

#endif
