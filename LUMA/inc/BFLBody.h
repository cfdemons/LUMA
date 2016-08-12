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

#pragma once
#include "BFLMarker.h"
#include "Body.h"
#include "PCpts.h"

// Class representing return structure for marker query
class MarkerData {

public:

	MarkerData(int i, int j, int k, double x, double y, double z, int ID) {
	
		// Custom constructor
		this->i = i;
		this->j = j;
		this->k = k;
		this->x = x;
		this->y = y;
		this->z = z;
		this->ID = ID;	
	
	};
	MarkerData(void) {
	
		// Essentially null a double in the store making it invalid
		this->x = std::numeric_limits<double>::quiet_NaN();
	
	};
	~MarkerData(void) {};

	// Voxel indices
	int i;
	int j;
	int k;

	// Marker ID (position in array of markers)
	int ID;

	// Marker position
	double x;
	double y;
	double z;

};



/************************************************
 *********** BFL Body implementation ************
 ***********************************************/

// Definition of the BFL body class representing a BFL body made up of BFLMarkers //
class BFLBody :
	public Body<BFLMarker>
{
	
	// Allow markers, Grid and ObjectManager access to static utilities
	friend class BFLMarker;
	friend class GridObj;
	friend class ObjectManager;

public:
	// Default constructor and destructor
	BFLBody(void);
	~BFLBody(void);
	// Custom constructor which takes pointer to point cloud data and a pointer to the grid hierarchy for the labelling
	BFLBody(PCpts* _PCpts, GridObj* g);

protected:

	/*
	***************************************************************************************************************
	********************************************* Member Data *****************************************************
	***************************************************************************************************************
	*/

	// Q values (append store 2 onto store 1)
	std::vector< std::vector<double> > Q;


	/*
	***************************************************************************************************************
	********************************************* Member Methods **************************************************
	***************************************************************************************************************
	*/

	// Marker adder (pass the marker index and the counter by reference so they are updated in the calling scope)
	void bflMarkerAdder(double x, double y, double z, int& curr_mark, std::vector<int>& counter);

	// Compute Q routine + overload
	void computeQ(int i, int j, int k, int N_lim, int M_lim, int K_lim, GridObj* g);
	void computeQ(int i, int j, int N_lim, int M_lim, GridObj* g);

	// Utility functions (all static)
	static std::vector<int> getVoxInd(double x, double y, double z);
	static bool isInVoxel(double x, double y, double z, int curr_mark, BFLBody* body);
	static bool isVoxelBflVoxel(double x, double y, double z, BFLBody* body);
	static MarkerData* getMarkerData(double x, double y, double z, BFLBody* body);
public :
	static int getVoxInd(double p);

};