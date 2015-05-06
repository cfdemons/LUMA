#ifndef ENSIGHTGOLD_H
#define ENSIGHTGOLD_H

class EnsightGold {
	
	// Use default contructors and destructors

	public :

	// Generate case file
	void genCase(int nsteps, int saveEvery);
	// Generate geometry file
	void genGeo(int r);
	// Generate vectors file
	void genVec(int fileNum, int r);
	// Generate scalars file
	void genScal(int fileNum, int r);
} ;

#endif
