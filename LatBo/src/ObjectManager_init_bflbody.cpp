#include "..\inc\stdafx.h"
#include "..\inc\ObjectManager.h"

// Prefab body building routine
void ObjectManager::bfl_build_body(int body_type) {

	// This builder is like the IB body builder where we can use it to create prefab objects


}

// Overloaded builder for building from point cloud data
void ObjectManager::bfl_build_body(PCpts* _PCpts) {
	

	// Call body constructor and pass on pointer to hierarchy
	pBody.emplace_back(_PCpts, _Grids);

}