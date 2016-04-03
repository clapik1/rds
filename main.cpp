#include <iostream>
#include <fstream>
#include "mesh.h"

int main(int argc, char *argv[]) {
	std::ifstream ifs(argv[1], std::ifstream::in);
	//std::ifstream ifs("D:/studia/przejsciowka/meshes/grid_8.msh2", std::ifstream::in);
	mesh mMesh;
	if (!mMesh.init(ifs)) {
		std::cout << "error initializing mesh" << std::endl;
		return 0;
	}
	double area = 0;
	for (int i = 0; i < mMesh.getTrianglesCount(); ++i) {
		area += mMesh.getTriangleArea(i);
	}
	std::cout << area << std::endl;
	mMesh.points[2].value = 1;
	/*for (int i = 0; i < 10; ++i) {
		for (int j = 0; j < mMesh.pointsCount; ++j) {
			std::cout << mMesh.points[j].value << ' ';
		}
		std::cout << std::endl;
		//iterate
	}*/
	return 0;
}
