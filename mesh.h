#pragma once

#include <string>
#include <vector>
#include "point2D.h"

struct wall2D {
	int nr;
	int vertices[2];
};

struct triangle2D {
	int vertices[3];
};

class mesh {
public:
	bool init(std::istream &ifs);
	int getPointsCount();
	int getTrianglesCount();
	int getWallsCount();
	double getTriangleArea(int triangleNr);
	int pointsCount;
	std::vector<point2D> points;
	int trianglesCount;
	std::vector<int> triangles[3];
	int wallsCount;
	std::vector<wall2D> walls;
};
