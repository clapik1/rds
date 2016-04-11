#pragma once

#include <string>
#include <vector>
#include "point2D.h"
#include "vector2D.h"
#include "triangle2D.h"

struct wall2D {
	int nr;
	int vertices[2];
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
	std::vector<triangle2D> triangles;
	int wallsCount;
	std::vector<wall2D> walls;
};
