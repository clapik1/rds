#include "mesh.h"
#include "point2D.h"
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>

bool mesh::init(std::istream &ifs) {
	std::string s;
	ifs >> s >> s >> s >> pointsCount;
	points.resize(pointsCount);
	for (int i = 0; i < pointsCount; ++i) {
		ifs >> points[i].x >> points[i].y;
	}
	ifs >> s >> s >> trianglesCount;
	triangles[0].resize(trianglesCount);
	triangles[1].resize(trianglesCount);
	triangles[2].resize(trianglesCount);
	for (int i = 0; i < trianglesCount; ++i) {
		ifs >> triangles[0][i] >> triangles[1][i] >> triangles[2][i];
		--triangles[0][i];
		--triangles[1][i];
		--triangles[2][i];
	}
	ifs >> s >> s >> wallsCount;
	walls.resize(wallsCount);
	for (int i = 0; i < wallsCount; ++i) {
		ifs >> walls[i].nr >> walls[i].vertices[0] >> walls[i].vertices[1] >> s;
		--walls[i].vertices[0];
		--walls[i].vertices[1];
	}
	return true;
}

int mesh::getPointsCount() {
	return pointsCount;
}

int mesh::getTrianglesCount() {
	return trianglesCount;
}

int mesh::getWallsCount() {
	return wallsCount;
}

double mesh::getTriangleArea(int triangleNr) {
	double a = sqrt(pow(points[triangles[0][triangleNr]].x - points[triangles[1][triangleNr]].x, 2) + pow(points[triangles[0][triangleNr]].y - points[triangles[1][triangleNr]].y, 2));
	double b = sqrt(pow(points[triangles[0][triangleNr]].x - points[triangles[2][triangleNr]].x, 2) + pow(points[triangles[0][triangleNr]].y - points[triangles[2][triangleNr]].y, 2));
	double c = sqrt(pow(points[triangles[2][triangleNr]].x - points[triangles[1][triangleNr]].x, 2) + pow(points[triangles[2][triangleNr]].y - points[triangles[1][triangleNr]].y, 2));
	double p = (a + b + c) / 2;
	//std::cout << sqrt(p * (p - a) * (p - b) * (p - c)) << std::endl;
	return sqrt(p * (p - a) * (p - b) * (p - c));
}
