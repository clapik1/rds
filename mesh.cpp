#include "mesh.h"
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
    triangles.resize(trianglesCount);
    for (int i = 0; i < trianglesCount; ++i) {
        ifs >> triangles[i].vertices[0] >> triangles[i].vertices[1] >> triangles[i].vertices[2];
        --triangles[i].vertices[0];
        --triangles[i].vertices[1];
        --triangles[i].vertices[2];
        for (int j = 0; j < 3; ++j) {
            triangles[i].lengths[j] = sqrt(pow(points[triangles[i].vertices[(j + 2) % 3]].x - points[triangles[i].vertices[(j + 1) % 3]].x, 2) + pow(points[triangles[i].vertices[(j + 2) % 3]].y - points[triangles[i].vertices[(j + 1) % 3]].y, 2));
            double a1 = points[triangles[i].vertices[(j + 2) % 3]].x - points[triangles[i].vertices[(j + 1) % 3]].x;
            double a2 = points[triangles[i].vertices[(j + 2) % 3]].y - points[triangles[i].vertices[(j + 1) % 3]].y;
            double b1 = points[triangles[i].vertices[j]].x - points[triangles[i].vertices[(j + 1) % 3]].x;
            double b2 = points[triangles[i].vertices[j]].y - points[triangles[i].vertices[(j + 1) % 3]].y;
            triangles[i].norm[j].x = a2 * (a2 * b1 - a1 * b2);
            triangles[i].norm[j].y = a1 * (a1 * b2 - a2 * b1);
            triangles[i].norm[j].normalize();
            triangles[i].norm[j] *= triangles[i].lengths[j];
        }
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

size_t mesh::getPointsCount() const {
    return pointsCount;
}

size_t mesh::getTrianglesCount() const {
    return trianglesCount;
}

size_t mesh::getWallsCount() const {
    return wallsCount;
}

const triangle2D &mesh::triangle(size_t nr) const {
    return triangles[nr];
}

const point2D &mesh::point(size_t nr) const {
    return points[nr];
}
