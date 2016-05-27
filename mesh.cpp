#include "mesh.h"
#include <fstream>
#include <iostream>
#include <cmath>

mesh::mesh(std::istream &ifs) {
    std::string s;
    ifs >> s >> s >> s >> pointsCount;
    points.resize(pointsCount);
    for (size_t i = 0; i < pointsCount; ++i) {
        ifs >> points[i].x >> points[i].y;
    }
    ifs >> s >> s >> trianglesCount;
    triangles.resize(trianglesCount);
    for (size_t i = 0; i < trianglesCount; ++i) {
        ifs >> triangles[i].vertices[0] >> triangles[i].vertices[1] >> triangles[i].vertices[2];
        --triangles[i].vertices[0];
        --triangles[i].vertices[1];
        --triangles[i].vertices[2];
        for (size_t j = 0; j < 3; ++j) {
            triangles[i].lengths[j] = sqrt(pow(points[triangles[i].vertices[(j + 2) % 3]].x - points[triangles[i].vertices[(j + 1) % 3]].x, 2) + pow(points[triangles[i].vertices[(j + 2) % 3]].y - points[triangles[i].vertices[(j + 1) % 3]].y, 2));
        }
    }
    ifs >> s >> s >> wallsCount;
    walls.resize(wallsCount);
    for (size_t i = 0; i < wallsCount; ++i) {
        ifs >> walls[i].nr >> walls[i].vertices[0] >> walls[i].vertices[1] >> s;
        --walls[i].vertices[0];
        --walls[i].vertices[1];
        //return 1;
    }
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
