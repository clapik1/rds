#include "mesh.h"
#include <fstream>
#include <iostream>
#include <cmath>

mesh::mesh(std::istream &ifs) {
    std::string s;
    size_t pointsCount, trianglesCount, wallElemsCount;

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
    ifs >> s >> s >> wallElemsCount;
    wallElems.resize(wallElemsCount);
    for (size_t i = 0; i < wallElemsCount; ++i) {
        ifs >> wallElems[i].nr >> wallElems[i].vertices[0] >> wallElems[i].vertices[1] >> s;
        --wallElems[i].nr;
        --wallElems[i].vertices[0];
        --wallElems[i].vertices[1];
    }
}

const std::vector<triangle2D>& mesh::getTriangles() const {
    return triangles;
}

const std::vector<point2D>& mesh::getPoints() const {
    return points;
}

const std::vector<wall2D>& mesh::getWallElems() const {
    return wallElems;
}