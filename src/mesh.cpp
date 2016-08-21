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
        std::array<point2D, 3> coords;
        for(size_t j = 0; j < 3; ++j) {
            ifs >> triangles[i].vertices[j];
            --triangles[i].vertices[j];
            coords[j] = points[triangles[i].vertices[j]];
        }
        triangles[i].updateArea(coords);
    }
    ifs >> s >> s >> wallElemsCount;
    wallElems.resize(wallElemsCount);
    for (size_t i = 0; i < wallElemsCount; ++i) {
        ifs >> wallElems[i].wallNr >> wallElems[i].vertices[1] >> wallElems[i].vertices[0] >> s;
        --wallElems[i].wallNr;
        --wallElems[i].vertices[1];
        --wallElems[i].vertices[0];
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