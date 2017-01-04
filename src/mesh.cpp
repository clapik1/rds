#include "mesh.h"
#include <fstream>
#include <iostream>
#include <cmath>

mesh::mesh(std::istream &ifs) {
    std::string s;
    size_t pointsCount, trianglesCount, wallElemsCount;

    ifs >> s >> s >> s >> pointsCount;
    points.resize(pointsCount);
    values.resize(pointsCount);
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

const std::vector<double>& mesh::getValues() const {
    return values;
}

void mesh::setValue(size_t i, double value) {
    values[i] = value;
}

void mesh::addToValue(size_t i, double value) {
    values[i] += value;
}

void mesh::toTecplot(std::ostream &os) const {
    os << "TITLE = \"Advection Solution\"\nVARIABLES = \"X\", \"Y\", \"VALUE\"\nZONE N = " << getPoints().size() << ", E = " << getTriangles().size() << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE\n";
    for(size_t i = 0; i < getPoints().size(); ++i) {
        os << getPoints()[i].x << ' ' << getPoints()[i].y << ' ' << values[i] << '\n';
    }
    for(size_t i = 0; i < getTriangles().size(); ++i) {
        os << getTriangles()[i].vertices[0] + 1 << ' ' << getTriangles()[i].vertices[1] + 1 << ' ' << getTriangles()[i].vertices[2] + 1 << '\n';
    }
}

void mesh::toTecplotAnimationHeaderAndFirstZone(std::ostream &os) const {
    os << "TITLE = \"Advection Solution\"\nVARIABLES = \"X\", \"Y\", \"VALUE\"\nZONE STRANDID = 1, SOLUTIONTIME = 0, N = " << getPoints().size() << ", E = " << getTriangles().size() << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE\n";
    for(size_t i = 0; i < getPoints().size(); ++i) {
        os << getPoints()[i].x << ' ' << getPoints()[i].y << ' ' << values[i] << '\n';
    }
    for(size_t i = 0; i < getTriangles().size(); ++i) {
        os << getTriangles()[i].vertices[0] + 1 << ' ' << getTriangles()[i].vertices[1] + 1 << ' ' << getTriangles()[i].vertices[2] + 1 << '\n';
    }
}

void mesh::toTecplotAnimationNextZone(std::ostream &os, double t) const {
    os << "ZONE STRANDID = 1, SOLUTIONTIME = " << t << ", N = " << getPoints().size() << ", E = " << getTriangles().size() << ", DATAPACKING = POINT, ZONETYPE = FETRIANGLE, VARSHARELIST = ([1-2]=1), CONNECTIVITYSHAREZONE = 1\n";
    for(size_t i = 0; i < getPoints().size(); ++i) {
        os << values[i] << '\n';
    }
}
