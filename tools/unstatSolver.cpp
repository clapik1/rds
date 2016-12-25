#include "unstatSolver.h"
#include <set>
#include <iostream>
#include <cmath>
#include "../src/constants.h"

unstatSolver::unstatSolver(std::istream &meshStream, vector2D &advection, methodUnstat method) : mMesh(meshStream), advection(advection), method(method) {
    values.resize(mMesh.getPoints().size());

    matrix.resize(values.size());
    for(size_t i = 0; i < values.size(); ++i) {
        matrix[i].resize(values.size());
    }

    std::vector<std::set<size_t>> adjList(values.size());
    for(size_t i = 0; i < mMesh.getTriangles().size(); ++i) {
        for(size_t j = 0; j < 3; ++j) {
            adjList[mMesh.getTriangles()[i].vertices[j]].insert(mMesh.getTriangles()[i].vertices[j]);
            adjList[mMesh.getTriangles()[i].vertices[j]].insert(mMesh.getTriangles()[i].vertices[(j + 1) % 3]);
            adjList[mMesh.getTriangles()[i].vertices[j]].insert(mMesh.getTriangles()[i].vertices[(j + 2) % 3]);
        }
    }

    for(size_t i = 0; i < mMesh.getTriangles().size(); ++i) {
        std::array<point2D, 3> coords;
        for(size_t j = 0; j < 3; ++j) {
            coords[j] = mMesh.getPoints()[mMesh.getTriangles()[i].vertices[j]];
        }
        triangle2D tr;
        tr.updateArea(coords);

        std::array<double, 3> k = calcK(coords, advection);
        std::array<double, 3> kTilde = calcKTilde(tr.getArea(), k);
        std::array<double, 3> beta = calcUnstatBeta(tr.getArea(), k, method);

        //std::cout << i << ' ';
        for(size_t j = 0; j < 3; ++j) {
            std::cout << k[j] << ' ' << kTilde[j] << ' ' << beta[j] << ' ' << std::endl;
            for(auto&& m : adjList[mMesh.getTriangles()[i].vertices[j]]) {
                std::cout << m << ' ' << mMesh.getTriangles()[i].vertices[j] << ' ' << beta[j] * kTilde[j] << std::endl;
                matrix[m][mMesh.getTriangles()[i].vertices[j]] += beta[j] * kTilde[j];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    for(size_t i = 0; i < values.size(); ++i) {
        for(size_t j = 0; j < values.size(); ++j) {
            std::cout << matrix[i][j] << ' ';
        }
        std::cout << std::endl;
    }

    for(size_t i = 0; i < values.size(); ++i) {
        if(std::pow(mMesh.getPoints()[i].x - (-0.25), 2) + std::pow(mMesh.getPoints()[i].y - (-0.25), 2) < 0.01)
            values[i] = 1.;
    } //fixme delete this shit
}

void unstatSolver::unstatSolve(double t, double (*wallElemValue)(double, double)) {
    std::vector<double> prevValues(values.size());
    for(size_t i = 0; i < values.size(); ++i) {
        prevValues[i] = values[i];
    }

    for(auto p = 0.; p < t; p += dt) {

#pragma omp parallel for
        for(size_t i = 0; i < mMesh.getTriangles().size(); ++i) {
            std::array<double, 3> localPrevValues, localValues, delta;
            std::array<point2D, 3> coords;

            for(size_t j = 0; j < 3; ++j) {
                coords[j] = mMesh.getPoints()[mMesh.getTriangles()[i].vertices[j]];
                localPrevValues[j] = prevValues[mMesh.getTriangles()[i].vertices[j]];
                localValues[j] = values[mMesh.getTriangles()[i].vertices[j]];
            }


        }
        /*
        #pragma omp parallel for
        for(size_t i = 0; i < mMesh.getWallElems().size(); ++i) {
            std::array<double, 3> localValues, delta;
            std::array<point2D, 3> coords;

            for(size_t j = 0; j < 2; ++j) {
                coords[j + 1] = mMesh.getPoints()[mMesh.getWallElems()[i].vertices[j]];
                localValues[j + 1] = values[mMesh.getWallElems()[i].vertices[j]];
            }
            for(size_t j = 0; j < 2; ++j) {
                coords[0] = coords[j + 1];
                localValues[0] = wallElemValue(coords[0].x, coords[0].y);
                coords[0].x += ghostHeight * ghost_dx[mMesh.getWallElems()[i].wallNr];
                coords[0].y += ghostHeight * ghost_dy[mMesh.getWallElems()[i].wallNr];

                delta = statSolveTriangle(coords, localValues, method);

                #pragma omp atomic
                nu[mMesh.getWallElems()[i].vertices[j]] += delta[j + 1];
            }
        }*/

        for(size_t i = 0; i < mMesh.getPoints().size(); ++i) {
            prevValues[i] = values[i];
            values[i] -= 0;
        }
    }
}

void unstatSolver::toTecplot(std::ostream &os) const {
    os << "TITLE = \"\"\nVARIABLES = \"X\", \"Y\", \"VALUE\"\nZONE T=\"\", N=" << mMesh.getPoints().size() << ", E=" << mMesh.getTriangles().size() << ", F=FEPOINT, ET=QUADRILATERAL\n";
    for(size_t i = 0; i < mMesh.getPoints().size(); ++i) {
        os << mMesh.getPoints()[i].x << ' ' << mMesh.getPoints()[i].y << ' ' << values[i] << '\n';
    }
    for(size_t i = 0; i < mMesh.getTriangles().size(); ++i) {
        os << mMesh.getTriangles()[i].vertices[0] + 1 << ' ' << mMesh.getTriangles()[i].vertices[1] + 1 << ' ' << mMesh.getTriangles()[i].vertices[2] + 1 << ' ' << mMesh.getTriangles()[i].vertices[2] + 1 << '\n';
    }
}
