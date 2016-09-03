#include <set>
#include <iostream>
#include <cmath>
#include "constants.h"
#include "unstatSolverExplicit.h"

unstatSolverExplicit::unstatSolverExplicit(std::istream &meshStream, vector2D &advection, methodUnstat method) : mMesh(meshStream), advection(advection), method(method) {
    values.resize(mMesh.getPoints().size());

    for(size_t i = 0; i < values.size(); ++i) {
        if(std::pow(mMesh.getPoints()[i].x - (-0.25), 2) + std::pow(mMesh.getPoints()[i].y - (-0.25), 2) < 0.01)
            values[i] = 1.;
    } //fixme delete this shit
}

void unstatSolverExplicit::unstatSolve(double t, double (*wallElemValue)(double, double)) {
    std::vector<double> prevValues(values.size());
    for(size_t i = 0; i < values.size(); ++i) {
        prevValues[i] = values[i];
    }

    for(auto p = 0.; p < t; p += dt) {
        std::vector<double> nu(mMesh.getPoints().size()), si(mMesh.getPoints().size()), ti(mMesh.getPoints().size()), tem(mMesh.getPoints().size());

        #pragma omp parallel for
        for(size_t i = 0; i < mMesh.getTriangles().size(); ++i) {
            std::array<double, 3> localPrevValues, localValues, delta;
            std::array<point2D, 3> coords;

            for(size_t j = 0; j < 3; ++j) {
                coords[j] = mMesh.getPoints()[mMesh.getTriangles()[i].vertices[j]];
                localPrevValues[j] = prevValues[mMesh.getTriangles()[i].vertices[j]];
                localValues[j] = values[mMesh.getTriangles()[i].vertices[j]];
            }

            triangle2D tr;
            tr.updateArea(coords);

            std::array<double, 3> k = calcK(coords, advection);
            std::array<double, 3> kTilde = calcKTilde(tr.getArea(), k);
            std::array<double, 3> beta = calcUnstatBeta(tr.getArea(), k, method);

            double fi = 0.;
            for(size_t j = 0; j < 3; ++j) {
                fi += kTilde[j] * localValues[j] + (kTilde[j] - 2 * tr.getArea() / 3) * localPrevValues[j];
            }

            for(size_t j = 0; j < 3; ++j) {
                #pragma omp atomic
                nu[mMesh.getTriangles()[i].vertices[j]] += beta[j] * fi;
                #pragma omp atomic
                si[mMesh.getTriangles()[i].vertices[j]] += tr.getArea() / 3;
                #pragma omp atomic
                ti[mMesh.getTriangles()[i].vertices[j]] += dt * std::max(0., k[j]) / 2 + tr.getArea() / 3;
                #pragma omp atomic
                tem[mMesh.getTriangles()[i].vertices[j]] += kTilde[j];
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
            values[i] -= nu[i] / tem[i];
            //values[i] -= nu[i] * ti[i] / si[i];
            std::cout << 1. / tem[i] << ' ' << ti[i] / si[i] << std::endl;
        }
    }
}

void unstatSolverExplicit::toTecplot(std::ostream &os) const {
    os << "TITLE = \"\"\nVARIABLES = \"X\", \"Y\", \"VALUE\"\nZONE T=\"\", N=" << mMesh.getPoints().size() << ", E=" << mMesh.getTriangles().size() << ", F=FEPOINT, ET=QUADRILATERAL\n";
    for(size_t i = 0; i < mMesh.getPoints().size(); ++i) {
        os << mMesh.getPoints()[i].x << ' ' << mMesh.getPoints()[i].y << ' ' << values[i] << '\n';
    }
    for(size_t i = 0; i < mMesh.getTriangles().size(); ++i) {
        os << mMesh.getTriangles()[i].vertices[0] + 1 << ' ' << mMesh.getTriangles()[i].vertices[1] + 1 << ' ' << mMesh.getTriangles()[i].vertices[2] + 1 << ' ' << mMesh.getTriangles()[i].vertices[2] + 1 << '\n';
    }
}
