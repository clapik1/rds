#include <set>
#include <iostream>
#include <cmath>
#include <iomanip>
#include "../../constants.h"
#include "unstatSolverImplicit.h"
#include "functions/statFunctions.h"

unstatSolverImplicit::unstatSolverImplicit(mesh &mesh, vector2D &advection, methodUnstat method, double dt)
        : mMesh(mesh), advection(advection), method(method), dt(dt) {}

std::array<double, 3> unstatSolverImplicit::unstatSolveTriangle(std::array<point2D, 3> &coords, std::array<double, 3> &localValues, std::array<double, 3> &localPrevValues) const {
    std::array<double, 3> delta;

    triangle2D tr;
    tr.updateArea(coords);

    std::array<double, 3> k = calcK(coords, advection);
    std::array<double, 3> kTilde = calcKTilde(dt, tr.getArea(), k);
    std::array<double, 3> kHat = calcKHat(dt, tr.getArea(), k);
    double N = calcN(k);
    double NST = calcNST(kTilde, kHat);
    double inflowST = calcInflowST(NST, kTilde, kHat, localValues, localPrevValues);
    double outflowST = calcOutflowST(NST, kTilde, kHat, localValues, localPrevValues);

    double fi = (outflowST - inflowST) / NST;

    switch(method) {
        case methodUnstat::N:break;
        case methodUnstat::N_ST:break;
        case methodUnstat::LDA:
            for(size_t i = 0; i < 3; ++i) {
                delta[i] = std::max(0., k[i]) * N * fi;
            }
            break;
        case methodUnstat::LDA_ST:
            for(size_t i = 0; i < 3; ++i) {
                delta[i] = std::max(0., kTilde[i]) * NST * fi;
            }
            break;
    }

    return delta;
}

void unstatSolverImplicit::unstatSolve(double t, double (*wallElemValue)(double, double)) {
    std::vector<double> prevValues(mMesh.getValues());

    for(auto p = 0.; p < t; p += dt) {
        double change = 1.;
        int it = 0;

        while(std::abs(change) > 1e-6) {
            std::vector<double> nu(mMesh.getPoints().size()), tem(mMesh.getPoints().size());

            #pragma omp parallel for
            for (size_t i = 0; i < mMesh.getTriangles().size(); ++i) {
                std::array<double, 3> localPrevValues, localValues, delta;
                std::array<point2D, 3> coords;

                for (size_t j = 0; j < 3; ++j) {
                    coords[j] = mMesh.getPoints()[mMesh.getTriangles()[i].vertices[j]];
                    localPrevValues[j] = prevValues[mMesh.getTriangles()[i].vertices[j]];
                    localValues[j] = mMesh.getValues()[mMesh.getTriangles()[i].vertices[j]];
                }

                triangle2D tr;
                tr.updateArea(coords);

                std::array<double, 3> k = calcK(coords, advection);
                std::array<double, 3> kTilde = calcKTilde(dt, tr.getArea(), k);

                delta = unstatSolveTriangle(coords, localValues, localPrevValues);

                for (size_t j = 0; j < 3; ++j) {
                    #pragma omp atomic
                    nu[mMesh.getTriangles()[i].vertices[j]] += delta[j];
                    #pragma omp atomic
                    tem[mMesh.getTriangles()[i].vertices[j]] += kTilde[j];
                }
            }

            #pragma omp parallel for
            for(size_t i = 0; i < mMesh.getWallElems().size(); ++i) {
                std::array<double, 3> localPrevValues, localValues, delta;
                std::array<point2D, 3> coords;

                for(size_t j = 0; j < 2; ++j) {
                    coords[j + 1] = mMesh.getPoints()[mMesh.getWallElems()[i].vertices[j]];
                    localValues[j + 1] = mMesh.getValues()[mMesh.getWallElems()[i].vertices[j]];
                    localPrevValues[j + 1] = prevValues[mMesh.getWallElems()[i].vertices[j]];
                }
                vector2D side(coords[2].x - coords[1].x, coords[2].y - coords[1].y);
                double len = side.length();
                for(size_t j = 0; j < 2; ++j) {
                    coords[0] = coords[j + 1];
                    localValues[0] = wallElemValue(coords[0].x, coords[0].y);
                    localPrevValues[0] = localValues[0];
                    coords[0].x += ghostHeight * len * ghost_dx[mMesh.getWallElems()[i].wallNr];
                    coords[0].y += ghostHeight * len * ghost_dy[mMesh.getWallElems()[i].wallNr];

                    delta = unstatSolveTriangle(coords, localValues, localPrevValues);

                    #pragma omp atomic
                    nu[mMesh.getWallElems()[i].vertices[j]] += delta[j + 1];
                }
            }

            change = 0;
            for (size_t i = 0; i < mMesh.getPoints().size(); ++i) {
                change += std::abs(nu[i] / tem[i]);
                mMesh.addToValue(i, -nu[i] / tem[i]);
            }
            change /= mMesh.getPoints().size();

            //std::cout << it << ' ' << change << std::endl;
            ++it;
        }

        prevValues = mMesh.getValues();
        std::cout << "t = " << std::fixed << std::setprecision(6) << p << "\tit = " << it << std::endl;
    }
}
