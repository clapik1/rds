#include <cmath>
#include <iostream>
#include <stdexcept>
#include "steadySolver.h"
#include "constants.h"

statSolver::statSolver(mesh &mMesh, vector2D &advection, steadyMethods method) : mMesh(mMesh), advection(advection), method(method) {}

std::array<double, 3> statSolver::solveTriangle(std::array<point2D, 3> &coords, std::array<double, 3> &localValues) const {
    std::array<double, 3> delta;
    std::array<double, 3> k = calcK(coords, advection);
    double N = calcN(k);
    double inflow = calcInflow(N, k, localValues);
    double outflow = calcOutflow(N, k, localValues);

    switch(method) {
        case steadyMethods::N:
            delta = distributeN(k, localValues, inflow);
            break;
        case steadyMethods::LDA:
            delta = distributeLDA(k, inflow, outflow);
            break;
        case steadyMethods::Blended: {
            std::array<double, 3> NDist = distributeN(k, localValues, inflow);
            std::array<double, 3> LDADist = distributeLDA(k, inflow, outflow);
            double fi = (outflow - inflow) / N;
            double theta = 0.;
            for (size_t i = 0; i < 3; ++i) {
                theta += std::abs(NDist[i]);
            }
            if (theta != 0.)
                theta = std::abs(fi) / theta;
            else
                theta = 1.;
            for (size_t i = 0; i < 3; ++i) {
                delta[i] = (1 - theta) * LDADist[i] + theta * NDist[i];
            }
            break;
        }
        case steadyMethods::LimitedN: {
            double fi = (outflow - inflow) / N;
            std::array<double, 3> NDist = distributeN(k, localValues, inflow);
            delta = applyMapping(fi, NDist);
            break;
        }
    }
    return delta;
}

void statSolver::solve(double (*wallElemValue)(double, double)) {
    double change, dt = 0.001; // fixme: not tested!
    do {
        std::vector<double> nu(mMesh.getPoints().size()), si(mMesh.getPoints().size());

        #pragma omp parallel for
        for(size_t i = 0; i < mMesh.getTriangles().size(); ++i) {
            std::array<double, 3> localValues, delta;
            std::array<point2D, 3> coords;

            for(size_t j = 0; j < 3; ++j) {
                coords[j] = mMesh.getPoints()[mMesh.getTriangles()[i].vertices[j]];
                localValues[j] = mMesh.getValues()[mMesh.getTriangles()[i].vertices[j]];
            }

            delta = solveTriangle(coords, localValues);
            for(size_t j = 0; j < 3; ++j) {
                #pragma omp atomic
                nu[mMesh.getTriangles()[i].vertices[j]] += delta[j];
                #pragma omp atomic
                si[mMesh.getTriangles()[i].vertices[j]] += mMesh.getTriangles()[i].getArea() / 3;
            }
        }

        #pragma omp parallel for
        for(size_t i = 0; i < mMesh.getWallElems().size(); ++i) {
            std::array<double, 3> localValues, delta;
            std::array<point2D, 3> coords;

            for(size_t j = 0; j < 2; ++j) {
                coords[j + 1] = mMesh.getPoints()[mMesh.getWallElems()[i].vertices[j]];
                localValues[j + 1] = mMesh.getValues()[mMesh.getWallElems()[i].vertices[j]];
            }
            vector2D side(coords[2].x - coords[1].x, coords[2].y - coords[1].y);
            double len = side.length();
            for(size_t j = 0; j < 2; ++j) {
                coords[0] = coords[j + 1];
                localValues[0] = wallElemValue(coords[0].x, coords[0].y);
                coords[0].x += ghostHeight * len * ghost_dx[mMesh.getWallElems()[i].wallNr];
                coords[0].y += ghostHeight * len * ghost_dy[mMesh.getWallElems()[i].wallNr];

                delta = solveTriangle(coords, localValues);

                triangle2D tr;
                tr.updateArea(coords);

                #pragma omp atomic
                nu[mMesh.getWallElems()[i].vertices[j]] += delta[j + 1];
                #pragma omp atomic
                si[mMesh.getWallElems()[i].vertices[j]] += tr.getArea() / 3;
            }
        }

        change = 0.;
        for(size_t i = 0; i < mMesh.getPoints().size(); ++i) {
            change += std::abs(dt * nu[i] / si[i]);
            mMesh.addToValue(i, -dt * nu[i] / si[i]);
        }
        change /= mMesh.getPoints().size();
        std::cout << change << std::endl;
    }
    while(change > 1e-10);
}

double statSolver::calcError(double (*wallElemValue)(double, double)) {
    double x, y, x1, y1, squares = 0.;
    if(advection.x >= 0. && advection.y >= 0.) {
        for (size_t i = 0; i < mMesh.getPoints().size(); ++i) {
            x = mMesh.getPoints()[i].x;
            y = mMesh.getPoints()[i].y;
            if(x * advection.y > y * advection.x) {
                y1 = -0.5;
                x1 = x + (y1 - y) * advection.x / advection.y;
            }
            else {
                x1 = -0.5;
                y1 = y + (x1 - x) * advection.y / advection.x;
            }
            squares += std::pow(std::abs(wallElemValue(x1, y1) - mMesh.getValues()[i]), 2);
        }
        return std::sqrt(squares / mMesh.getPoints().size());
    }
    throw std::invalid_argument("negative advection");
}
