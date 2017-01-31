#include <set>
#include <iostream>
#include <cmath>
#include <iomanip>
#include "../../constants.h"
#include "timeDepSolverExplicit.h"
#include "functions/steadyFunctions.h"

timeDepSolverExplicit::timeDepSolverExplicit(mesh &mMesh, vector2D &advection, timeDepMethods method)
        : mMesh(mMesh), advection(advection), method(method) {}

std::array<double, 3> timeDepSolverExplicit::solveTriangle(std::array<point2D, 3> &coords, std::array<double, 3> &localValues,
                                                           std::array<double, 3> &localPrevValues, double dt) const {
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
        case timeDepMethods::N: {
            double inflow = calcInflow(N, k, localValues);
            double prevInflow = calcInflow(N, k, localPrevValues);
            for(size_t i = 0; i < 3; ++i) {
                delta[i] = (localValues[i] - localPrevValues[i]) * tr.getArea() / 3 + (localValues[i] + localPrevValues[i] - inflow - prevInflow) * std::max(0., k[i]) * dt / 2;
            }
            break;
        }
        case timeDepMethods::N_ST:
            for(size_t i = 0; i < 3; ++i) {
                delta[i] = std::max(0., kTilde[i]) * (localValues[i] - inflowST);
            }
            break;
        case timeDepMethods::LDA:
            for(size_t i = 0; i < 3; ++i) {
                delta[i] = std::max(0., k[i]) * N * fi;
            }
            break;
        case timeDepMethods::LDA_ST:
            for(size_t i = 0; i < 3; ++i) {
                delta[i] = std::max(0., kTilde[i]) * NST * fi;
            }
            break;
        case timeDepMethods::LN: {
            double inflow = calcInflow(N, k, localValues);
            double prevInflow = calcInflow(N, k, localPrevValues);
            for(size_t i = 0; i < 3; ++i) {
                delta[i] = (localValues[i] - localPrevValues[i]) * tr.getArea() / 3 + (localValues[i] + localPrevValues[i] - inflow - prevInflow) * std::max(0., k[i]) * dt / 2;
            }
            delta = applyMapping(fi, delta);
            break;
        }
        case timeDepMethods::LN_ST:
            for(size_t i = 0; i < 3; ++i) {
                delta[i] = std::max(0., kTilde[i]) * (localValues[i] - inflowST);
            }
            delta = applyMapping(fi, delta);
            break;
    }

    return delta;
}

void
timeDepSolverExplicit::timeStep(std::vector<double> &prevValues, double (*wallElemValue)(double, double), double dt) {
    double convergence = 1.;
    int iter_count = 0;

    while(convergence > 1e-11) {
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

            delta = solveTriangle(coords, localValues, localPrevValues, dt);

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

                triangle2D tr;
                tr.updateArea(coords);

                std::array<double, 3> k = calcK(coords, advection);
                std::array<double, 3> kTilde = calcKTilde(dt, tr.getArea(), k);

                delta = solveTriangle(coords, localValues, localPrevValues, dt);

                #pragma omp atomic
                nu[mMesh.getWallElems()[i].vertices[j]] += delta[j + 1];
                #pragma omp atomic
                tem[mMesh.getWallElems()[i].vertices[j]] += kTilde[j + 1];
            }
        }

        convergence = 0;
        for (size_t i = 0; i < mMesh.getPoints().size(); ++i) {
            convergence += std::abs(nu[i]);
            mMesh.addToValue(i, -0.9 * nu[i] / tem[i]);
        }
        convergence /= mMesh.getPoints().size();

        //std::cout << iter_count << ' ' << convergence << std::endl;
        ++iter_count;
    }

    prevValues = mMesh.getValues();
    //std::cout << "\titer_count = " << std::fixed << std::setprecision(6) << iter_count << '\n';
}

void timeDepSolverExplicit::solve(double t, double (*wallElemValue)(double, double), unsigned int iterTotal) {
    std::vector<double> prevValues(mMesh.getValues());
    auto dt = t / iterTotal;
    for(auto i = 0; i < iterTotal; ++i) {
        //std::cout << "t = " << std::fixed << std::setprecision(6) << t * (i + 1) / iterTotal;
        timeStep(prevValues, wallElemValue, dt);
    }
}

void timeDepSolverExplicit::animate(double t, double (*wallElemValue)(double, double), unsigned int iterTotal, std::ostream &solutionStream) {
    std::vector<double> prevValues(mMesh.getValues());

    mMesh.toTecplotAnimationHeaderAndFirstZone(solutionStream);
    auto dt = t / iterTotal;
    for(auto i = 0; i < iterTotal; ++i) {
        //std::cout << "t = " << std::fixed << std::setprecision(6) << t * (i + 1) / iterTotal;
        timeStep(prevValues, wallElemValue, dt);
        mMesh.toTecplotAnimationNextZone(solutionStream, t * (i + 1) / iterTotal);
    }
}

double timeDepSolverExplicit::calcError(double t, double (*initialValue)(double, double), double (*wallElemValue)(double, double)) {
    double x, y, x1, y1, squares = 0.;
    for (size_t i = 0; i < mMesh.getPoints().size(); ++i) {
        x = mMesh.getPoints()[i].x;
        y = mMesh.getPoints()[i].y;
        x1 = x - t * advection.x;
        y1 = y - t * advection.y;
        if(-0.5 <= x1 && x1 <= 0.5 && -0.5 <= y1 && y1 <= 0.5) {
            squares += std::pow(initialValue(x1, y1) - mMesh.getValues()[i], 2);
        }
        else {
            if (x * advection.y > y * advection.x) {
                y1 = -0.5;
                x1 = x + (y1 - y) * advection.x / advection.y;
            } else {
                x1 = -0.5;
                y1 = y + (x1 - x) * advection.y / advection.x;
            }
            squares += std::pow(wallElemValue(x1, y1) - mMesh.getValues()[i], 2);
        }
    }
    return std::sqrt(squares / mMesh.getPoints().size());
}
