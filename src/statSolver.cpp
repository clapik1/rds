#include <cmath>
#include <iostream>
#include "statSolver.h"
#include "constants.h"
#include "functions.h"

statSolver::statSolver(std::istream &meshStream, vector2D &advection, methodStat method) : mMesh(meshStream), advection(advection), method(method) {
    values.resize(mMesh.getPoints().size());
}

std::array<double, 3> statSolver::statSolveTriangle(std::array<point2D, 3> &coords, std::array<double, 3> &localValues) const {
    std::array<double, 3> delta;

    std::vector<vector2D> norm(3);
    for(size_t i = 0; i < 3; ++i) {
        double ax = coords[(i + 2) % 3].x;
        double ay = coords[(i + 2) % 3].y;
        double bx = coords[(i + 1) % 3].x;
        double by = coords[(i + 1) % 3].y;
        norm[i].x = by - ay;
        norm[i].y = ax - bx;
    }

    double k[3], kp[3], km[3], kpsum = 0, fi = 0, theta = 0;
    for(size_t i = 0; i < 3; ++i) {
        k[i] = dotProduct(advection, norm[i]) / 2;
        kp[i] = std::max(0., k[i]);
        km[i] = std::min(0., k[i]);
        kpsum += kp[i];
        fi += k[i] * localValues[i];
    }
    if(fi != 0.) {
        if (method == methodStat::N) {
            double s = 0;
            for (size_t i = 0; i < 3; ++i) {
                s += km[i] * localValues[i] / kpsum;
            }
            for (size_t i = 0; i < 3; ++i) {
                delta[i] = kp[i] * (localValues[i] + s);
            }
        }
        else if (method == methodStat::LDA) {
            for (size_t i = 0; i < 3; ++i) {
                delta[i] = (kp[i] / kpsum) * fi;
            }
        }
        else if (method == methodStat::Blended) {
            double s = 0;
            for (size_t i = 0; i < 3; ++i) {
                s += km[i] * localValues[i] / kpsum;
            }
            for (size_t i = 0; i < 3; ++i) {
                theta += std::abs(kp[i] * (localValues[i] + s));
            }
            if (theta != 0.)
                theta = std::abs(fi) / theta;
            else
                theta = 1.;
            for (size_t i = 0; i < 3; ++i) {
                delta[i] = (1 - theta) * ((kp[i] / kpsum) * fi) + theta * (kp[i] * (localValues[i] + s)); //fixme: fix this shit
            }
        }
        else if (method == methodStat::LimitedN) {
            double s = 0, beta[3];
            for (size_t i = 0; i < 3; ++i) {
                s += km[i] * localValues[i] / kpsum;
            }
            for (size_t i = 0; i < 3; ++i) {
                beta[i] = kp[i] * (localValues[i] + s) / fi;
            }
            s = 0.;
            for (size_t i = 0; i < 3; ++i) {
                s += std::max(0., beta[i]);
            }
            for (size_t i = 0; i < 3; ++i) {
                delta[i] = (std::max(0., beta[i]) + 1e-10) * fi / (s + 3e-10);
            }
        }
    }
    else {
        for (size_t i = 0; i < 3; ++i) {
            delta[i] = 0.;
        }
    }
    return delta;
}

void statSolver::statSolve(double (*wallElemValue)(double, double)) {
    double change, dt = 0.001; // fixme: not tested!
    do {
        std::vector<double> nu(mMesh.getPoints().size()), si(mMesh.getPoints().size());

        #pragma omp parallel for
        for(size_t i = 0; i < mMesh.getTriangles().size(); ++i) {
            std::array<double, 3> localValues, delta;
            std::array<point2D, 3> coords;

            for(size_t j = 0; j < 3; ++j) {
                coords[j] = mMesh.getPoints()[mMesh.getTriangles()[i].vertices[j]];
                localValues[j] = values[mMesh.getTriangles()[i].vertices[j]];
                #pragma omp atomic
                si[mMesh.getTriangles()[i].vertices[j]] += mMesh.getTriangles()[i].getArea() / 3;
            }

            delta = statSolveTriangle(coords, localValues);
            for(size_t j = 0; j < 3; ++j) {
                #pragma omp atomic
                nu[mMesh.getTriangles()[i].vertices[j]] += delta[j];
            }
        }

        #pragma omp parallel for
        for(size_t i = 0; i < mMesh.getWallElems().size(); ++i) {
            std::array<double, 3> localValues, delta;
            std::array<point2D, 3> coords;

            for(size_t j = 0; j < 2; ++j) {
                coords[j + 1] = mMesh.getPoints()[mMesh.getWallElems()[i].vertices[j]];
                localValues[j + 1] = values[mMesh.getWallElems()[i].vertices[j]];
            }
            vector2D side(coords[2].x - coords[1].x, coords[2].y - coords[1].y);
            double len = side.length();
            for(size_t j = 0; j < 2; ++j) {
                coords[0] = coords[j + 1];
                localValues[0] = wallElemValue(coords[0].x, coords[0].y);
                coords[0].x += ghostHeight * len * ghost_dx[mMesh.getWallElems()[i].wallNr];
                coords[0].y += ghostHeight * len * ghost_dy[mMesh.getWallElems()[i].wallNr];

                delta = statSolveTriangle(coords, localValues);

                #pragma omp atomic
                nu[mMesh.getWallElems()[i].vertices[j]] += delta[j + 1];
            }
        }

        change = 0.;
        for(size_t i = 0; i < mMesh.getPoints().size(); ++i) {
            change += std::abs(dt * nu[i] / si[i]);
            values[i] -= dt * nu[i] / si[i];
        }
        change /= mMesh.getPoints().size();
        std::cout << change << std::endl;
    }
    while(change > 1e-10);
}

double statSolver::statCheck(double (*wallElemValue)(double, double)) {
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
            squares += std::pow(std::abs(wallElemValue(x1, y1) - values[i]), 2);
        }
        return std::sqrt(squares / mMesh.getPoints().size());
    }
    return 5;
}

void statSolver::toTecplot(std::ostream &os) const {
    os << "TITLE = \"\"\nVARIABLES = \"X\", \"Y\", \"VALUE\"\nZONE T=\"\", N=" << mMesh.getPoints().size() << ", E=" << mMesh.getTriangles().size() << ", F=FEPOINT, ET=QUADRILATERAL\n";
    for(size_t i = 0; i < mMesh.getPoints().size(); ++i) {
        os << mMesh.getPoints()[i].x << ' ' << mMesh.getPoints()[i].y << ' ' << values[i] << '\n';
    }
    for(size_t i = 0; i < mMesh.getTriangles().size(); ++i) {
        os << mMesh.getTriangles()[i].vertices[0] + 1 << ' ' << mMesh.getTriangles()[i].vertices[1] + 1 << ' ' << mMesh.getTriangles()[i].vertices[2] + 1 << ' ' << mMesh.getTriangles()[i].vertices[2] + 1 << '\n';
    }
}
