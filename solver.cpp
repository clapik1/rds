#include "solver.h"

solver::solver(std::istream &meshStream, vector2D &advection) : mMesh(meshStream), advection(advection) {
    values.resize(mMesh.getPoints().size());
}

std::array<double, 3> solver::solveTriangle(std::array<point2D, 3>& coords, std::array<double, 3>& localValues, methodRDS method) const {
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

    double k[3], kp[3], km[3], beta[3], kpsum = 0, fi = 0;
    for(size_t i = 0; i < 3; ++i) {
        k[i] = dotProduct(advection, norm[i]) / 2;
        kp[i] = std::max(0., k[i]);
        km[i] = std::min(0., k[i]);
        kpsum += kp[i];
        fi += k[i] * localValues[i];
    }
    if(method == methodRDS::N) {
        double s = 0;
        for(size_t i = 0; i < 3; ++i) {
            s += km[i] * localValues[i] / kpsum;
        }
        for(size_t i = 0; i < 3; ++i) {
            delta[i] = kp[i] * (localValues[i] + s);
        }
    }
    else if(method == methodRDS::LDA) {
        for(size_t i = 0; i < 3; ++i) {
            beta[i] = kp[i] / kpsum;
            delta[i] = beta[i] * fi;
        }
    }
    return delta;
}

void solver::solve(double t, double dt, double ghostHeight, methodRDS method) {
    for (auto it = 0; it < 200; ++it) {
        std::vector<double> nu(mMesh.getPoints().size()), si(mMesh.getPoints().size());
        for(size_t i = 0; i < mMesh.getTriangles().size(); ++i) {
            std::array<double, 3> localValues, delta;
            std::array<point2D, 3> coords;

            for(size_t j = 0; j < 3; ++j) {
                localValues[j] = values[mMesh.getTriangles()[i].vertices[j]];
                coords[j] = point2D(mMesh.getPoints()[mMesh.getTriangles()[i].vertices[j]]);
                si[mMesh.getTriangles()[i].vertices[j]] += mMesh.getTriangles()[i].getArea() / 3;
            }

            delta = solveTriangle(coords, localValues, method);
            for(size_t j = 0; j < 3; ++j) {
                nu[mMesh.getTriangles()[i].vertices[j]] += delta[j];
            }
        }
        for(size_t i = 0; i < mMesh.getWallElems().size(); ++i) {

        }

        for(size_t i = 0; i < mMesh.getPoints().size(); ++i) {
            values[i] -= dt * nu[i] / si[i];
        }
    }
}

void solver::toTecplot(std::ostream &os) const {
    os << "TITLE = \"\"\nVARIABLES = \"X\", \"Y\", \"VALUE\"\nZONE T=\"\", N=" << mMesh.getPoints().size() << ", E=" << mMesh.getTriangles().size() << ", F=FEPOINT, ET=QUADRILATERAL\n";
    for(size_t i = 0; i < mMesh.getPoints().size(); ++i) {
        os << mMesh.getPoints()[i].x << ' ' << mMesh.getPoints()[i].y << ' ' << values[i] << '\n';
    }
    for(size_t i = 0; i < mMesh.getTriangles().size(); ++i) {
        os << mMesh.getTriangles()[i].vertices[0] + 1 << ' ' << mMesh.getTriangles()[i].vertices[1] + 1 << ' ' << mMesh.getTriangles()[i].vertices[2] + 1 << ' ' << mMesh.getTriangles()[i].vertices[2] + 1 << '\n';
    }
}
