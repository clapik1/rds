#include "solver.h"

solver::solver(mesh &mMesh, vector2D &advection) {
    this->mMesh = mMesh;
    this->advection = advection;
    values.resize(mMesh.getPointsCount());
}

void solver::solveTriangle(const std::vector<vector2D> &norm, const std::vector<double> &localValues, methodRDS method, std::vector<double> &delta) const {
    double k[3], kp[3], km[3], beta[3], kpsum = 0, fi = 0;
    for(size_t i = 0; i < 3; ++i) {
        k[i] = dot(advection, norm[i]) / 2;
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
}

void solver::solve(double t, double dt, methodRDS method) {
    for (int it = 0; it < 200; ++it) {
        std::vector<double> nu(mMesh.getPointsCount()), si(mMesh.getPointsCount());
        for(size_t i = 0; i < mMesh.getTrianglesCount(); ++i) {
            std::vector<double> localValues(3), delta(3);
            std::vector<vector2D> norm(3);

            for(size_t j = 0; j < 3; ++j) {
                localValues[j] = values[mMesh.triangle(i).vertices[j]];
                double a1 = mMesh.point(mMesh.triangle(i).vertices[(j + 2) % 3]).x - mMesh.point(mMesh.triangle(i).vertices[(j + 1) % 3]).x;
                double a2 = mMesh.point(mMesh.triangle(i).vertices[(j + 2) % 3]).y - mMesh.point(mMesh.triangle(i).vertices[(j + 1) % 3]).y;
                double b1 = mMesh.point(mMesh.triangle(i).vertices[j]).x - mMesh.point(mMesh.triangle(i).vertices[(j + 1) % 3]).x;
                double b2 = mMesh.point(mMesh.triangle(i).vertices[j]).y - mMesh.point(mMesh.triangle(i).vertices[(j + 1) % 3]).y;
                norm[j].x = a2 * (a2 * b1 - a1 * b2);
                norm[j].y = a1 * (a1 * b2 - a2 * b1);
                norm[j].normalize();
                norm[j] *= mMesh.triangle(i).lengths[j];

                si[mMesh.triangle(i).vertices[j]] += mMesh.triangle(i).getArea() / 3;
            }

            solveTriangle(norm, localValues, method, delta);
            for(size_t j = 0; j < 3; ++j) {
                nu[mMesh.triangle(i).vertices[j]] += delta[j];
            }
        }
        for(size_t i = 0; i < mMesh.getPointsCount(); ++i) {
            values[i] -= dt * nu[i] / si[i];
        }
    }
}

void solver::toTecplot(std::ostream &os) const {
    os << "TITLE = \"\"\nVARIABLES = \"X\", \"Y\", \"VALUE\"\nZONE T=\"\", N=" << mMesh.getPointsCount() << ", E=" << mMesh.getTrianglesCount() << ", F=FEPOINT, ET=QUADRILATERAL\n";
    for(size_t i = 0; i < mMesh.getPointsCount(); ++i) {
        os << mMesh.point(i).x << ' ' << mMesh.point(i).y << ' ' << values[i] << '\n';
    }
    for(size_t i = 0; i < mMesh.getTrianglesCount(); ++i) {
        os << mMesh.triangle(i).vertices[0] + 1 << ' ' << mMesh.triangle(i).vertices[1] + 1 << ' ' << mMesh.triangle(i).vertices[2] + 1 << ' ' << mMesh.triangle(i).vertices[2] + 1 << '\n';
    }
}
