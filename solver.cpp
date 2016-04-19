#include "solver.h"

solver::solver(mesh &mMesh, vector2D &advection) {
    this->mMesh = mMesh;
    this->advection = advection;
    values.resize(mMesh.getPointsCount());
}

void solver::solve(double t, double dt, methodRDS method) {
    for (int it = 0; it < 500; ++it) {
        std::vector<double> nu(mMesh.getPointsCount()), si(mMesh.getPointsCount());
        for(size_t i = 0; i < mMesh.getTrianglesCount(); ++i) {
            double k[3], kp[3], km[3], beta[3], kpsum = 0, fi = 0;
            for(size_t j = 0; j < 3; ++j) {
                k[j] = dot(advection, mMesh.triangle(i).norm[j]) / 2;
                kp[j] = std::max(0., k[j]);
                km[j] = std::min(0., k[j]);
                kpsum += kp[j];
                fi += k[j] * values[mMesh.triangle(i).vertices[j]];
                si[mMesh.triangle(i).vertices[j]] += mMesh.triangle(i).getArea() / 3;
            }
            for(size_t j = 0; j < 3; ++j) {
                if(method == methodRDS::N) {
                    double s = 0;
                    for(size_t m = 0; m < 3; ++m) {
                        s += km[m] * values[mMesh.triangle(i).vertices[m]] / kpsum;
                    }
                    nu[mMesh.triangle(i).vertices[j]] += kp[j] * (values[mMesh.triangle(i).vertices[j]] + s);
                }
                else if(method == methodRDS::LDA) {
                    beta[j] = kp[j] / kpsum;
                    nu[mMesh.triangle(i).vertices[j]] += beta[j] * fi;
                }
            }
        }
        for(size_t i = 0; i < mMesh.getPointsCount(); ++i) {
            values[i] -= dt * nu[i] / si[i]; // (very) temporary fix
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
