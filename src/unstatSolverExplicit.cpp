#include <set>
#include <iostream>
#include <cmath>
#include <iomanip>
#include "constants.h"
#include "unstatSolverExplicit.h"

unstatSolverExplicit::unstatSolverExplicit(std::istream &meshStream, vector2D &advection, methodUnstat method) : mMesh(meshStream), advection(advection), method(method) {
    values.resize(mMesh.getPoints().size());
    double r;
    for(size_t i = 0; i < values.size(); ++i) {
        r = std::pow(mMesh.getPoints()[i].x - (-0.25), 2.) + std::pow(mMesh.getPoints()[i].y - (-0.25), 2.);
        if(r < 0.01)
            values[i] = std::cos(2 * std::acos(-1.) * r * 25) / 25;
    } //fixme delete this shit

    //values[8] = 1.;
}

void unstatSolverExplicit::unstatSolve(double t, double (*wallElemValue)(double, double)) {
    std::vector<double> prevValues(values.size());
    for(size_t i = 0; i < values.size(); ++i) {
        prevValues[i] = values[i];
    }

    for(auto p = 0.; p < t; p += dt) {
        double ddt = 100, ddti, ddtj, change = 1.;
        std::vector<double> nu(mMesh.getPoints().size()), si(mMesh.getPoints().size()), ti(mMesh.getPoints().size()), tem(mMesh.getPoints().size());

        int it = 0;
        while(std::abs(change) > 1e-6) {
            #pragma omp parallel for
            for (size_t i = 0; i < mMesh.getTriangles().size(); ++i) {
                std::array<double, 3> localPrevValues, localValues, delta;
                std::array<point2D, 3> coords;

                for (size_t j = 0; j < 3; ++j) {
                    coords[j] = mMesh.getPoints()[mMesh.getTriangles()[i].vertices[j]];
                    localPrevValues[j] = prevValues[mMesh.getTriangles()[i].vertices[j]];
                    localValues[j] = values[mMesh.getTriangles()[i].vertices[j]];
                }

                triangle2D tr;
                tr.updateArea(coords);

                std::array<double, 3> k = calcK(coords, advection);
                std::array<double, 3> kTilde = calcKTilde(tr.getArea(), k);
                std::array<double, 3> beta = calcUnstatBeta(tr.getArea(), k, method);

                /*double ktp[3], kp[3], Nt = 0, N = 0;
                for(size_t i = 0; i < 3; ++i) {
                    ktp[i] = std::max(0., kTilde[i]); // k tilde plus
                    Nt += ktp[i];
                    kp[i] = std::max(0., k[i]);
                    N += kp[i];
                }*/

                double fi = 0.;
                for (size_t j = 0; j < 3; ++j) {
                    fi += kTilde[j] * localValues[j] + (kTilde[j] - 2 * tr.getArea() / 3) * localPrevValues[j];
                }

                //std::cout << "Trojkat " << i + 1 << "\n\tfi = " << fi << "\n\n";
                for (size_t j = 0; j < 3; ++j) {
#pragma omp atomic
                    nu[mMesh.getTriangles()[i].vertices[j]] += beta[j] * fi;
#pragma omp atomic
                    si[mMesh.getTriangles()[i].vertices[j]] += tr.getArea() / 3;
#pragma omp atomic
                    ti[mMesh.getTriangles()[i].vertices[j]] += dt * std::max(0., k[j]) / 2 + tr.getArea() / 3;
#pragma omp atomic
                    tem[mMesh.getTriangles()[i].vertices[j]] += kTilde[j];

                    //if(mMesh.getTriangles()[i].vertices[j] == 208)
                    //if(i == 2)
                    //std::cout << "\tPunkt lokalny " << j << ", globalny " << mMesh.getTriangles()[i].vertices[j] + 1 << "\n\t\t" << "k = " << std::setw(5) << k[j] << '\t' << "k~ = " << kTilde[j] << '\t' << "k^ = " << kTilde[j] - 2 * tr.getArea() / 3 << "  \t" << "beta = " << beta[j] << "  \tdistributed fi = " << beta[j] * fi << std::endl;
                }
                //std::cout << std::endl;

                for (size_t j = 0; j < 3; ++j) {
                    if (k[j] > 0.) {
                        double a = 2 * tr.getArea() / (3 * k[j]);
#pragma omp critical
                        ddt = std::min(a, ddt);
                    }
                }
            }

            change = 0;
            for (size_t i = 0; i < mMesh.getPoints().size(); ++i) {
                //prevValues[i] = values[i];
                change += std::abs(nu[i] / tem[i]);
                values[i] -= nu[i] / tem[i];
                //values[i] -= nu[i] * ti[i] / si[i];
                //if(i == 209)
                //std::cout << values[i] << ' ' << nu[i] << ' ' << 1. / tem[i] << ' ' << ti[i] / si[i] << std::endl;
            }
            change /= mMesh.getPoints().size();
            //std::cout << it << ' ' << change << std::endl;
            //std::cout << 0.9 * ddt << ' ' << change << std::endl;
            ++it;
        }

        for (size_t i = 0; i < mMesh.getPoints().size(); ++i) {
            prevValues[i] = values[i];
        }
        std::cout << p << ' ' << it << std::endl;
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
