#include <set>
#include <iostream>
#include <cmath>
#include <iomanip>
#include "constants.h"
#include "unstatSolverImplicit.h"

unstatSolverImplicit::unstatSolverImplicit(std::istream &meshStream, vector2D &advection, methodUnstat method) : mMesh(meshStream), advection(advection), method(method) {
    values.resize(mMesh.getPoints().size());
    double r, ma = -1;
    for(size_t i = 0; i < values.size(); ++i) {
        r = std::sqrt(std::pow(mMesh.getPoints()[i].x - (-0.25), 2.) + std::pow(mMesh.getPoints()[i].y - (-0.25), 2.));
        if(r < 0.1) {
            values[i] = std::pow(std::cos(2 * std::acos(-1.) * r * 2.5), 2.) * 0.4;
            ma = std::max(ma, values[i]);
        }
    } //fixme delete this shit

    std::cout << "max value: " << ma << '\n';

    //values[8] = 1.;
}

std::array<double, 3> unstatSolverImplicit::unstatSolveTriangle(std::array<point2D, 3> &coords, std::array<double, 3> &localValues, std::array<double, 3> &localPrevValues) const {
    std::array<double, 3> delta;

    triangle2D tr;
    tr.updateArea(coords);

    std::array<double, 3> k = calcK(coords, advection);
    std::array<double, 3> kTilde = calcKTilde(dt, tr.getArea(), k);

    double fi = 0.;
    for (size_t i = 0; i < 3; ++i) {
        fi += kTilde[i] * localValues[i] + (kTilde[i] - 2 * tr.getArea() / 3) * localPrevValues[i];
    }
    if(fi == 0.) {
        for (size_t i = 0; i < 3; ++i) {
            delta[i] = 0.;
        }
    }
    else {
        std::array<double, 3> beta = calcUnstatBeta(dt, tr.getArea(), k, method);

        for (size_t i = 0; i < 3; ++i) {
            delta[i] = beta[i] * fi;
        }
    }

    return delta;
}

void unstatSolverImplicit::unstatSolve(double t, double (*wallElemValue)(double, double)) {
    std::vector<double> prevValues(values.size());
    dt = t;
    for(size_t i = 0; i < values.size(); ++i) {
        prevValues[i] = values[i];
    }

    for (size_t i = 0; i < mMesh.getTriangles().size(); ++i) {
        std::array<point2D, 3> coords;
        for (size_t j = 0; j < 3; ++j) {
            coords[j] = mMesh.getPoints()[mMesh.getTriangles()[i].vertices[j]];
        }
        triangle2D tr;
        tr.updateArea(coords);
        std::array<double, 3> k = calcK(coords, advection);
        for (size_t j = 0; j < 3; ++j) {
            if (k[j] > 0.) {
                double a = 2 * tr.getArea() / (3 * k[j]);
                dt = std::min(a, dt);
            }
        }
    }
    dt *= 0.9;
    std::cout << "dt: " << dt << '\n';

    for(auto p = 0.; p < t; p += dt) {
        double change = 1.;
        int it = 0;

        while(std::abs(change) > 1e-6) {
            std::vector<double> nu(mMesh.getPoints().size()), si(mMesh.getPoints().size()), tem(mMesh.getPoints().size());

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
                std::array<double, 3> kTilde = calcKTilde(dt, tr.getArea(), k);

                /*double ktp[3], kp[3], Nt = 0, N = 0;
                for(size_t i = 0; i < 3; ++i) {
                    ktp[i] = std::max(0., kTilde[i]); // k tilde plus
                    Nt += ktp[i];
                    kp[i] = std::max(0., k[i]);
                    N += kp[i];
                }*/

                delta = unstatSolveTriangle(coords, localValues, localPrevValues);

                //std::cout << "Trojkat " << i + 1 << "\n\tfi = " << fi << "\n\n";
                for (size_t j = 0; j < 3; ++j) {
                    #pragma omp atomic
                    nu[mMesh.getTriangles()[i].vertices[j]] += delta[j];
                    #pragma omp atomic
                    si[mMesh.getTriangles()[i].vertices[j]] += tr.getArea() / 3;
                    //#pragma omp atomic
                    //ti[mMesh.getTriangles()[i].vertices[j]] += dt * std::max(0., k[j]) / 2 + tr.getArea() / 3;
                    #pragma omp atomic
                    tem[mMesh.getTriangles()[i].vertices[j]] += kTilde[j];

                    //std::cout << "\tPunkt lokalny " << j << ", globalny " << mMesh.getTriangles()[i].vertices[j] + 1 << "\n\t\t" << "k = " << std::setw(5) << k[j] << '\t' << "k~ = " << kTilde[j] << '\t' << "k^ = " << kTilde[j] - 2 * tr.getArea() / 3 << "  \t" << "beta = " << beta[j] << "  \tdistributed fi = " << beta[j] * fi << std::endl;
                }
                //std::cout << std::endl;
            }

            #pragma omp parallel for
            for(size_t i = 0; i < mMesh.getWallElems().size(); ++i) {
                std::array<double, 3> localPrevValues, localValues, delta;
                std::array<point2D, 3> coords;

                for(size_t j = 0; j < 2; ++j) {
                    coords[j + 1] = mMesh.getPoints()[mMesh.getWallElems()[i].vertices[j]];
                    localValues[j + 1] = values[mMesh.getWallElems()[i].vertices[j]];
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

            /*for(size_t i = 0; i < mMesh.getWallElems().size(); ++i) {
                for(size_t j = 0; j < 2; ++j) {
                    nu[mMesh.getWallElems()[i].vertices[j]] = 0;
                }
            }*/

            change = 0;
            for (size_t i = 0; i < mMesh.getPoints().size(); ++i) {
                change += std::abs(nu[i] / tem[i]);
            }
            change /= mMesh.getPoints().size();

            for (size_t i = 0; i < mMesh.getPoints().size(); ++i) {
                values[i] -= nu[i] / tem[i];
                //values[i] -= nu[i] * ti[i] / si[i];
            }

            //std::cout << it << ' ' << change << std::endl;
            ++it;
        }

        /*for(size_t i = 0; i < mMesh.getWallElems().size(); ++i) {
            for(size_t j = 0; j < 2; ++j) {
                values[mMesh.getWallElems()[i].vertices[j]] = prevValues[mMesh.getWallElems()[i].vertices[j]];
            }
        }*/

        for (size_t i = 0; i < mMesh.getPoints().size(); ++i) {
            prevValues[i] = values[i];
        }
        std::cout << p << ' ' << it << std::endl;
    }
}

void unstatSolverImplicit::toTecplot(std::ostream &os) const {
    os << "TITLE = \"\"\nVARIABLES = \"X\", \"Y\", \"VALUE\"\nZONE T=\"\", N=" << mMesh.getPoints().size() << ", E=" << mMesh.getTriangles().size() << ", F=FEPOINT, ET=QUADRILATERAL\n";
    for(size_t i = 0; i < mMesh.getPoints().size(); ++i) {
        os << mMesh.getPoints()[i].x << ' ' << mMesh.getPoints()[i].y << ' ' << values[i] << '\n';
    }
    for(size_t i = 0; i < mMesh.getTriangles().size(); ++i) {
        os << mMesh.getTriangles()[i].vertices[0] + 1 << ' ' << mMesh.getTriangles()[i].vertices[1] + 1 << ' ' << mMesh.getTriangles()[i].vertices[2] + 1 << ' ' << mMesh.getTriangles()[i].vertices[2] + 1 << '\n';
    }
}
