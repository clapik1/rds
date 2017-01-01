#include <array>
#include <vector>
#include <limits>
#include "mesh.h"
#include "timeDepFunctions.h"
#include "steadyFunctions.h"

double calcMaxDt(mesh& mMesh, vector2D advection) {
    double dt = std::numeric_limits<double>::max();
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
    return dt;
}

std::array<double, 3> calcKTilde(double dt, double area, std::array<double, 3>& k) {
    std::array<double, 3> kt;

    for(size_t i = 0; i < 3; ++i) {
        kt[i] = dt * k[i] / 2 + area / 3; // k tilde
    }

    return kt;
};

std::array<double, 3> calcKHat(double dt, double area, std::array<double, 3>& k) {
    std::array<double, 3> kt;

    for(size_t i = 0; i < 3; ++i) {
        kt[i] = dt * k[i] / 2 - area / 3; // k tilde
    }

    return kt;
};

double calcNST(std::array<double, 3> &kTilde, std::array<double, 3> &kHat) {
    double NSTInv = 0.;
    for(size_t i = 0; i < 3; ++i) {
        NSTInv += std::max(0., kTilde[i]) + std::max(0., kHat[i]);
    }
    return 1. / NSTInv;
}

double calcInflowST(double NST, std::array<double, 3>& kTilde, std::array<double, 3>& kHat, std::array<double, 3> &localValues, std::array<double, 3> &localPrevValues) {
    double inflowST = 0.;
    for(size_t i = 0; i < 3; ++i) {
        inflowST += std::min(0., kTilde[i]) * localValues[i] + std::min(0., kHat[i]) * localPrevValues[i];
    }
    inflowST *= -NST;
    return inflowST;
}

double calcOutflowST(double NST, std::array<double, 3>& kTilde, std::array<double, 3>& kHat, std::array<double, 3> &localValues, std::array<double, 3> &localPrevValues) {
    double outflowST = 0.;
    for(size_t i = 0; i < 3; ++i) {
        outflowST += std::max(0., kTilde[i]) * localValues[i] + std::max(0., kHat[i]) * localPrevValues[i];
    }
    outflowST *= NST;
    return outflowST;
}
