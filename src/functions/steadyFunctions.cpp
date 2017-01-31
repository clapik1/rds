#include <vector>
#include "steadyFunctions.h"

std::array<double, 3> calcK(std::array<point2D, 3>& coords, vector2D advection) {
    std::array<double, 3> k;
    std::vector<vector2D> norm(3);

    for(size_t i = 0; i < 3; ++i) {
        double ax = coords[(i + 2) % 3].x;
        double ay = coords[(i + 2) % 3].y;
        double bx = coords[(i + 1) % 3].x;
        double by = coords[(i + 1) % 3].y;
        norm[i].x = by - ay;
        norm[i].y = ax - bx;
    }

    for(size_t i = 0; i < 3; ++i) {
        k[i] = dotProduct(advection, norm[i]) / 2;
    }

    return k;
}

double calcN(std::array<double, 3> &k) {
    double NInv = 0.;
    for(size_t i = 0; i < 3; ++i) {
        NInv += std::max(0., k[i]);
    }
    return 1. / NInv;
}

double calcInflow(double N, std::array<double, 3>& k, std::array<double, 3> &localValues) {
    double inflow = 0.;
    for(size_t i = 0; i < 3; ++i) {
        inflow += std::min(0., k[i]) * localValues[i];
    }
    inflow *= -N;
    return inflow;
}

double calcOutflow(double N, std::array<double, 3>& k, std::array<double, 3> &localValues) {
    double outflow = 0.;
    for(size_t i = 0; i < 3; ++i) {
        outflow += std::max(0., k[i]) * localValues[i];
    }
    outflow *= N;
    return outflow;
}

std::array<double, 3> distributeN(std::array<double, 3> &k, std::array<double, 3> &localValues, double inflow) {
    std::array<double, 3> delta;
    for (size_t i = 0; i < 3; ++i) {
        delta[i] = std::max(0., k[i]) * (localValues[i] - inflow);
    }
    return delta;
}

std::array<double, 3> distributeLDA(std::array<double, 3> &k, double inflow, double outflow) {
    std::array<double, 3> delta;
    for (size_t i = 0; i < 3; ++i) {
        delta[i] = std::max(0., k[i]) * (outflow - inflow);
    }
    return delta;
}

std::array<double, 3> applyMapping(double fi, std::array<double, 3> &old_dist) {
    std::array<double, 3> dist;
    if(fi == 0.) {
        for (size_t i = 0; i < 3; ++i) {
            dist[i] = 0.;
        }
    }
    else {
        double beta_plus[3];
        double beta_sum = 0.;
        int sign = (fi > 0.)? 1 : -1;
        for (size_t i = 0; i < 3; ++i) {
            beta_plus[i] = std::max(0., sign * old_dist[i]) + sign * 1e-10 * fi;
            beta_sum += beta_plus[i];
        }
        for (size_t i = 0; i < 3; ++i) {
            dist[i] = fi * beta_plus[i] / beta_sum;
        }
    }
    return dist;
}
