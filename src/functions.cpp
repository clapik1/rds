#include <array>
#include <vector>
#include "functions.h"
#include "unstatSolver.h"
#include "constants.h"

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
};

std::array<double, 3> calcKTilde(double area, std::array<double, 3>& k) {
    std::array<double, 3> kt;

    for(size_t i = 0; i < 3; ++i) {
        kt[i] = dt * k[i] / 2 + area / 3; // k tilde
    }

    return kt;
};

std::array<double, 3> calcUnstatBeta(double area, std::array<double, 3>& k, methodUnstat method) {
    std::array<double, 3> beta, kt;
    double ktp[3], kp[3], Nt = 0, N = 0;

    kt = calcKTilde(area, k);
    for(size_t i = 0; i < 3; ++i) {
        ktp[i] = std::max(0., kt[i]); // k tilde plus
        Nt += ktp[i];
        kp[i] = std::max(0., k[i]);
        N += kp[i];
    }
    Nt = 1 / Nt; // N tilde
    N = 1 / N;

    switch(method) {
        case methodUnstat::N:
            break;
        case methodUnstat::N_ST:
            break;
        case methodUnstat::LDA:
            for(size_t i = 0; i < 3; ++i) {
                beta[i] = kp[i] * N;
            }
            break;
        case methodUnstat::LDA_ST:
            for(size_t i = 0; i < 3; ++i) {
                beta[i] = ktp[i] * Nt;
            }
            break;
    }
    return beta;
};
