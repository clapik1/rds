#ifndef RDS_STEADYFUNCTIONS_H
#define RDS_STEADYFUNCTIONS_H


#include <array>
#include "geo/point2D.h"
#include "geo/vector2D.h"

enum class steadyMethods {
    N,
    LDA,
    Blended,
    LimitedN
};

std::array<double, 3> calcK(std::array<point2D, 3>& coords, vector2D advection);
double calcN(std::array<double, 3> &k);
double calcInflow(double N, std::array<double, 3>& k, std::array<double, 3> &localValues);
double calcOutflow(double N, std::array<double, 3>& k, std::array<double, 3> &localValues);
std::array<double, 3> distributeN(std::array<double, 3> &k, std::array<double, 3> &localValues, double inflow);
std::array<double, 3> distributeLDA(std::array<double, 3> &k, double inflow, double outflow);
std::array<double, 3> applyMapping(double fi, std::array<double, 3> &old_dist);


#endif //RDS_STEADYFUNCTIONS_H
