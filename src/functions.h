#ifndef RDS_FUNCTIONS_H
#define RDS_FUNCTIONS_H


#include <array>
#include "geo/point2D.h"
#include "geo/vector2D.h"

enum class methodStat {
    N,
    LDA,
    Blended,
    LimitedN
};

enum class methodUnstat {
    N,
    N_ST,
    LDA,
    LDA_ST
};

std::array<double, 3> calcK(std::array<point2D, 3>& coords, vector2D advection);
std::array<double, 3> calcKTilde(double area, std::array<double, 3>& k);
std::array<double, 3> calcUnstatBeta(double area, std::array<double, 3>& k, methodUnstat method);


#endif //RDS_FUNCTIONS_H
