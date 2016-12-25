#ifndef RDS_FUNCTIONS_H
#define RDS_FUNCTIONS_H


#include <array>
#include "geo/point2D.h"
#include "geo/vector2D.h"

enum class methodUnstat {
    N,
    N_ST,
    LDA,
    LDA_ST
};

double calcMaxDt(mesh mMesh, vector2D advection);
std::array<double, 3> calcKTilde(double dt, double area, std::array<double, 3>& k);
std::array<double, 3> calcKHat(double dt, double area, std::array<double, 3>& k);
double calcNST(std::array<double, 3> &kTilde, std::array<double, 3> &kHat);
double calcInflowST(double NST, std::array<double, 3>& kTilde, std::array<double, 3>& kHat, std::array<double, 3> &localValues, std::array<double, 3> &localPrevValues);
double calcOutflowST(double NST, std::array<double, 3>& kTilde, std::array<double, 3>& kHat, std::array<double, 3> &localValues, std::array<double, 3> &localPrevValues);

#endif //RDS_FUNCTIONS_H
