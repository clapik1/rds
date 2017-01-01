#ifndef RDS_STEADYSOLVER_H
#define RDS_STEADYSOLVER_H


#include <ostream>
#include <istream>
#include <array>
#include "mesh.h"
#include "geo/vector2D.h"
#include "functions/steadyFunctions.h"

class statSolver {
public:
    statSolver(mesh &mesh, vector2D &advection, steadyMethods method);
    void solve(double (*wallElemValue)(double, double));
    double calcError(double (*wallElemValue)(double, double));
private:
    mesh &mMesh;
    vector2D advection;
    steadyMethods method;
    std::array<double, 3> solveTriangle(std::array<point2D, 3> &coords, std::array<double, 3> &localValues) const;
};


#endif //RDS_STEADYSOLVER_H
