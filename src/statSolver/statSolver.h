#ifndef RDS_SOLVER_H
#define RDS_SOLVER_H


#include <ostream>
#include <istream>
#include <array>
#include "mesh.h"
#include "geo/vector2D.h"
#include "functions/statFunctions.h"

class statSolver {
public:
    statSolver(mesh &mesh, vector2D &advection, methodStat method);
    void statSolve(double (*wallElemValue)(double, double));
    double statCheck(double (*wallElemValue)(double, double));
private:
    mesh &mMesh;
    vector2D advection;
    methodStat method;
    std::array<double, 3> statSolveTriangle(std::array<point2D, 3> &coords, std::array<double, 3> &localValues) const;
};


#endif //RDS_SOLVER_H
