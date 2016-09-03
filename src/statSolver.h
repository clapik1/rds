#ifndef RDS_SOLVER_H
#define RDS_SOLVER_H


#include <ostream>
#include <istream>
#include <array>
#include "mesh.h"
#include "geo/vector2D.h"
#include "functions.h"

class statSolver {
public:
    statSolver(std::istream &meshStream, vector2D &advection, methodStat method);
    void statSolve(double (*wallElemValue)(double, double));
    double statCheck(double (*wallElemValue)(double, double));
    void toTecplot(std::ostream &os) const;
    std::vector<double> values;
private:
    mesh mMesh;
    vector2D advection;
    methodStat method;
    std::array<double, 3> statSolveTriangle(std::array<point2D, 3> &coords, std::array<double, 3> &localValues) const;
};


#endif //RDS_SOLVER_H
