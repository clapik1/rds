#ifndef RDS_UNSTATSOLVER_H
#define RDS_UNSTATSOLVER_H


#include <array>
#include <vector>
#include "geo/point2D.h"
#include "geo/vector2D.h"
#include "mesh.h"
#include "functions.h"

class unstatSolverImplicit {
public:
    unstatSolverImplicit(std::istream &meshStream, vector2D &advection, methodUnstat method);

    void unstatSolve(double t, double (*wallElemValue)(double, double));

    void toTecplot(std::ostream &os) const;

    std::vector<double> values;
private:
    methodUnstat method;
    mesh mMesh;
    vector2D advection;
    double dt;

    std::array<double, 3> unstatSolveTriangle(std::array<point2D, 3> &coords, std::array<double, 3> &localValues,
                                              std::array<double, 3> &localPrevValues) const;
};

#endif //RDS_UNSTATSOLVER_H
