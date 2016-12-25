#ifndef RDS_UNSTATSOLVER_H
#define RDS_UNSTATSOLVER_H


#include <array>
#include <vector>
#include "../../geo/point2D.h"
#include "../../geo/vector2D.h"
#include "../../mesh.h"
#include "functions/unstatFunctions.h"

class unstatSolverImplicit {
public:
    unstatSolverImplicit(mesh &mesh, vector2D &advection, methodUnstat method, double dt);
    void unstatSolve(double t, double (*wallElemValue)(double, double));
private:
    methodUnstat method;
    mesh& mMesh;
    vector2D advection;
    double dt;

    std::array<double, 3> unstatSolveTriangle(std::array<point2D, 3> &coords, std::array<double, 3> &localValues,
                                              std::array<double, 3> &localPrevValues) const;
};

#endif //RDS_UNSTATSOLVER_H
