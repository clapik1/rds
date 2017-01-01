#ifndef RDS_TIMEDEPSOLVEREXPLICIT_H
#define RDS_TIMEDEPSOLVEREXPLICIT_H


#include <array>
#include <vector>
#include "../../geo/point2D.h"
#include "../../geo/vector2D.h"
#include "../../mesh.h"
#include "functions/timeDepFunctions.h"

class timeDepSolverExplicit {
public:
    timeDepSolverExplicit(mesh &mMesh, vector2D &advection, timeDepMethods method, double dt);
    void solve(double t, double (*wallElemValue)(double, double));
    double calcError(double t, double (*initialValue)(double, double),
                         double (*wallElemValue)(double, double));
private:
    mesh& mMesh;
    vector2D advection;
    timeDepMethods method;
    double dt;

    std::array<double, 3> solveTriangle(std::array<point2D, 3> &coords, std::array<double, 3> &localValues,
                                        std::array<double, 3> &localPrevValues) const;
};

#endif //RDS_TIMEDEPSOLVEREXPLICIT_H
