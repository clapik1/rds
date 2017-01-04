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
    timeDepSolverExplicit(mesh &mMesh, vector2D &advection, timeDepMethods method);
    void solve(double t, double (*wallElemValue)(double, double), unsigned int iterTotal);
    void animate(double t, double (*wallElemValue)(double, double), unsigned int iterTotal, std::ostream &solutionStream);
    double calcError(double t, double (*initialValue)(double, double),
                         double (*wallElemValue)(double, double));
private:
    mesh& mMesh;
    vector2D advection;
    timeDepMethods method;

    std::array<double, 3> solveTriangle(std::array<point2D, 3> &coords, std::array<double, 3> &localValues,
                                            std::array<double, 3> &localPrevValues, double dt) const;

    void timeStep(std::vector<double> &prevValues, double (*wallElemValue)(double, double), double dt);
};

#endif //RDS_TIMEDEPSOLVEREXPLICIT_H
