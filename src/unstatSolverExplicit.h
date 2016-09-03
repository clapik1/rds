#ifndef RDS_UNSTATSOLVER_H
#define RDS_UNSTATSOLVER_H


#include <array>
#include <vector>
#include "geo/point2D.h"
#include "geo/vector2D.h"
#include "mesh.h"
#include "functions.h"

class unstatSolverExplicit {
public:
    unstatSolverExplicit(std::istream &meshStream, vector2D &advection, methodUnstat method);
    void unstatSolve(double t, double (*wallElemValue)(double, double));
    void toTecplot(std::ostream &os) const;
    std::vector<double> values;
private:
    methodUnstat method;
    mesh mMesh;
    vector2D advection;
};


#endif //RDS_UNSTATSOLVER_H
