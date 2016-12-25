#ifndef RDS_UNSTATSOLVER_H
#define RDS_UNSTATSOLVER_H


#include <array>
#include <vector>
#include "../src/geo/point2D.h"
#include "../src/geo/vector2D.h"
#include "../src/mesh.h"
#include "../src/functions.h"

class unstatSolver {
public:
    unstatSolver(std::istream &meshStream, vector2D &advection, methodUnstat method);
    void unstatSolve(double t, double (*wallElemValue)(double, double));
    void toTecplot(std::ostream &os) const;
    std::vector<double> values;
private:
    methodUnstat method;
    mesh mMesh;
    vector2D advection;
    std::vector<std::vector<double>> matrix;
};


#endif //RDS_UNSTATSOLVER_H
