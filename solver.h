#ifndef RDS_SOLVER_H
#define RDS_SOLVER_H


#include <ostream>
#include <istream>
#include <array>
#include "mesh.h"
#include "geo/vector2D.h"

enum class methodRDS {
    N,
    LDA,
    Blended,
    LimitedN
};

class solver {
public:
    solver(std::istream &meshStream, vector2D &advection);
    void statSolve(double dt, double (*wallElemValue)(double, double), double ghostHeight, methodRDS method);
    void toTecplot(std::ostream &os) const;
    std::vector<double> values;
private:
    mesh mMesh;
    vector2D advection;
    std::array<double, 3> solveTriangle(std::array<point2D, 3>& coords, std::array<double, 3>& localValues, methodRDS method) const;
};


#endif //RDS_SOLVER_H
