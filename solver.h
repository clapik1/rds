#ifndef RDS_SOLVER_H
#define RDS_SOLVER_H


#include <ostream>
#include <istream>
#include "mesh.h"
#include "geo/vector2D.h"

enum class methodRDS {
    N,
    LDA
};

class solver {
public:
    solver(std::istream &meshStream, vector2D &advection);
    void solve(double t, double dt, double ghostHeight, methodRDS method);
    void toTecplot(std::ostream &os) const;
    std::vector<double> values;
    std::vector<double> borderValues;
private:
    mesh mMesh;
    vector2D advection;
    std::vector<double> solveTriangle(const std::vector<vector2D> &norm, const std::vector<double> &localValues, methodRDS method) const;
};


#endif //RDS_SOLVER_H
