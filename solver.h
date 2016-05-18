#ifndef RDS_SOLVER_H
#define RDS_SOLVER_H


#include <ostream>
#include "mesh.h"
#include "geo/vector2D.h"

enum class methodRDS {
    N,
    LDA
};

class solver {
public:
    solver(mesh &mMesh, vector2D &advection);
    void solveTriangle(const triangle2D &tr, methodRDS method) const;
    void solve(double t, double dt, methodRDS method);
    void toTecplot(std::ostream &os) const;
    std::vector<double> values;
private:
    mesh mMesh;
    vector2D advection;
};


#endif //RDS_SOLVER_H
