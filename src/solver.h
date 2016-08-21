#ifndef RDS_SOLVER_H
#define RDS_SOLVER_H


#include <ostream>
#include <istream>
#include <array>
#include "mesh.h"
#include "geo/vector2D.h"

enum class methodStat {
    N,
    LDA,
    Blended,
    LimitedN
};

enum class methodUnstat {
    N,
    N_ST,
    LDA,
    LDA_ST
};

class solver {
public:
    solver(std::istream &meshStream, vector2D &advection);
    void statSolve(double (*wallElemValue)(double, double), methodStat method);
    void unstatSolve(double t, double (*wallElemValue)(double, double), methodUnstat method);
    double statCheck(double (*wallElemValue)(double, double));
    void toTecplot(std::ostream &os) const;
    std::vector<double> values;
private:
    mesh mMesh;
    vector2D advection;
    std::array<double, 3> statSolveTriangle(std::array<point2D, 3> &coords, std::array<double, 3> &localValues, methodStat method) const;
    std::array<double, 3> unstatSolveTriangle(std::array<point2D, 3>& coords, std::array<double, 3>& prevValues, std::array<double, 3>& localValues, methodUnstat method) const;
};


#endif //RDS_SOLVER_H
