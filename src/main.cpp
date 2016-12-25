#include <iostream>
#include <fstream>
#include <cmath>
#include "mesh.h"
#include "unstatSolver/implicit/unstatSolverImplicit.h"
#include "statSolver/statSolver.h"

void setValues(mesh& mMesh) {
    double r, ma = -1;
    for(size_t i = 0; i < mMesh.getValues().size(); ++i) {
        r = std::sqrt(std::pow(mMesh.getPoints()[i].x - (-0.25), 2.) + std::pow(mMesh.getPoints()[i].y - (-0.25), 2.));
        if(r < 0.1) {
            mMesh.setValue(i, std::pow(std::cos(2 * std::acos(-1.) * r * 2.5), 2.) * 0.4);
            ma = std::max(ma, mMesh.getValues()[i]);
        }
    }

    std::cout << "max value: " << ma << '\n';
}

int main(int argc, char *argv[]) {
    if(argc != 3) {
        std::cout << "Bad usage" << std::endl;
        return 1;
    }

    std::ifstream ifs(argv[1], std::ifstream::in);
    mesh mMesh(ifs);
    ifs.close();

    vector2D advection(0.5, 0.5);
    setValues(mMesh);
    double dt = 0.9 * calcMaxDt(mMesh, advection);
    std::cout << "dt = " << dt << '\n';
    //statSolver mSolver(mMesh, advection, methodStat::LimitedN);
    unstatSolverImplicit mSolver(mMesh, advection, methodUnstat::LDA, dt);

    //auto sinLambda = [](double x, double y) { return (std::sin(std::acos(-1.) * (x - y)) + 1) / 2; };
    //auto stepLambda = [](double x, double y) { return x > -0.5 ? 1. : 0.; };

    //mSolver.statSolve(sinLambda);
    //mSolver.statSolve(stepLambda);
    mSolver.unstatSolve(1., [](double x, double y) { return 0.; });


    std::ofstream ofs(argv[2], std::ofstream::out);
    mMesh.toTecplot(ofs);
    ofs.close();

    //std::cout << mSolver.statCheck(sinLambda) << std::endl;
    //std::cout << mSolver.statCheck(stepLambda) << std::endl;

    return 0;
}
