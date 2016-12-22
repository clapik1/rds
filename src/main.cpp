#include <iostream>
#include <fstream>
#include <cmath>
#include "mesh.h"
#include "constants.h"
#include "unstatSolverExplicit.h"

int main(int argc, char *argv[]) {
    if(argc != 3) {
        std::cout << "Bad usage" << std::endl;
        return 1;
    }

    std::ifstream ifs(argv[1], std::ifstream::in);
    vector2D advection(1, 1);
    unstatSolverExplicit mSolver(ifs, advection, methodUnstat::LDA_ST);
    ifs.close();

    //mSolver.values[159] = 1;

    auto sinLambda = [](double x, double y) { return (std::sin(std::acos(-1.) * (x - y)) + 1) / 2; };
    //auto stepLambda = [](double x, double y) { return x > -0.5 ? 1. : 0.; };

    //mSolver.statSolve(sinLambda);
    //mSolver.statSolve(stepLambda);
    mSolver.unstatSolve(0.04, sinLambda);


    std::ofstream ofs(argv[2], std::ofstream::out);
    mSolver.toTecplot(ofs);
    ofs.close();

    //std::cout << mSolver.statCheck(sinLambda) << std::endl;
    //std::cout << mSolver.statCheck(stepLambda) << std::endl;

    return 0;
}
