#include <iostream>
#include <fstream>
#include <cmath>
#include "mesh.h"
#include "solver.h"

int main(int argc, char *argv[]) {
    if(argc != 3) {
        std::cout << "Bad usage" << std::endl;
        return 1;
    }
    std::ifstream ifs(argv[1], std::ifstream::in);
    //std::ifstream ifs("/home/clapik/workspace/meshes/grid_8334.msh2", std::ifstream::in);
    vector2D advection(1,1);
    solver mSolver(ifs, advection);
    ifs.close();

    auto dt = 0.001;
    auto ghostHeight = 0.01;
    mSolver.values[159] = 1;

    auto sinLambda = [](double x, double y) { return std::cos(2. * std::acos(-1.) * (x - y)); };
    //auto stepLambda = [](double x, double y) { return x > -0.5 ? 1. : 0.; };
    mSolver.statSolve(dt, sinLambda, ghostHeight, methodRDS::Blended);
    //mSolver.statSolve(dt, stepLambda, ghostHeight, methodRDS::LDA);

    std::ofstream ofs(argv[2], std::ofstream::out);
    //std::ofstream ofs("/home/clapik/workspace/temp/out_dev.dat", std::ofstream::out);
    mSolver.toTecplot(ofs);
    ofs.close();

    return 0;
}
