#include <iostream>
#include <fstream>
#include "mesh.h"
#include "solver.h"

int main(int argc, char *argv[]) {
    //std::ifstream ifs(argv[1], std::ifstream::in);
    std::ifstream ifs("/home/clapik/workspace/meshes/grid_524.msh2", std::ifstream::in);
    vector2D advection(1,1);
    solver mSolver(ifs, advection);
    ifs.close();

    auto dt = 0.001;
    auto ghostHeight = 0.01;
    mSolver.values[159] = 1;
    mSolver.solve(1, dt, ghostHeight, methodRDS::N);

    std::ofstream ofs("/home/clapik/workspace/temp/out2.dat", std::ofstream::out);
    mSolver.toTecplot(ofs);
    ofs.close();

    return 0;
}
