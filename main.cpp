#include <iostream>
#include <fstream>
#include "mesh.h"
#include "solver.h"

int main(int argc, char *argv[]) {
    //std::ifstream ifs(argv[1], std::ifstream::in);
    std::ifstream ifs("/home/clapik/workspace/meshes/grid_524.msh2", std::ifstream::in);
    mesh mMesh;
    if (!mMesh.init(ifs)) {
        std::cout << "error initializing mesh" << std::endl;
        return 0;
    }
    ifs.close();
    vector2D advection(1,1);
    double dt = 0.001;

    solver mSolver(mMesh, advection);
    mSolver.values[159] = 1;
    mSolver.solve(1, dt, methodRDS::N);
    std::ofstream ofs("/home/clapik/workspace/temp/out2.dat", std::ofstream::out);
    mSolver.toTecplot(ofs);
    ofs.close();

    /*for(int i = 0; i < mMesh.trianglesCount; ++i) {
        for(int j = 0; j < 3; ++j) {
            std::cout << mMesh.triangles[i].norm[j].x << ' ' << mMesh.triangles[i].norm[j].y << std::endl;
        }
        std::cout << std::endl;
    }*/
    return 0;
}
