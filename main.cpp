#include <iostream>
#include <fstream>
#include "mesh.h"

int main(int argc, char *argv[]) {
    //std::ifstream ifs(argv[1], std::ifstream::in);
    std::ifstream ifs("/home/clapik/workspace/meshes/grid_8.msh2", std::ifstream::in);
    mesh mMesh;
    vector2D a(1,1);
    double dt = 0.001;
    if (!mMesh.init(ifs)) {
        std::cout << "error initializing mesh" << std::endl;
        return 0;
    }
    double area = 0;
    for (int i = 0; i < mMesh.getTrianglesCount(); ++i) {
        area += mMesh.triangles[i].getArea();
    }
    std::cout << area << std::endl;
    mMesh.points[2].value = 1;
    for (int it = 0; it < 10; ++it) {
        std::vector<double> nu(mMesh.getPointsCount()), si(mMesh.getPointsCount());
        for(int i = 0; i < mMesh.trianglesCount; ++i) {
            double k[3], kp[3], beta[3], kpsum = 0, fi = 0;
            for(int j = 0; j < 3; ++j) {
                k[j] = dot(a, mMesh.triangles[i].norm[j]) / 2;
                kp[j] = std::max(0., k[j]);
                kpsum += kp[j];
                fi += k[j] * mMesh.points[mMesh.triangles[i].vertices[j]].value;
                si[mMesh.triangles[i].vertices[j]] += mMesh.triangles[i].getArea() / 3;
            }
            for(int j = 0; j < 3; ++j) {
                beta[j] = kp[j] / kpsum;
                nu[mMesh.triangles[i].vertices[j]] += beta[j] * fi;
            }
        }
        for(int i = 0; i < mMesh.pointsCount; ++i) {
            mMesh.points[i].value += dt * nu[i] / si[i];
        }
        for (int j = 0; j < mMesh.pointsCount; ++j) {
            std::cout << mMesh.points[j].value << ' ';
        }
        std::cout << std::endl;
    }
    /*for(int i = 0; i < mMesh.trianglesCount; ++i) {
        for(int j = 0; j < 3; ++j) {
            std::cout << mMesh.triangles[i].norm[j].x << ' ' << mMesh.triangles[i].norm[j].y << std::endl;
        }
        std::cout << std::endl;
    }*/
    return 0;
}
