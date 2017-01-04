#include <iostream>
#include <fstream>
#include <cmath>
#include "helpers/argsParser.h"
#include "mesh.h"
#include "timeDepSolver/explicit/timeDepSolverExplicit.h"
#include "steadySolver/steadySolver.h"

double bump(double x, double y) {
    double r = std::sqrt(std::pow(x - (-0.25), 2.) + std::pow(y - (-0.25), 2.));
    if(r < 0.2)
        return std::pow(std::cos(2 * std::acos(-1.) * r * 1.25), 2.);
    return 0.;
}

void setValues(mesh& mMesh) {
    double mv = -1;
    for(size_t i = 0; i < mMesh.getValues().size(); ++i) {
        mMesh.setValue(i, bump(mMesh.getPoints()[i].x, mMesh.getPoints()[i].y));
        mv = std::max(mv, mMesh.getValues()[i]);
    }

    //std::cout << "max value: " << mv << '\n';
}

int main(int argc, char *argv[]) {
    argsParser args(argc, argv);
    std::string inputPath = args.getOption("-i");
    std::string outputPath = args.getOption("-o");
    std::string methodArg = args.getOption("-m");

    if(inputPath == "" || methodArg == "" || outputPath == "") {
        std::cout << "Bad usage\n";
        return 1;
    }

    std::ifstream ifs(inputPath, std::ifstream::in);
    mesh mMesh(ifs);
    ifs.close();

    vector2D advection(0.5, 0.5);

    //steadySolver mSolver(mMesh, advection, steadyMethods::LimitedN);

    setValues(mMesh);
    auto maxDt = 0.9 * calcMaxDt(mMesh, advection);
    auto t = 1.;
    auto iterTotal = static_cast<uint32_t>(std::ceil(t / maxDt));
    std::cout << "maxDt = " << maxDt << "    dt = " << iterTotal << '\n';
    timeDepSolverExplicit mSolver(mMesh, advection, static_cast<timeDepMethods>(std::stoi(methodArg)));

    //auto sinLambda = [](double x, double y) { return (std::sin(std::acos(-1.) * (x - y)) + 1) / 2; };
    //auto stepLambda = [](double x, double y) { return x > -0.5 ? 1. : 0.; };
    auto zeroLambda = [](double x, double y) { return 0.; };

    //mSolver.solve(sinLambda);
    //mSolver.solve(stepLambda);

    std::ofstream ofs(outputPath, std::ofstream::out);
    //mMesh.toTecplot(ofs);
    mSolver.animate(t, zeroLambda, iterTotal, ofs);
    ofs.close();

    //std::cout << mSolver.calcError(sinLambda) << '\n';
    //std::cout << mSolver.calcError(stepLambda) << '\n';
    std::cout << mSolver.calcError(t, bump, zeroLambda) << '\n';

    return 0;
}
