//
// Created by clapik on 11.04.16.
//

#include "triangle2D.h"
#include <cmath>

double triangle2D::getArea() {
    double p = (lengths[0] + lengths[1] + lengths[2]) / 2;
    return sqrt(p * (p - lengths[0]) * (p - lengths[1]) * (p - lengths[2]));
}
