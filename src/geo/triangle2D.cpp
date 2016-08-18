#include "triangle2D.h"
#include <cmath>

double triangle2D::getArea() const {
    double p = (lengths[0] + lengths[1] + lengths[2]) / 2;
    return sqrt(p * (p - lengths[0]) * (p - lengths[1]) * (p - lengths[2]));
}
