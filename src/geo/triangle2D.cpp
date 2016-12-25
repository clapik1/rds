#include "triangle2D.h"
#include <cmath>

void triangle2D::updateArea(std::array<point2D, 3>& coords) {
    double lengths[3];
    for (size_t j = 0; j < 3; ++j) {
        lengths[j] = std::sqrt(std::pow(coords[(j + 2) % 3].x - coords[(j + 1) % 3].x, 2.) + std::pow(coords[(j + 2) % 3].y - coords[(j + 1) % 3].y, 2.));
    }
    double p = (lengths[0] + lengths[1] + lengths[2]) / 2;
    double temp = p * (p - lengths[0]) * (p - lengths[1]) * (p - lengths[2]);
    if(temp < 0.) area = 0.;
    else area = std::sqrt(temp);
}

double triangle2D::getArea() const {
    return area;
}
