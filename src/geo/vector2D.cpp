#include "vector2D.h"
#include <cmath>
#include <iostream>

vector2D::vector2D(const vector2D &v) : x(v.x), y(v.y) {}
vector2D::vector2D(double x, double y) : x(x), y(y) {}

void vector2D::normalize() {
    double len = length();
    x /= len;
    y /= len;
}

vector2D & vector2D::operator *=(double m) {
    x *= m;
    y *= m;
    return *this;
}

double vector2D::length() {
    return std::sqrt(x * x + y * y);
}

double dotProduct(const vector2D &a, const vector2D &b) {
    return a.x * b.x + a.y * b.y;
}
