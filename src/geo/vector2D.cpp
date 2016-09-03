#include "vector2D.h"
#include <cmath>

vector2D::vector2D(const vector2D &v) : x(v.x), y(v.y) {}
vector2D::vector2D(double x, double y) : x(x), y(y) {}

void vector2D::normalize() {
    double len = std::sqrt(x * x + y * y);
    x /= len;
    y /= len;
}

vector2D & vector2D::operator *=(double m) {
    x *= m;
    y *= m;
}

double dotProduct(const vector2D &a, const vector2D &b) {
    return a.x * b.x + a.y * b.y;
}
