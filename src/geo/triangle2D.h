#ifndef RDS_TRIANGLE2D_H
#define RDS_TRIANGLE2D_H


#include <cstddef>
#include <array>
#include "vector2D.h"
#include "point2D.h"

class triangle2D {
public:
    size_t vertices[3];
    void updateArea(std::array<point2D, 3>& coords);
    double getArea() const;
private:
    double area;
};


#endif //RDS_TRIANGLE2D_H
