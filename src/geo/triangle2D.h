#ifndef RDS_TRIANGLE2D_H
#define RDS_TRIANGLE2D_H


#include <cstddef>
#include "vector2D.h"

class triangle2D {
public:
    size_t vertices[3];
    double lengths[3];
    double getArea() const;
};


#endif //RDS_TRIANGLE2D_H
