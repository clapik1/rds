//
// Created by clapik on 11.04.16.
//

#ifndef RDS_TRIANGLE2D_H
#define RDS_TRIANGLE2D_H

#include "vector2D.h"

class triangle2D {
public:
    int vertices[3];
    vector2D norm[3];
    double lengths[3];
    double getArea();
};


#endif //RDS_TRIANGLE2D_H
