#ifndef RDS_GHOST2D_H
#define RDS_GHOST2D_H


#include <cstddef>

enum class ghostType {
    left,
    right
};

class ghost2D {
    size_t wallNr;
    double height;
    ghostType type;
};


#endif //RDS_GHOST2D_H
