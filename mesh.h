#ifndef RDS_MESH_H
#define RDS_MESH_H


#include <string>
#include <vector>
#include "geo/point2D.h"
#include "geo/vector2D.h"
#include "geo/triangle2D.h"
#include "geo/ghost2D.h"

struct wall2D {
    int wallNr;
    int vertices[2];
};

class mesh {
public:
    mesh(std::istream &ifs);
    const std::vector<triangle2D>& getTriangles() const;
    const std::vector<point2D>& getPoints() const;
    const std::vector<wall2D>& getWallElems() const;
private:
    std::vector<point2D> points;
    std::vector<triangle2D> triangles;
    std::vector<wall2D> wallElems;
};


#endif //RDS_MESH_H
