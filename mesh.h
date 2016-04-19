#ifndef RDS_MESH_H
#define RDS_MESH_H


#include <string>
#include <vector>
#include "geo/point2D.h"
#include "geo/vector2D.h"
#include "geo/triangle2D.h"

struct wall2D {
    int nr;
    int vertices[2];
};

class mesh {
public:
    bool init(std::istream &ifs);
    size_t getPointsCount() const;
    size_t getTrianglesCount() const;
    size_t getWallsCount() const;
    const triangle2D& triangle(size_t nr) const;
    const point2D& point(size_t nr) const;
private:
    size_t pointsCount;
    std::vector<point2D> points;
    size_t trianglesCount;
    std::vector<triangle2D> triangles;
    size_t wallsCount;
    std::vector<wall2D> walls;
};


#endif //RDS_MESH_H
