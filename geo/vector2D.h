#ifndef RDS_VECTOR2D_H
#define RDS_VECTOR2D_H


class vector2D
{
public:
    double x;
    double y;
    vector2D() = default;
    vector2D(vector2D &v);
    vector2D(double x, double y);
    void normalize();
    vector2D & operator *=(double m);
};

double dot(const vector2D &a, const vector2D &b);


#endif //RDS_VECTOR2D_H
