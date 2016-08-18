#ifndef RDS_POINT2D_H
#define RDS_POINT2D_H


class point2D
{
public:
	point2D() = default;
	point2D(double x, double y);
	double x;
	double y;

	/*point2D &operator+=(point2D p) {
		this->x += p.x;
		this->y += p.y;
		return *this;
	}

	point2D &operator+(point2D p) {
		return point2D(this->x + p.x, this->y + p.y);
	}*/
};


#endif //RDS_POINT2D_H
