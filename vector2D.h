#pragma once
class vector2D
{
public:
	double x;
	double y;
	vector2D() = default;
	vector2D(double x, double y);
	void normalize();
    vector2D & operator *=(double m);
};

double dot(const vector2D &a, const vector2D &b);
