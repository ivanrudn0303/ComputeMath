#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H
#include <vector>

struct Point
{
	double x, y;
};

class Interpolator
{
private:
	std::vector<Point> Points;
	bool Ready;
	double* Coefficients;
	uint32_t Size;
	uint32_t Find(double) const;
public:
	Interpolator(const std::vector<Point>&);
	~Interpolator();
	Interpolator();
	int SetUp();
	int SetUp(const std::vector<Point>&);
	double operator()(double) const;
};
#endif // ! INTERPOLATOR_H
