#include "Interpolator.h"
#include "Gauss.h"

uint32_t Interpolator::Find(double x) const
{
	if(x < Points[0].x)
		return 0;
	if (x > Points[Size - 1].x)
		return Size - 2;
	uint32_t min = 0;
	uint32_t max = Size - 1;
	while (max - min > 1)
	{
		uint32_t mid = (max + min) / 2;
		if (x > Points[mid].x)
		{
			min = mid;
		}
		else
		{
			max = mid;
		}
	}
	return min;

}

Interpolator::Interpolator(const std::vector<Point>& Array): Coefficients(nullptr), Points(Array), Size(Array.size()), Ready(false)
{}

Interpolator::~Interpolator()
{
	delete[] Coefficients;
}

Interpolator::Interpolator():Ready(false), Coefficients(nullptr), Size(0)
{}

int Interpolator::SetUp()
{
	if(Ready)
		return 0;
	if (Size < 2)
		return 1;
	uint32_t N = 4 * (Size - 1);
	std::vector<double> matrix(N * (N + 1), 0.0);
	double h = 0.0;
	delete[] Coefficients;
	Coefficients = nullptr;
	for (uint32_t i = 0; i < Size - 1; ++i)
	{
		//1
		matrix[4 * i * N + 4 * i] = 1;
		matrix[N * N + 4 * i] = Points[i].y;
		//2
		h = Points[i + 1].x - Points[i].x;
		matrix[4 * i * N + 4 * i + 1] = 1;
		matrix[(4 * i + 1) * N + 4 * i + 1] = h;
		matrix[(4 * i + 2) * N + 4 * i + 1] = h * h;
		matrix[(4 * i + 3) * N + 4 * i + 1] = h * h * h;
		matrix[N * N + 4 * i + 1] = Points[i + 1].y;
		//3
		if (i != Size - 2)
		{
			matrix[(4 * i + 1) * N + 4 * i + 2] = 1;
			matrix[(4 * i + 2) * N + 4 * i + 2] = 2.0 * h;
			matrix[(4 * i + 3) * N + 4 * i + 2] = 3.0 * h * h;
			matrix[(4 * i + 1 + 4) * N + 4 * i + 2] = -1;
			//			matrix[(4 * i + 2 + 4) * N + 4 * i + 2] = -2.0 * h;
			//			matrix[(4 * i + 3 + 4) * N + 4 * i + 2] = -3.0 * h * h;
			matrix[N * N + 4 * i + 2] = 0.0;
		}
		//4
		matrix[(4 * i + 2) * N + 4 * i + 3] = 2.0;
		matrix[(4 * i + 3) * N + 4 * i + 3] = 6.0 * h;
		if (i != Size - 2)
		{
			matrix[(4 * i + 2 + 4) * N + 4 * i + 3] = -2.0;
//			matrix[(4 * i + 3 + 4) * N + 4 * i + 3] = -6.0 * h;
		}
		matrix[N * N + 4 * i + 3] = 0.0;
	}
	//boundary condition
	matrix[2 * N + N - 2] = 2.0;
	matrix[N * N + N - 2] = 0.0;

	Coefficients = GaussMatr(matrix.data(), N);
//	Coefficients = ZeidelMatr(matrix.data(), N);
	if (Coefficients)
	{
		Ready = true;
		return 0;
	}
	else
	{
		Ready = false;
		return 1;
	}
}

int Interpolator::SetUp(const std::vector<Point>&array)
{
	Points = array;
	Size = array.size();
	if (Ready)
		Ready = false;
	delete[] Coefficients;
	Coefficients = nullptr;
	return SetUp();
}

double Interpolator::operator()(double x) const
{
	if(!Ready)
		return NAN;
	uint32_t i = Find(x);
	double h = x - Points[i].x;
	return Coefficients[i * 4] + Coefficients[i * 4 + 1] * h + Coefficients[i * 4 + 2] * h * h + Coefficients[i * 4 + 3] * h * h * h;
}
