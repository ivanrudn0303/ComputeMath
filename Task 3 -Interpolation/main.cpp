#include "Interpolator.h"
#include <iostream>


int main()
{
	Interpolator f;
	std::vector<Point> func = { {-0.4, 1.9823}, {-0.1, 1.6710}, {0.2, 1.3694}, {0.5, 1.0472}, {0.8, 0.64350} };
	std::vector<Point> func2 = { {-0.1, 1.6710}, {0.2, 1.3694} };
	f.SetUp(func);
	std::cout << "Interpolated on all span f(0.1) = " << f(0.1) << std::endl;
	f.SetUp(func2);
	std::cout << "Interpolated on [-0.1; 0.2] span f(0.1) = " << f(0.1) << std::endl;
	system("PAUSE");
	return 0;
}
