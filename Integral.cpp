#include <iostream>
#include <cmath>
#include <cstdint>
#define F4MAX 249.0/32.0
#define EPS 1e-3

int main()
{
	double a = 1.0, b = 16.0 / EPS / EPS;
	double step = 1.0;
	double Integral = 0.0;
	step = pow(EPS / 2.0 * 2880.0 / F4MAX / (b - a), 0.25);
	std::cout << "Integrational Step = " << step << "\n";
	for (uint64_t i = 0; i < uint64_t((b - a) / step); ++i)
	{	
		double x = a + i * step;
		Integral += step / 6.0 * (1.0 / sqrt(x) / (x + 1.0) + 4.0 / sqrt(x + step / 2.0) / (x + 1.0 + step / 2.0) + 1.0 / sqrt(x + step) / (x + 1.0 + step));
		std::cout << (x - a) / (b - a) * 100.0 << "\n";
	}
	std::cout << "Integral = " << Integral << "\n";
	system("PAUSE");
	return 0;
}
