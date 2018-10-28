#include <iostream>
#include <cinttypes>
#include <cmath>

#define ACCURACY_ROOT 1e-6
#define ITERATIONS_ROOT 2000

#define ITERATIONS_DERIVATIVE 30
#define START_STEP 1e-3

#define INF DBL_MAX * DBL_MAX
//  any error returns inf
//  df/dx(x);
double Derivative(double (*)(double), double);
//  f(x) = x, [left; right]
double IterationRoot(double(*)(double), double, double);
//  f(x) = 0, [left; right]
double NewTonMethod(double(*)(double), double, double);

int main()
{
	double(*func)(double) = [](double x){return (1.0 / sqrt(x + 1.0)); };
	Derivative(func, 1.0);
	IterationRoot(func, 0.1, 1.0);
	NewTonMethod(func, 0.1, 1.0);
	return 0;
}


double Derivative(double(*f)(double), double x)
{
	double delta = START_STEP;
	double res = (f(x + delta) - f(x)) / delta;
	delta /= 2.0;
	double var = (f(x + delta) - f(x)) / delta;
	bool Greater = var > res;
	for (int i = 0; i < ITERATIONS_DERIVATIVE; ++i)
	{
		delta /= 2.0;
		res = (f(x + delta) - f(x)) / delta;
		delta /= 2.0;
		var = (f(x + delta) - f(x)) / delta;
		if ((var == res) || ((var > res) != Greater))
			break;
		var = res;
	}
	return res;
}

double IterationRoot(double(*f)(double), double left, double right)
{
	if (left < right)
	{
		double first = (right - left) / 2.0;
		double next = f(first);
		for (int i = 0; i < ITERATIONS_ROOT; ++i)
		{
			if (abs(next - f(next)) < ACCURACY_ROOT)
				break;
			first = next;
			next = f(first);
		}
		if ((next >= left) && (next <= right))
			return next;
		else
			return INF;
	}
	else if ((left == right) && (f(left) == left))
		return left;
	else
		return INF;
}

double NewTonMethod(double(*f)(double), double left, double right)
{
	if (left < right)
	{
		double first = (right - left) / 2.0;
		double next = first - (f(first) - first)/ (Derivative(f, first) -1.0);
		for (int i = 0; i < ITERATIONS_ROOT; ++i)
		{
			if (abs(next - f(next)) < ACCURACY_ROOT)
				break;
			first = next;
			next = first - (f(first) - first) / (Derivative(f, first) - 1.0);
		}
		if ((next >= left) && (next <= right))
			return next;
		else
			return IterationRoot(f, left, right);
	}
	else if ((left == right) && (f(left) == left))
		return left;
	else
		return INF;
}
