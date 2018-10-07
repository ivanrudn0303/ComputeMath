#include <iostream>
#include <cinttypes>
#define ACCURACY 1e-3
#define ITERATIONS 20

double* GaussMatr(double *, uint32_t);
double* ZeidelMatr(const double *, uint32_t);

int main()
{
	uint16_t N = 0;
	std::cin >> N;
	double *data = new double[N * (N + 1)];
	for (uint16_t i = 0; i < N + 1; i++)
	{
		std::cout <<"Vvod vectora\n";
		for (uint16_t j = 0; j < N; j++)
			std::cin >> data[i * N + j];
	}
//	float da[12] = {1, 2, 3, 1, 3, 7, 8, 2, 1, 0, 0, 0};
	double* zeidresult = ZeidelMatr(data, N);
	if (zeidresult)
	{
		for (uint16_t j = 0; j < N; j++)
			std::cout << zeidresult[j] << '\n';
		delete[] zeidresult;
	}
	else
		std::cout << "Error Zeidel\n";
	double* result = GaussMatr(data, N);
	if (result)
	{
		for (uint16_t j = 0; j < N; j++)
			std::cout << result[j]<<'\n';
		delete[] result;
	}
	else
		std::cout << "Error Gauss\n";
	system("PAUSE");
	delete[] data;
	return 0;
}

double * GaussMatr(double * Data, uint32_t Dimension)
{
	if ((Data != nullptr) && Dimension)
	{
		for (uint32_t i = 0; i < Dimension; i++)
		{
			for (uint32_t j = Dimension - 1; j > i; j--)
			{
				double coeff = -Data[Dimension * i + j] / Data[Dimension * i + i];
				for (uint32_t k = i; k < Dimension + 1; k++)
					Data[k * Dimension + j] += Data[k * Dimension + i] * coeff;
			}
		}
		double *res = new double[Dimension], sum = 0.0;
		for (int64_t i = Dimension - 1; i >= 0; i--)
		{
			sum = Data[Dimension * Dimension + i];
			for (uint32_t j = 1; j < Dimension - i; j++)
			{
				sum -= Data[(Dimension - j) * Dimension + i] * res[Dimension - j];
			}
			if (Data[i * Dimension + i] == 0)
			{
				delete[] res;
				return nullptr;
			}
			res[i] = sum / Data[i * Dimension + i];
		}
		return res;
	}
	else
		return nullptr;
}

double * ZeidelMatr(const double *Data, uint32_t Dimension)
{
	if ((Data != nullptr) && Dimension)
	{
		double *res = new double[Dimension], *p = new double[Dimension];
		for (uint32_t i = 0; i < Dimension; i++)
			res[i] = 0.0;
		double SumSquare = 0.0;
		uint32_t counter = 0;
		do
		{
			SumSquare = 0.0;
			for (int i = 0; i < Dimension; i++)
				p[i] = res[i];
			for (uint32_t i = 0; i < Dimension; i++)
			{
				double var = 0;
				for (int j = 0; j < i; j++)
					var += (Data[j * Dimension + i] * res[j]);
				for (int j = i + 1; j < Dimension; j++)
					var += (Data[j * Dimension + i] * p[j]);
				res[i] = (Data[Dimension * Dimension + i] - var) / Data[i * Dimension + i];
				SumSquare += (res[i] - p[i]) * (res[i] - p[i]);
				counter++;
			}

		} while ((SumSquare > ACCURACY) && (counter < ITERATIONS));
		delete[] p;
		return res;
	}
	else
	{
		return nullptr;
	}
}
