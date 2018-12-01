#ifndef GAUSS_H
#define GAUSS_H
#include <cstdint>
#define ACCURACY 1e-6
#define ITERATIONS 2000
double* GaussMatr(double *, uint32_t);
double* ZeidelMatr(const double *, uint32_t);
#endif
