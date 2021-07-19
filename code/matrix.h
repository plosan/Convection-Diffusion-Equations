#include <random>
#include <cmath>
#include <algorithm>

void printMatrix(const double* mat, const unsigned int rows, const unsigned int cols);
void printReversedRowMatrix(const double* mat, const unsigned int rows, const unsigned int cols);
void getRandomMatrix(double* mat, const unsigned int rows, const unsigned int cols, const int min_range, const int max_range);
void getSDDMatrix(double* mat, const unsigned int rows, const int min_range, const int max_range);
void gaussSeidel(const double* mat, const double* vec, double* sol, const unsigned int rows, double* error);

double pNorm(double* vec, unsigned int size, double p);
double infNorm(double* vec, unsigned int size);
