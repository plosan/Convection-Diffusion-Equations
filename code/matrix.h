#include <random>
#include <cmath>
#include <algorithm>

void printMatrix(double* mat, unsigned int rows, unsigned int cols);
void getRandomMatrix(double* mat, unsigned int rows, unsigned int cols, int min_range, int max_range);
void getSDDMatrix(double* mat, unsigned int rows, int min_range, int max_range);
void gaussSeidel(double* mat, double* vec, double* sol, unsigned int rows, double error);
double pNorm(double* vec, unsigned int size, double p);
double infNorm(double* vec, unsigned int size);
