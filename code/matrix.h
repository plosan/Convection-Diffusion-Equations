#ifndef CDE_MATRIX_H_
#define CDE_MATRIX_H_

#include <random>
#include <cmath>
#include <algorithm>

void printMatrix(const int* mat, const unsigned int rows, const unsigned int cols);
void printMatrix(const double* mat, const unsigned int rows, const unsigned int cols);
void printMatrix(int** mat, const unsigned int rows, const unsigned int cols);
void printMatrix(double** mat, const unsigned int rows, const unsigned int cols);

void printReversedRowMatrix(const int* mat, const unsigned int rows, const unsigned int cols);
void printReversedRowMatrix(const double* mat, const unsigned int rows, const unsigned int cols);


void getRandomMatrix(double* mat, const unsigned int rows, const unsigned int cols, const int min_range, const int max_range);
void getRandomMatrix(int* mat, const unsigned int rows, const unsigned int cols, const int lower, const int upper);
void getRandomMatrix(double** mat, const unsigned int rows, const unsigned int cols, const int min_range, const int max_range);
void getRandomMatrix(int** mat, const unsigned int rows, const unsigned int cols, const int lower, const int upper);


void getSDDMatrix(double* mat, const unsigned int rows, const int min_range, const int max_range);
void gaussSeidel(const double* A, const double* b, double* x, const unsigned int n, const double tol, const unsigned int maxIt, int &exitCode);
void pdGaussSeidel(const double* A, const double* b, double* x, const unsigned int n, const double tol, const unsigned int maxIt, int &exitCode);

double pNorm(double* vec, unsigned int size, double p);
double infNorm(double* vec, unsigned int size);

void swapRows(double** A, int* perm, int i, int k, int& permSign);
int factorLU(double** A, int* perm, const int n, const double tol);
void solveLUP(double** A, double* b, double* x, int* perm, const int n);

void swapRows(double* A, const int n, int* perm, int i, int k, int& permSign);
int factorLU(double* A, int* perm, const int n, const double tol);
void solveLUP(double* A, double* b, double* x, int* perm, const int n);

#endif //CDE_MATRIX_H_
