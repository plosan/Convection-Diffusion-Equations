#ifndef CDE_MATRIX_H_
#define CDE_MATRIX_H_

#include <random>
#include <cmath>
#include <algorithm>

void printMatrix(const int* mat, const int rows, const int cols);
void printMatrix(const double* mat, const int rows, const int cols);
void printMatrix(int** mat, const int rows, const int cols);
void printMatrix(double** mat, const int rows, const int cols);
void printNonZeroElements(const double* mat, const int rows, const int cols, const double tol);

void printReversedRowMatrix(const int* mat, const int rows, const int cols);
void printReversedRowMatrix(const double* mat, const int rows, const int cols);


void getRandomMatrix(double* mat, const int rows, const int cols, const int min_range, const int max_range);
void getRandomMatrix(int* mat, const int rows, const int cols, const int lower, const int upper);
void getRandomMatrix(double** mat, const int rows, const int cols, const int min_range, const int max_range);
void getRandomMatrix(int** mat, const int rows, const int cols, const int lower, const int upper);


void getSDDMatrix(double* mat, const int rows, const int min_range, const int max_range);
void gaussSeidel(const double* A, const double* b, double* x, const int n, const double tol, const int maxIt, int &exitCode);
void pdGaussSeidel(const double* A, const double* b, double* x, const int n, const double tol, const int maxIt, int &exitCode);

double pNorm(double* vec, int size, double p);
double infNorm(double* vec, int size);

void swapRows(double** A, int* perm, int i, int k, int& permSign);
int factorLU(double** A, int* perm, const int n, const double tol);
void solveLUP(double** A, const double* b, double* x, int* perm, const int n);

void swapRows(double* A, const int n, int* perm, int i, int k, int& permSign);
int factorLU(double* A, int* perm, const int n, const double tol);
void solveLUP(const double* A, const double* b, double* x, int* perm, const int n);

#endif //CDE_MATRIX_H_
