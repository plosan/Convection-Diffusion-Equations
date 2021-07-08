#include <iostream>
#include <cstdlib>
#include <ctime>
#include <random>
#include <string>
#include "matrix.h"

void computeGeometry();

int main(int arg, char* argv[]) {

    srand((unsigned) time(0));

    int rows = 8;
    int cols = 4;

    double* mat = (double*) malloc(rows * cols * sizeof(double*));
    getRandomMatrix(mat, rows, cols, -50, 50);
    printMatrix(mat, rows, cols);
    free(mat);

    double* vec = (double*) malloc(rows*sizeof(double*));
    getRandomMatrix(vec, rows, 1, -50, 50);
    printMatrix(vec, rows, 1);

    printf("%.5f\n", infNorm(vec, rows));


    double p = 2.45678;
    printf("%.5f\n", pNorm(vec, rows, p));

    return 0;
}
