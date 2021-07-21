#include "matrix.h"

void printMatrix(const double* mat, const unsigned int rows, const unsigned int cols) {
    /*
    printMatrix: prints a double array. The array must be given in "vector" format.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - mat: matrix to be printed     [const double*]
        - rows: matrix rows             [const unsigned int]
        - cols: matrix columns          [const double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: none
    */

    if(mat) {
        for(unsigned int i = 0; i < rows; ++i) {
            for(unsigned int j = 0; j < cols; ++j)
                printf("%10.4f", mat[i*cols+j]);
            printf("\n");
        }
        printf("\n");
    }
}

void printReversedRowMatrix(const double* mat, const unsigned int rows, const unsigned int cols) {
    /*
    printReversedRowMatrix: prints a double array with rows in reversed order, that is, from last row to first.
    The array must be given in "vector" format.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - mat: matrix to be printed     [const double*]
        - rows: matrix rows             [const unsigned int]
        - cols: matrix columns          [const double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: none
    */

    if(mat) {
        for(int i = rows-1; i > -1; --i) {
            for(unsigned int j = 0; j < cols; ++j)
                printf("%10.4f", mat[i*cols+j]);
            printf("\n");
        }
        printf("\n");
    }
}


void getRandomMatrix(double* mat, const unsigned int rows, const unsigned int cols, const int lower, const int upper) {
    /*
    getRandomMatrix: returns a double array filled with random integers in the interval [lower, upper]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - mat: matrix to be filled                      [double*]
        - rows: matrix rows                             [const unsigned int]
        - cols: matrix columns                          [const unsigned int]
        - lower: lower bound for the random integers    [const int]
        - upper: upper bound for the random integers    [const int]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - mat: double array filled with random integers in the interval [lower, upper]  [double*]
    */

    if(mat) {
        for(unsigned int i = 0; i < rows; ++i)
            for(unsigned int j = 0; j < cols; ++j)
                mat[i*cols+j] = rand()%(upper - lower + 1) + lower;
    }
}


void getSDDMatrix(double* mat, const unsigned int rows, const int lower, const int upper) {
    /*
    getSDDMatrix: returns a double array filled with random integers in the interval [lower, upper] which make the matrix be a strictly diagonally
    dominant (SDD) matrix. This is useful to the Gauss-Seidel method to solve linear systems.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - mat: matrix to be filled                      [double*]
        - rows: matrix rows                             [const unsigned int]
        - lower: lower bound for the random integers    [const int]
        - upper: upper bound for the random integers    [const int]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - mat: double array filled with random integers in the interval [lower, upper] which turn it into a SDD matrix  [double*]
    */
    if(mat) {
        for(unsigned int i = 0; i < rows; ++i) {
            int sum = 0;
            for(unsigned int j = 0; j < rows; ++j) {
                mat[i*rows+j] = rand()%(upper - lower + 1) + lower;
                sum += abs(mat[i*rows+j]);
            }
            mat[i*(rows+1)] = (rand()%2 == 0 ? 1 : -1)*sum;
        }
    }
}


void gaussSeidel(const double* A, const double* b, double* x, const unsigned int n, const double tol, int* exitCode) {
    // if(mat && vec && sol) {
    //
    // } else {
    //     error = -1;
    // }
}

double pNorm(double* vec, unsigned int size, double p) {
    if(!vec)
        return -1;

    if(p < 1)
        return -1;

    double norm = 0;
    for(unsigned int i = 0; i < size; ++i)
        norm += pow(abs(vec[i]), p);
    return pow(norm, 1/p);
}

double infNorm(double* vec, unsigned int size) {
    if(!vec)
        return -1;

    double norm = -1;
    for(unsigned int i = 0; i < size; ++i)
        norm = (norm > abs(vec[i]) ? norm : abs(vec[i]));
    return norm;
}
