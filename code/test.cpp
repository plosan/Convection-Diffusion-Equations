#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include "matrix.h"

void printIntegerMatrix(const int* mat, const unsigned int rows, const unsigned int cols) {
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
                printf("%10d", mat[i*cols+j]);
            printf("\n");
        }
        printf("\n");
    }
}

void getRandomIntegerMatrix(int* mat, const unsigned int rows, const unsigned int cols, const int lower, const int upper) {
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

void gaussSeidel(const double* A, const double* b, double* x, const unsigned int n, const double tol, const unsigned int maxIt, int &exitCode) {
    /*
    gaussSeidel: solves the linear system with matrix A and vector b using Gauss-Seidel method
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - A: linear system matrix                           [const double*]
        - b: linear system vector                           [const double*]
        - x: initial approximation to the solution          [double*]
        - n: linear system dimension                        [const unsigned int]
        - tol: tolerance criterion to stop the iteration.   [const double]
        - maxIt: max number of iterations to be made.       [const unsigned int]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - x: solution to the linear system or last iteration.   [double*]
        - exitCode: code expressing how the algorithm finished. [int&]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Explanation: Gauss-Seidel's algorithm is an iterative algorithm to compute the solution of a linear system. From an initial solution, it
    computes a sequence of vectors which ideally converge to the actual solution of the system. If the matrix is strictly diagonally dominant, then
    the algorithm is guaranteed to converge. Moreover, it needs non zero elements in the matrix diagonal in order to work properly.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    This implementation does not check whether the sequence of vectors diverges or not.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    exitCode:
        - It is initialized to -maxIt.
        - As iterations happen, its value is increased. If at some point exitCode == 0, the iteration is halted. Then the exitCode is 0.
        - Each iteration the convergence condition is checked. If x, y in R^n are two consecutive vectors of the sequence, the convergence checked
            consists on calculating the maximum of the differences abs(x[i]-y[i]). In the case that the maximum difference is less than the tolerance,
            the convergence has been achieved. The exitCode in this case is 1.
        - If the algorithm cannot be started as a consequence of A, b or x being null pointers, the exitCode is 2.
    */

    if(A && b && x) { // Check if pointers are non-null
        exitCode = -maxIt; // exitCode following the aforementioned
        // Covergence of the solution and the iteration can be controlled with the exitCode. Just initialize exitCode to some value and execute the
        // while loop as long as it is that value. Change the value if convergence or divergence is detected.
        while(exitCode < 0) {
            double maxDiff = -1;    // Maximum of the differences abs(x[i]-y[i]), where x and y are two consecutive elements in the sequence of vectors
            for(unsigned int i = 0; i < n; i++) {
                double aux = x[i];  // Previous value
                // Compute x[i] using Gauss-Seidel
                x[i] = b[i];
                for(unsigned int j = 0; j < i; j++)
                    x[i] -= A[i*n+j] * x[j];
                for(unsigned int j = i + 1; j < n; j++)
                    x[i] -= A[i*n+j] * x[j];
                x[i] /= A[i*n+i];
                // Compute difference and store the maximum
                double diff = std::abs(x[i] - aux);
                maxDiff = (maxDiff > diff ? maxDiff : diff);
            }
            exitCode++; // Increase iteration counter
            // Check convergence
            if(maxDiff < tol)
                exitCode = 1;
        }
    } else // Error when some pointer is not initialized
        exitCode = 2;
}

int main(void) {

    const unsigned int n = 5;
    double* A = (double*) malloc(n * n * sizeof(double*));
    double* b = (double*) malloc(n * sizeof(double*));
    double* x = (double*) malloc(n * sizeof(double*));
    const double tol = 1e-15;
    const unsigned int maxIt = 500;

    getSDDMatrix(A, n, -20, 20);
    getRandomMatrix(b, n, 1, -20, 20);

    printf("A = \n");
    printMatrix(A, n, n);

    printf("b = \n");
    printMatrix(b, n, 1);

    int exitCode = 0;

    gaussSeidel(A, b, x, n, tol, maxIt, exitCode);

    printf("exitCode = %d\n\n", exitCode);

    printf("A = \n");
    printMatrix(A, n, n);

    printf("b = \n");
    printMatrix(b, n, 1);

    printf("x = \n");
    printMatrix(x, n, 1);



    double* prod = (double*) malloc(n * sizeof(double*));
    for(int i = 0; i < n; i++) {
        prod[i] = 0;
        for(int j = 0; j < n; j++)
            prod[i] += A[i*n+j] * x[j];
    }

    double maxDiff = -1;
    for(int i = 0; i < n; i++) {
        double diff = std::abs(prod[i] - b[i]);
        maxDiff = (maxDiff > diff ? maxDiff : diff);
    }

    printf("%10s : %10.5e\n", "maxDiff", maxDiff);

    free(A);
    free(b);
    free(x);
    free(prod);
}
