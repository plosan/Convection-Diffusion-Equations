#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <cstring>

#include "Mesh.h"
#include "matrix.h"

#define TOL 1e-12

int n;
int* mat;

int at(int i) {
    return mat[i];
}

int add(int x, int y) {
    return x + y;
}

int prod(int x, int y) {
    return x*y;
}

void printFunction(int x, int y, int (*f)(int, int), int (*g)(int,int)) {
    int z = (*f)(x,y);
    int p = (*g)(x,y);
    printf("%d + %d = %d\n", x, y, z);
    printf("%d * %d = %d\n", x, y, p);
}

void testLUP1(int min_n, int max_n, int max_m) {
    for(int n = min_n; n < max_n; n++) {
        for(int m = 0; m < max_m; m++) {
            // Linear system matrix
            double** A = (double**) malloc(n * sizeof(double*));
            for(int i = 0; i < n; i++)
                A[i] = (double*) malloc(n * sizeof(double*));
            getRandomMatrix(A, n, n, -20, 20);
            // Linear system matrix backup
            double** B = (double**) malloc(n * sizeof(double*));
            for(int i = 0; i < n; i++)
                B[i] = (double*) malloc(n * sizeof(double*));
            for(int i = 0; i < n; i++)
                memcpy(B[i], A[i], n*sizeof(double*));
            // Permutation vector
            int* perm = (int*) malloc(n * sizeof(int*));
            // LUP factorization of A
            factorLU(A, perm, n, TOL);
            // Vector of independent terms
            double* b = (double*) malloc(n * sizeof(double*));
            getRandomMatrix(b, n, 1, -20, 20);
            // Solve linear system
            double* x = (double*) malloc(n * sizeof(double*));
            solveLUP(A, b, x, perm, n);
            // Compute infinity norm of A x - b
            double maxDiff = 0;
            for(int i = 0; i < n; i++) {
                double sum = 0;
                for(int j = 0; j < n; j++)
                    sum += B[i][j] * x[j];
                maxDiff = std::max(maxDiff, std::abs(sum - b[i]));
            }
            printf("%10d %10d %20.5e %s\n", n, m, maxDiff, (maxDiff < 1e-9 ? "" : "Error"));
        }
    }
}

void testLUP2(int min_n, int max_n, int max_m) {
    for(int n = min_n; n < max_n; n++) {
        for(int m = 0; m < max_m; m++) {
            // Linear system matrix
            double* A = (double*) malloc(n * n * sizeof(double*));
            getRandomMatrix(A, n, n, -20, 20);
            // Linear system matrix backup
            double* B = (double*) malloc(n * n * sizeof(double*));
            memcpy(B, A, n*n*sizeof(double*));
            // Permutation vector
            int* perm = (int*) malloc(n * sizeof(int*));
            // LUP factorization of A
            factorLU(A, perm, n, TOL);
            // Vector of independent terms
            double* b = (double*) malloc(n * sizeof(double*));
            getRandomMatrix(b, n, 1, -20, 20);
            // Solve linear system
            double* x = (double*) malloc(n * sizeof(double*));
            solveLUP(A, b, x, perm, n);
            // Compute infinity norm of A x - b
            double maxDiff = 0;
            for(int i = 0; i < n; i++) {
                double sum = 0;
                for(int j = 0; j < n; j++)
                    sum += B[i*n+j] * x[j];
                maxDiff = std::max(maxDiff, std::abs(sum - b[i]));
            }
            printf("%10d %10d %20.5e %s\n", n, m, maxDiff, (maxDiff < 1e-9 ? "" : "Error"));
        }
    }
}

int main(void) {

    srand(time(NULL));

    testLUP2(10, 200, 10);

    return 0;
}
