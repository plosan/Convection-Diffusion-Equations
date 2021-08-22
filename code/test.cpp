#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <cstring>
#include <fstream>
#include <sstream>

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

struct ProblemInput {
    int type;
    int N;
    double L;
    double lz;
    double x0;
    double y0;
    double rho;
    double gamma;
    double phi_low;
    double phi_high;
};

int parseProblemInput(const char* filename, ProblemInput &p) {

    std::ifstream file;
    file.open(filename);

    printf("filename: %s\n", filename);
    if(!file.is_open()) {
        printf("Error. Could not open the file %s\n", filename);
        return -1;
    }

    // Parse case
    std::string line;
    getline(file, line);

    std::istringstream iss(line);
    std::string word;
    iss >> word;
    std::string caseDiagonal("diagonal");
    std::string caseSmithHutton("smith-hutton");
    if(caseDiagonal.compare(word) == 0)
        p.type = 1;
    else if(caseSmithHutton.compare(word) == 0)
        p.type = 2;
    else {
        printf("Error. Case provided (%s) is not valid.\n", word);
        return -1;
    }

    file >> p.N;
    if(p.N < 2) {
        printf("Error. Number of control volumes provided (%d) must be at least 2.\n", p.N);
        return -1;
    }
    if(p.type == 2 && p.N % 2 == 0) {
        printf("Warning. Number of control volumes provided (%d) for the Smith-Hutton case has been increased in one.\n", p.N);
        p.N++;
    }

    file >> p.L;
    file >> p.lz;
    file >> p.x0;
    file >> p.y0;

    file >> p.rho;
    if(p.rho <= 0) {
        printf("Error. The density provided (%.5e) must be positive.\n", p.rho);
        return -1;
    }

    file >> p.gamma;
    if(p.gamma <= 0) {
        printf("Error. The diffusion coefficient provided (%.5e) must be positive.\n", p.gamma);
        return -1;
    }

    if(p.type == 1) {
        file >> p.phi_low;
        file >> p.phi_high;
    }

    file.close();
    return 0;
}


int main(int argc, char* argv[]) {

    srand(time(NULL));

    printf("argc : %d\n", argc);

    printf("%5s%5s%s\n", "i", "", "argv[i]");
    for(int i = 0; i < argc; i++)
        printf("%5d%5s%s\n", i, "", argv[i]);

    if(argc < 2) {
        return -1;
    }

    ProblemInput p;

    parseProblemInput(argv[1], p);



    return 0;
}
