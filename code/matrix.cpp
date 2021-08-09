#include "matrix.h"

void printMatrix(const double* mat, const int rows, const int cols) {
    /*
    printMatrix: prints a double array. The array must be given in "vector" format.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - mat: matrix to be printed     [const double*]
        - rows: matrix rows             [const int]
        - cols: matrix columns          [const double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: none
    */

    if(mat) {
        for(int i = 0; i < rows; ++i) {
            for(int j = 0; j < cols; ++j)
                printf("%10.4f", mat[i*cols+j]);
            printf("\n");
        }
        printf("\n");
    }
}

void printMatrix(const int* mat, const int rows, const int cols) {
    /*
    printMatrix: prints a double array. The array must be given in "vector" format.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - mat: matrix to be printed     [const double*]
        - rows: matrix rows             [const int]
        - cols: matrix columns          [const double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: none
    */

    if(mat) {
        for(int i = 0; i < rows; ++i) {
            for(int j = 0; j < cols; ++j)
                printf("%10d", mat[i*cols+j]);
            printf("\n");
        }
        printf("\n");
    }
}

void printMatrix(double** mat, const int rows, const int cols) {
    /*
    printMatrix: prints a double array. The array must be given in "vector" format.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - mat: matrix to be printed     [const double*]
        - rows: matrix rows             [const int]
        - cols: matrix columns          [const double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: none
    */

    if(mat) {
        for(int i = 0; i < rows; ++i) {
            for(int j = 0; j < cols; ++j)
                printf("%10.4f", mat[i][j]);
            printf("\n");
        }
        printf("\n");
    }
}

void printMatrix(int** mat, const int rows, const int cols) {
    /*
    printMatrix: prints a double array. The array must be given in "vector" format.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - mat: matrix to be printed     [const double*]
        - rows: matrix rows             [const int]
        - cols: matrix columns          [const double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: none
    */

    if(mat) {
        for(int i = 0; i < rows; ++i) {
            for(int j = 0; j < cols; ++j)
                printf("%10d", mat[i][j]);
            printf("\n");
        }
        printf("\n");
    }
}

void printReversedRowMatrix(const int* mat, const int rows, const int cols) {
    /*
    printReversedRowMatrix: prints a double array with rows in reversed order, that is, from last row to first.
    The array must be given in "vector" format.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - mat: matrix to be printed     [const double*]
        - rows: matrix rows             [const int]
        - cols: matrix columns          [const double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: none
    */
    if(mat) {
        for(int i = rows-1; i > -1; --i) {
            for(int j = 0; j < cols; ++j)
                printf("%10d", mat[i*cols+j]);
            printf("\n");
        }
        printf("\n");
    }
}

void printReversedRowMatrix(const double* mat, const int rows, const int cols) {
    /*
    printReversedRowMatrix: prints a double array with rows in reversed order, that is, from last row to first.
    The array must be given in "vector" format.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - mat: matrix to be printed     [const double*]
        - rows: matrix rows             [const int]
        - cols: matrix columns          [const double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: none
    */
    if(mat) {
        for(int i = rows-1; i > -1; --i) {
            for(int j = 0; j < cols; ++j)
                printf("%10.4f", mat[i*cols+j]);
            printf("\n");
        }
        printf("\n");
    }
}

void getRandomMatrix(double* mat, const int rows, const int cols, const int lower, const int upper) {
    /*
    getRandomMatrix: returns a double array filled with random integers in the interval [lower, upper]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - mat: matrix to be filled                      [double*]
        - rows: matrix rows                             [const int]
        - cols: matrix columns                          [const int]
        - lower: lower bound for the random integers    [const int]
        - upper: upper bound for the random integers    [const int]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - mat: double array filled with random integers in the interval [lower, upper]  [double*]
    */

    if(mat) {
        for(int i = 0; i < rows; ++i)
            for(int j = 0; j < cols; ++j)
                mat[i*cols+j] = rand()%(upper - lower + 1) + lower;
    }
}

void getRandomMatrix(int* mat, const int rows, const int cols, const int lower, const int upper) {
    /*
    getRandomMatrix: returns a double array filled with random integers in the interval [lower, upper]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - mat: matrix to be filled                      [double*]
        - rows: matrix rows                             [const int]
        - cols: matrix columns                          [const int]
        - lower: lower bound for the random integers    [const int]
        - upper: upper bound for the random integers    [const int]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - mat: double array filled with random integers in the interval [lower, upper]  [double*]
    */

    if(mat) {
        for(int i = 0; i < rows; ++i)
            for(int j = 0; j < cols; ++j)
                mat[i*cols+j] = rand()%(upper - lower + 1) + lower;
    }
}

void getRandomMatrix(double** mat, const int rows, const int cols, const int lower, const int upper) {
    /*
    getRandomMatrix: returns a double array filled with random integers in the interval [lower, upper]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - mat: matrix to be filled                      [double*]
        - rows: matrix rows                             [const int]
        - cols: matrix columns                          [const int]
        - lower: lower bound for the random integers    [const int]
        - upper: upper bound for the random integers    [const int]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - mat: double array filled with random integers in the interval [lower, upper]  [double*]
    */

    if(mat) {
        for(int i = 0; i < rows; ++i)
            for(int j = 0; j < cols; ++j)
                mat[i][j] = rand()%(upper - lower + 1) + lower;
    }
}

void getRandomMatrix(int** mat, const int rows, const int cols, const int lower, const int upper) {
    /*
    getRandomMatrix: returns a double array filled with random integers in the interval [lower, upper]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - mat: matrix to be filled                      [double*]
        - rows: matrix rows                             [const int]
        - cols: matrix columns                          [const int]
        - lower: lower bound for the random integers    [const int]
        - upper: upper bound for the random integers    [const int]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - mat: double array filled with random integers in the interval [lower, upper]  [double*]
    */

    if(mat) {
        for(int i = 0; i < rows; ++i)
            for(int j = 0; j < cols; ++j)
                mat[i][j] = rand()%(upper - lower + 1) + lower;
    }
}

void getSDDMatrix(double* mat, const int rows, const int lower, const int upper) {
    /*
    getSDDMatrix: returns a double array filled with random integers in the interval [lower, upper] which make the matrix be a strictly diagonally
    dominant (SDD) matrix. This is useful to the Gauss-Seidel method to solve linear systems.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - mat: matrix to be filled                      [double*]
        - rows: matrix rows                             [const int]
        - lower: lower bound for the random integers    [const int]
        - upper: upper bound for the random integers    [const int]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - mat: double array filled with random integers in the interval [lower, upper] which turn it into a SDD matrix  [double*]
    */
    if(mat) {
        for(int i = 0; i < rows; ++i) {
            int sum = 0;
            for(int j = 0; j < rows; ++j) {
                mat[i*rows+j] = rand()%(upper - lower + 1) + lower;
                sum += abs(mat[i*rows+j]);
            }
            mat[i*(rows+1)] = (rand()%2 == 0 ? 1 : -1)*sum;
        }
    }
}

void gaussSeidel(const double* A, const double* b, double* x, const int n, const double tol, const int maxIt, int &exitCode) {
    /*
    gaussSeidel: solves the linear system with matrix A and vector b using Gauss-Seidel method
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - A: linear system matrix                           [const double*]
        - b: linear system vector                           [const double*]
        - x: initial approximation to the solution          [double*]
        - n: linear system dimension                        [const int]
        - tol: tolerance criterion to stop the iteration.   [const double]
        - maxIt: max number of iterations to be made.       [const int]
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
            for(int i = 0; i < n; i++) {
                double aux = x[i];  // Previous value
                // Compute x[i] using Gauss-Seidel
                x[i] = b[i];
                for(int j = 0; j < i; j++)
                    x[i] -= A[i*n+j] * x[j];
                for(int j = i + 1; j < n; j++)
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

double pNorm(double* vec, int size, double p) {
    if(!vec)
        return -1;

    if(p < 1)
        return -1;

    double norm = 0;
    for(int i = 0; i < size; ++i)
        norm += pow(abs(vec[i]), p);
    return pow(norm, 1/p);
}

double infNorm(double* vec, int size) {
    if(!vec)
        return -1;

    double norm = -1;
    for(int i = 0; i < size; ++i)
        norm = (norm > abs(vec[i]) ? norm : abs(vec[i]));
    return norm;
}



void swapRows(double** A, int* perm, int i, int k, int& permSign) {
    /*
    swapRows: swaps rows i and k of the matrix A, writes the permutation on perm and changes the permutation sign
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - A         Matrix                                                              [double**]
        - perm      Permutation vector. If perm[a] = b, it means row b goes to row a    [int*]
        - i         First row to be swapped                                             [int]
        - k         Second row to be swapped                                            [int]
        - permSign  Sign of the permutation                                             [int&]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - A         Matrix with swapped rows                        [double**]
        - perm      Permutation vector with the new permutation     [int*]
        - permSign  Permutation sign times -1                       [int&]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    */
    // Swap rows
    double* aux = A[i];
    A[i] = A[k];
    A[k] = aux;
    // Change permutation sign
    permSign *= -1;
    // Write permutation
    int z = perm[i];
    perm[i] = perm[k];
    perm[k] = z;
}

int factorLU(double** A, int* perm, const int n, const double tol) {
    /*
    lup: performs the LUP factorization with scaled partial pivoting. Matrices L and U are stored in the original matrix A.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - A         Matrix                                                              [double**]
        - perm      Permutation vector. If perm[a] = b, it means row b goes to row a    [int*]
        - permSign  Sign of the permutation                                             [int&]
        - n         Matrix dimension                                                    [const int int]
        - tol       Tolerance to decide numerical zero                                  [cont double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - A         Matrix with L below the diagonal and U on the diagonal and above    [double**]
        - permSign  Sign of the permutation                                             [int&]
        - perm      Resulting permutation vector of the factorization                   [int*]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    If A is singular, its issues an error message and permSign is set to 0.
    */
    int permSign = 1;
    // Initialize permutation vector
    for(int i = 0; i < n; ++i)
        perm[i] = i;
    // LUP factorization
    for(int k = 0; k < n; ++k) {
        // Scaled partial pivoting
        double bestCoef = 0;
        int bestRow = 0;
        for(int i = k; i < n; i++) {
            double s = 0;
            for(int j = k; j < n; j++)
                s = std::max(s, std::abs(A[i][j]));
            double coef = std::abs(A[i][k] / s);
            if(coef > bestCoef) {
                bestCoef = coef;
                bestRow = i;
            }
        }
        // Row swap
        if(k != bestRow)
            swapRows(A, perm, k, bestRow, permSign);
        // Check if the pivot is non null
        if(std::abs(A[k][k]) < tol) {
            printf("\tError\n");
            return 0;
        }
        // Gauss' elimination
        for(int i = k+1; i < n; i++) {
            double m = A[i][k] / A[k][k];
            for(int j = k; j < n; j++)
                A[i][j] -= m * A[k][j];
            A[i][k] = m;
        }
    }
    return permSign;
}

void solveLUP(double** A, const double* b, double* x, int* perm, const int n) {
    /*
    solveLUP: solves the linear system given by matrix A and vector b. Matrix A has to be already in its LU form, that is to say, L below the
    diagonal and U on the diagonal and above.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - A         Linear system matrix already in LU form     [double**]
        - b         Linear system vector of independent terms   [double*]
        - x         Linear system solution                      [double*]
        - perm      Permutation vector                          [int*]
        - n         Matrix dimension                            [const int]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - x         Linear system solution                      [double*]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    */
    // Solve L y = P b
    double* y = (double*) malloc(n * sizeof(double*));
    for(int i = 0; i < n; i++) {
        double sum = 0;
        for(int j = 0; j < i; j++)
            sum += A[i][j] * y[j];
        y[i] = b[perm[i]] - sum;
    }
    // Solve U x = y
    for(int i = n-1; i >= 0; i--) {
        double sum = 0;
        for(int j = i+1; j < n; j++)
            sum += A[i][j] * x[j];
        x[i] = (y[i] - sum) / A[i][i];
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void swapRows(double* A, const int n, int* perm, int i, int k, int& permSign) {
    /*
    swapRows: swaps rows i and k of the matrix A, writes the permutation on perm and changes the permutation sign
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - A         Matrix given in vector form                                         [double**]
        - n         Matrix dimension                                                    [const int]
        - perm      Permutation vector. If perm[a] = b, it means row b goes to row a    [int*]
        - i         First row to be swapped                                             [int]
        - k         Second row to be swapped                                            [int]
        - permSign  Sign of the permutation                                             [int&]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - A         Matrix with swapped rows                        [double**]
        - perm      Permutation vector with the new permutation     [int*]
        - permSign  Permutation sign times -1                       [int&]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    */
    // Swap matrix rows
    for(int j = 0; j < n; j++) {
        double aux = A[i*n+j];
        A[i*n+j] = A[k*n+j];
        A[k*n+j] = aux;
    }
    // Change permutation sign
    permSign *= -1;
    // Write permutation
    int z = perm[i];
    perm[i] = perm[k];
    perm[k] = z;
}

int factorLU(double* A, int* perm, const int n, const double tol) {
    /*
    lup: performs the LUP factorization of matrix A with scaled partial pivoting. Matrices L and U are stored in the original matrix A. Matrix A is
    in diagonal form.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - A         Matrix in diagonal form                                             [double**]
        - perm      Permutation vector. If perm[a] = b, it means row b goes to row a    [int*]
        - n         Matrix dimension                                                    [const int int]
        - tol       Tolerance to decide numerical zero                                  [cont double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - A         Matrix with L below the diagonal and U on the diagonal and above    [double**]
        - permSign  Sign of the permutation                                             [int&]
        - perm      Resulting permutation vector of the factorization                   [int*]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    If A is singular, its issues an error message and permSign is set to 0.
    */
    int permSign = 1;
    // Initialize permutation vector
    for(int i = 0; i < n; ++i)
        perm[i] = i;
    // LUP factorization
    for(int k = 0; k < n; ++k) {
        // Scaled partial pivoting
        double bestCoef = 0;
        int bestRow = 0;
        for(int i = k; i < n; i++) {
            double s = 0;
            for(int j = k; j < n; j++)
                s = std::max(s, std::abs(A[i*n+j]));
            double coef = std::abs(A[i*n+k] / s);
            if(coef > bestCoef) {
                bestCoef = coef;
                bestRow = i;
            }
        }
        // Row swap
        if(k != bestRow)
            swapRows(A, n, perm, k, bestRow, permSign);
        // Check if the pivot is non null
        if(std::abs(A[k*n+k]) < tol) {
            printf("\tError\n");
            return 0;
        }
        // Gauss' elimination
        for(int i = k+1; i < n; i++) {
            double m = A[i*n+k] / A[k*n+k];
            for(int j = k; j < n; j++)
                A[i*n+j] -= m * A[k*n+j];
            A[i*n+k] = m;
        }
    }
    return permSign;
}

void solveLUP(const double* A, const double* b, double* x, int* perm, const int n) {
    /*
    solveLUP: solves the linear system given by matrix A and vector b. Matrix A has to be already in its LU form, that is to say, L below the
    diagonal and U on the diagonal and above. In addition, matrix A  (in its LU form) is given in vector form.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - A         Linear system matrix already in LU form (and vector form)   [double**]
        - b         Linear system vector of independent terms                   [double*]
        - x         Linear system solution                                      [double*]
        - perm      Permutation vector                                          [int*]
        - n         Matrix dimension                                            [const int]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - x         Linear system solution                                      [double*]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    */
    // Solve L y = P b
    double* y = (double*) malloc(n * sizeof(double*));
    for(int i = 0; i < n; i++) {
        double sum = 0;
        for(int j = 0; j < i; j++)
            sum += A[i*n+j] * y[j];
        y[i] = b[perm[i]] - sum;
    }
    // Solve U x = y
    for(int i = n-1; i >= 0; i--) {
        double sum = 0;
        for(int j = i+1; j < n; j++)
            sum += A[i*n+j] * x[j];
        x[i] = (y[i] - sum) / A[i*n+i];
    }
}
