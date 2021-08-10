#include <algorithm>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <string>

#include "matrix.h"
#include "meshC.h"
#include "Mesh.h"

#define V0 1
#define ALPHA 0.25*M_PI

#define TOL 1e-12
#define MAXIT 1000000





// Computation of internal nodes discretization coefficients
void computeSteadyStateDiscretizationCoefficientsInternalNodes(const Mesh m, const double* prop,
double (*vx)(double,double), double (*vy)(double, double), double (*source)(double, double), double* A, double* b, const int scheme);


// Diagonal case functions
void computeDiscCoefsBoundaryNodesDiagonal(const Mesh m, const double* phi_boundary, double* A, double* b);
double vxDiagonal(const double, const double);
double vyDiagonal(const double, const double);
double sourceDiagonal(const double, const double);

// Smith-Hutton case functions
void computeDiscretizationCoefficientsBoundaryNodesSmithHuttonCase(const Mesh m, const double* phi_boundary, double* A, double* b);
double vxSmithHutton(const double x, const double y);
double vySmithHutton(const double x, const double y);
double sourceSmithHutton(const double, const double);

// Pre-solving check functions
void checkSystemMatrix(const int nx, const int ny, const double tol, const double* A);

// Linear system solve functions
void solveSystem(const int nx, const int ny, const double* A, const double* b, double* phi, const int method);
void solveSystemGS(const int nx, const int ny, const double tol, const int maxIt, const double* A, const double* b, double* phi);
void solveSystemLUP(const int nx, const int ny, const double* A, const double* b, double* phi);
void assembleMatrix(const int nx, const int ny, const double* A, double** AA);
void assembleMatrix(const int nx, const int ny, const double* A, double* AA);

// Print results functions
void printToFile(const Mesh m, const double* phi, const char* filename, const int precision);
void plotSolution(const char* filename);

// Check solution functions
double computeSolutionDifference(const int nx, const int ny, const double* A, const double* b, const double* phi);
void verification(const Mesh m, const double* phi);

int main(int arg, char* argv[]) {

    // // Physical data
    // Sizes
    double x0 = 0;  // Lower left corner x coordinate for rectangular domain    [m]
    double y0 = 0;  // Lower left corner y coordinate for rectangular domain    [m]
    double L = 1;   // Domain size in x and y axis                              [m]
    double lz = 1;  // Domain size in z axis                                    [m]

    // Boundary conditions
    const double phi_low = -20;      // Minimum value for phi
    const double phi_high = 20;     // Maximum value for phi

    // Thermophysical properties for water at 20 ºC
    const double lambda = 0.5861;       // Thermal conductivity                         [W/(k*m)]
    const double cv = 4183;             // Specific heat at constant volume (pressure)  [J/(kg*K)]
    const double rho = 998.2;           // Density                                      [kg/m^3]
    const double gamma = 1e-9*lambda / cv;   // Diffusion coefficient
    double* prop = (double*) malloc(2 * sizeof(double*));
    prop[0] = rho;
    prop[1] = gamma;


    // // Numerical data
    int N = 20;
    int nx = N;      // Number of nodes in x axis
    int ny = N;      // Number of nodes in y axis
    const double phi0 = 1;      // Initial value to fill phi vector for linear system resolution

    double* phi_boundary = (double*) malloc(2 * sizeof(double*));
    phi_boundary[0] = phi_low;
    phi_boundary[1] = phi_high;


    Mesh m;
    printf("Building uniform cartesian mesh...\n\n");
    m.buildUniformMesh(x0, y0, L, L, lz, nx, ny);
    if(!m.isBuilt()) {
        printf("\tError. Could not build mesh\n");
        return -2;
    }

    double* A = (double*) malloc(5 * nx * ny * sizeof(double*));
    double* b = (double*) malloc(nx * ny * sizeof(double*));
    int scheme = 0;
    if(!A) {
        printf("Not A\n");
        return -2;
    }
    if(!b) {
        printf("Not b\n");
        return -2;
    }


    printf("Computing internal nodes discretization coefficients...\n\n");
    computeSteadyStateDiscretizationCoefficientsInternalNodes(m, prop, vxDiagonal, vyDiagonal, sourceDiagonal, A, b, scheme);

    printf("Computing boundary nodes discretization coefficients for the diagonal case...\n\n");
    computeDiscCoefsBoundaryNodesDiagonal(m, phi_boundary, A, b);

    checkSystemMatrix(nx, ny, TOL, A);

    // <NEW>
    double* phi1 = (double*) malloc(nx * ny * sizeof(double*));
    std::fill_n(phi1, nx*ny, phi0);

    solveSystem(nx, ny, A, b, phi1, 0);
    // solveSystemGS(nx, ny, TOL, MAXIT, A, b, phi1);

    double* phi2 = (double*) malloc(nx * ny * sizeof(double*));
    std::fill_n(phi2, nx*ny, 0);

    solveSystem(nx, ny, A, b, phi2, 1);

    double maxDiff = 0;
    for(int i = 0; i < nx*ny; i++)
        maxDiff = std::max(maxDiff, std::abs(phi1[i] - phi2[i]));
    printf("maxDiff : %.5e\n", maxDiff);


    const char* filename = "output/output.dat";
    printToFile(m, phi1,  filename, 5);
    plotSolution(filename);

    // Free memory allocated
    free(prop);
    free(phi_boundary);
    free(A);
    free(b);
    free(phi1);
    free(phi2);

    return 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// COMPUTATION OF INTERNAL NODES DISCRETIZATION coefficients
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void computeSteadyStateDiscretizationCoefficientsInternalNodes(const Mesh m, const double* prop,
double (*vx)(double,double), double (*vy)(double, double), double (*source)(double, double), double* A, double* b, const int scheme) {
    /*
    computeSteadyStateDiscretizationCoefficientsInternalNodes: computes the discretization coefficients for the internal nodes in a steady state
    convection diffusion problem in a 2D cartesian mesh.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - m         Mesh object                                                                     [const Mesh]
        - prop      Thermophysical properties. 0: rho, 1: gamma                                     [const double*]
        - *vx       Function that gives the velocity in the x axis provided the (x,y) coordinates   [returns double, needs (double,double)]
        - *vy       Function that gives the velocity in the y axis provided the (x,y) coordinates   [returns double, needs (double,double)]
        - *source   Function that gives the source term provided the (x,y) coordinates              [returns double, needs (double,double)]
        - A         Linear system matrix set to zero. Rows: nx*ny, Columns: 5                       [double*]
        - b         Vector of indenpendent terms set to zero. Rows: nx*ny, Columns: 1               [double*]
        - scheme    Scheme used to evaluate the convective terms. (In construction yet)             [const int]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - A         Linear system matrix. Rows: nx*ny, Columns: 5                            [double*]
        - b         Vector of indenpendent terms. Rows: nx*ny, Columns: 1                    [double*]
    */

    // Initialize matrix of discretization coefficients (A) and vector of independent terms (b) to zero
    std::fill_n(A, 5*m.getNX()*m.getNY(), 0);
    std::fill_n(b, m.getNX()*m.getNY(), 0);

    // Internal nodes
    if(scheme == 0) { // Upwind-Difference Scheme
        for(int j = 1; j < m.getNY()-1; j++) {
            for(int i = 1; i < m.getNX()-1; i++) {
                int node = j * m.getNX() + i;
                double x = m.atNodeX(i);
                double y = m.atNodeY(j);
                // South node
                double mf = -prop[0] * (*vy)(x,y) * m.atSurfY(i);
                double C = (std::abs(mf) != 0 ? (mf + std::abs(mf))/(2*mf) : 0);
                double D = prop[1] * m.atSurfY(i) / m.atDistY(j-1);
                A[5*node] = D - mf * C;
                // West node
                mf = -prop[0] * (*vx)(x,y) * m.atSurfX(j);
                C = (std::abs(mf) != 0 ? (mf + std::abs(mf))/(2*mf) : 0);
                D = prop[1] * m.atSurfX(j) / m.atDistX(i-1);
                A[5*node+1] = D - mf * C;
                // East node
                mf = prop[0] * (*vx)(x,y) * m.atSurfX(j);
                C = (std::abs(mf) != 0 ? (mf - std::abs(mf))/(2*mf) : 0);
                D = prop[1] * m.atSurfX(j) / m.atDistX(i);
                A[5*node+2] = D - mf * C;
                // North node
                mf = prop[0] * (*vy)(x,y) * m.atSurfY(i);
                C = (std::abs(mf) != 0 ? (mf - std::abs(mf))/(2*mf) : 0);
                D = prop[1] * m.atSurfY(i) / m.atDistY(j);
                A[5*node+3] = D - mf * C;
                // Central node
                A[5*node+4] = A[5*node] + A[5*node+1] + A[5*node+2] + A[5*node+3] - (*source)(x,y) * m.atVol(i,j);
                // Independent term
                b[node] = (*source)(x,y) * m.atVol(i,j);
            }
        }
    } else {
        for(int j = 1; j < m.getNY()-1; j++) {
            for(int i = 1; i < m.getNX()-1; i++) {
                int node = j * m.getNX() + i;
                double x = m.atNodeX(i);
                double y = m.atNodeY(j);
                // South node
                double mf = -prop[0] * (*vy)(x,y) * m.atSurfY(i);
                double C = m.atDistNFY(2*j) / m.atDistY(j-1);
                double D = prop[1] * m.atSurfY(i) / m.atDistY(j-1);
                A[5*node] = D - mf * C;
                // West node
                mf = -prop[0] * (*vx)(x,y) * m.atSurfX(j);
                C = m.atDistNFX(2*i) / m.atDistX(i-1);
                D = prop[1] * m.atSurfX(j) / m.atDistX(i-1);
                A[5*node+1] = D - mf * C;
                // East node
                mf = prop[0] * (*vx)(x,y) * m.atSurfX(j);
                C = m.atDistNFX(2*i+1) / m.atDistX(i);
                D = prop[1] * m.atSurfX(j) / m.atDistX(i);
                A[5*node+2] = D - mf * C;
                // North node
                mf = prop[0] * (*vy)(x,y) * m.atSurfY(i);
                C = m.atDistNFY(2*j+1) / m.atDistY(j);
                D = prop[1] * m.atSurfY(i) / m.atDistY(j);
                A[5*node+3] = D - mf * C;
                // Central node
                A[5*node+4] = A[5*node] + A[5*node+1] + A[5*node+2] + A[5*node+3] - (*source)(x,y) * m.atVol(i,j);
                // Independent term
                b[node] = (*source)(x,y) * m.atVol(i,j);
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DIAGONAL CASE FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void computeDiscCoefsBoundaryNodesDiagonal(const Mesh m, const double* phi_boundary, double* A, double* b) {
    /*
    computeDiscCoefsBoundaryNodesDiagonal: computes the discretization coefficients for the boundary nodes in the diagonal case
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - m                 Mesh object                                                                     [const Mesh]
        - phi_boundary      Boundary conditions for the diagonal case. 0: phi_low, 1: phi_high              [const double*]
        - A                 Linear system matrix set to zero. Rows: nx*ny, Columns: 5                       [double*]
        - b                 Vector of indenpendent terms set to zero. Rows: nx*ny, Columns: 1               [double*]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - A                 Linear system matrix. Rows: nx*ny, Columns: 5                                   [double*]
        - b                 Vector of indenpendent terms. Rows: nx*ny, Columns: 1                           [double*]
    */
    // printf("Computing boundary nodes discretization coefficients for the diagonal case...\n\n");
    // Lower row nodes: 0 <= i <= nx-1; j=0
    for(int i = 0; i < m.getNX(); i++) {
        A[5*i+4] = 1;
        b[i] = phi_boundary[0];
    }
    // Right column nodes: i = nx-1; 1 <= j <= ny-1
    for(int j = 1; j < m.getNY(); j++) {
        int node = (j + 1) * m.getNX() - 1;
        A[5*node+4] = 1;
        b[node] = phi_boundary[0];
    }
    // Left column nodes: i = 0; 1 <= j <= ny-1
    for(int j = 1; j < m.getNY(); j++) {
        int node = j * m.getNX();
        A[5*node+4] = 1;
        b[node] = phi_boundary[1];
    }
    // Upper row nodes: 1 <= i <= nx-2, j = ny-1
    for(int i = 1; i < m.getNX()-1; i++) {
        int node = (m.getNY() - 1) * m.getNX() + i;
        A[5*node+4] = 1;
        b[node] = phi_boundary[1];
    }
}

double vxDiagonal(const double x, const double y) {
    /*
    vxDiagonal: computes the x component of velocity for the diagonal case at the point (x,y)
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - x     x coordinate    [const double]
        - y     y coordinate    [const double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - vx    x component of the velocity for the diagonal case at point (x,y)    [double]
    */
    return V0*cos(ALPHA);
}

double vyDiagonal(const double x, const double y) {
    /*
    vyDiagonal: computes the y component of velocity for the diagonal case at the point (x,y)
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - x     x coordinate    [const double]
        - y     y coordinate    [const double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - vy    y component of the velocity for the diagonal case at point (x,y)    [double]
    */
    return V0*sin(ALPHA);
}

double sourceDiagonal(const double x, const double y) {
    /*
    sourceDiagonal: computes the SP coefficient of the source term for the diagonal case at the point (x,y)
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - x         x coordinate    [const double]
        - y         y coordinate    [const double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - source    SP coefficient of the source term for the diagonal case at the point (x,y)  [double]
    */
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SMITH-HUTTON CASE FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void computeDiscretizationCoefficientsBoundaryNodesSmithHuttonCase(const Mesh m, const double* phi_boundary, double* A, double* b) {

}

double vxSmithHutton(const double x, const double y) {
    /*
    vxSmithHutton: computes the x component of velocity for the Smith-Hutton case at the point (x,y)
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - x         x coordinate    [const double]
        - y         y coordinate    [const double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - vx        x component of the velocity for the Smith-Hutton case at point (x,y)    [double]
    */
    return 2*y*(1 - x*x);
}

double vySmithHutton(const double x, const double y) {
    /*
    vxSmithHutton: computes the y component of velocity for the Smith-Hutton case at the point (x,y)
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - x         x coordinate    [const double]
        - y         y coordinate    [const double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - vy        y component of the velocity for the Smith-Hutton case at point (x,y)    [double]
    */
    return -2*x*(1 - y*y);
}

double sourceSmithHutton(const double, const double) {
    /*
    sourceSmithHutton: computes the SP coefficient of the source term for the Smith-Hutton case at the point (x,y)
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - x         x coordinate    [const double]
        - y         y coordinate    [const double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - source    SP coefficient of the source term for the Smith-Hutton case at the point (x,y)  [double]
    */
    return 0;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PRE-SOLVING CHECK FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void checkSystemMatrix(const int nx, const int ny, const double tol, const double* A) {
    /*
    checkSystemMatrix: checks
        - whether the system matrix has a row full of zeros (sufficient condition for incompatible system but not necessary)
        - whether some coefficient aP (5th column) has a zero (prevent division by 0)
    In either case the function prints a message.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - nx: discretization nodes in X axis    [const int]
        - ny: discretization nodes in Y axis    [const int]
        - tol: tolerance                        [const double]
        - A: linear system matrix               [double*]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: none
    */
    printf("Checking linear system matrix...\n");
    // Check whether the system matrix has a row full of zeros (sufficient condition for incompatible system but not necessary)
    bool foundNullRow = false;  // Bool variable: true => a null row has been found; false => no null row has been found
    int node = 0;      // Node (row) being checked
    while(node < nx*ny && !foundNullRow) {
        bool foundNonZeroElement = false;   // Bool variable: true => an element different from zero has been found in the row; false => no element different from zero has been found in the row
        int col = 0;                        // Column that is being checked
        // Check columns
        while(col < 5 && !foundNonZeroElement) {
            if(std::abs(A[5*node+col]) > tol)
                foundNonZeroElement = true;
            col++;
        }
        if(!foundNonZeroElement) {
            foundNullRow = true;
            printf("Found null row: %d\n", node);
        }
        node++;
    }

    // Check whether some coefficient aP (5th column) has a zero (prevent division by 0)
    bool nullCoef = false;  // True => a null aP coefficient has been found; false => no null aP coefficient has been foudn
    node = 0;               // Node (row) being checked
    while(node < nx*ny && !nullCoef) {
        nullCoef = (std::abs(A[5*node+4]) < tol);
        node++;
    }
    if(nullCoef) {
        printf("Found null aP coefficient, row: %d\n", node-1);
    }
    printf("\n");
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// lINEAR SYSTEM SOLVING FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void solveSystem(const int nx, const int ny, const double* A, const double* b, double* phi, const int method) {
    /*
    solveSystem: solves the linear system resulting from a 2D convection-diffusion problem in a domain discretized with a cartesian mesh using the
    specified method.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - nx        Number of nodes in the X axis                                       [const int]
        - ny        Number of nodes in the Y axis                                       [const int]
        - A         Linear system matrix. Rows: nx*ny, Columns: 5                       [double*]
        - b         Vector of indenpendent terms. Rows: nx*ny, Columns: 1               [double*]
        - phi       Solution of the linar system. Rows: nx*ny, Columns: 1               [double*]
        - metho     Method to solve the linear system. 0 - Gauss-Seidel, Other - LUP    [int]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - A                 Linear system matrix. Rows: nx*ny, Columns: 5                                   [double*]
        - b                 Vector of indenpendent terms. Rows: nx*ny, Columns: 1                           [double*]
    */
    if(method == 0)  {// Solve using Gauss-Seidel
        printf("Solving linear system using Gauss-Seidel method...\n");
        solveSystemGS(nx, ny, TOL, MAXIT, A, b, phi);
    } else {          // Solve using LUP factorization
        printf("Solving linear system using LUP factorization...\n");
        solveSystemLUP(nx, ny, A, b, phi);
    }
}

void solveSystemGS(const int nx, const int ny, const double tol, const int maxIt, const double* A, const double* b, double* phi) {
    /*
    solveSystemGS: solves the linear system resulting from a 2D convection-diffusion problem in a domain discretized with a cartesian mesh using
    Gauss-Seidel algorithm. It has two criterion to stop the iteration:
        - Let phi* and phi be two consecutive vectors of the sequence produced by Gauss-Seidel algorithm. If the infinity norm of phi-phi* is less
        than tol, then the algorithm stops.
        - If the number of iterations reach maxIt, the algorithm stops.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - nx        Number of nodes in the X axis                               [const int]
        - ny        Number of nodes in the Y axis                               [const int]
        - tol       Tolerance to stop iteration                                 [const double]
        - maxIt     Maximum number of iterations                                [const int]
        - A         Linear system matrix. Rows: nx*ny, Columns: 5               [double*]
        - b         Vector of indenpendent terms. Rows: nx*ny, Columns: 1       [double*]
        - phi       Solution of the linar system. Rows: nx*ny, Columns: 1       [double*]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - A                 Linear system matrix. Rows: nx*ny, Columns: 5                                   [double*]
        - b                 Vector of indenpendent terms. Rows: nx*ny, Columns: 1                           [double*]
    */
    // printf("Solving linear system...\n");
    int it = 0;        // Current iteration
    bool convergence = false;   // Boolean variable to tell whether there is convergence or not. False: no convergence, True: convergence
    // Gauss-Seidel iteration
    while(it < maxIt && !convergence) {
        double maxDiff = -1;    // Infinity norm of the difference phi-phi*
        // Lower row nodes
        for(int i = 0; i < nx; i++) {
            int node = i;                                                       // Node whose phi is being computed
            double aux = phi[node];                                             // Previous value of phi[node]
            phi[node] = (b[node] + A[5*node+3] * phi[node+nx]) / A[5*node+4];   // Compute new value
            maxDiff = std::max(maxDiff, std::abs(aux - phi[node]));             // Update infinity norm
        }
        // Central rows nodes
        for(int j = 1; j < ny-1; j++) {
            for(int i = 0; i < nx; i++) {
                int node = j * nx + i;                                          // Node whose phi is being computed
                double aux = phi[node];                                         // Previous value of phi[node]
                phi[node] = (b[node] + A[5*node] * phi[node-nx] + A[5*node+1] * phi[node-1] + A[5*node+2] * phi[node+1] + A[5*node+3] * phi[node+nx]) / A[5*node+4];    // Compute new value
                maxDiff = std::max(maxDiff, std::abs(aux - phi[node]));         // Update infinity norm
            }
        }
        // Upper row nodes
        for(int i = 0; i < nx; i++) {
            int node = (ny - 1) * nx + i;                                       // Node whose phi is being computed
            double aux = phi[node];                                             // Previous value of phi[node]
            phi[node] = (b[node] + A[5*node] * phi[node-nx]) / A[5*node+4];     // Compute new value
            maxDiff = std::max(maxDiff, std::abs(aux - phi[node]));             // Update infinity norm
        }
        // Final checks of the current iteration
        convergence = (maxDiff < tol);  // Convergence condition
        it++;                           // Increase iteration counter
    }
    printf("\tIterations: %d\n\n", it);
}

void solveSystemLUP(const int nx, const int ny, const double* A, const double* b, double* phi) {
    /*
    solveSystemLUP: solves the linear system resulting from a 2D convection-diffusion problem in a domain discretized with a cartesian mesh using
    the LUP factorization.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - nx        Number of nodes in the X axis                               [const int]
        - ny        Number of nodes in the Y axis                               [const int]
        - A         Linear system matrix. Rows: nx*ny, Columns: 5               [double*]
        - b         Vector of indenpendent terms. Rows: nx*ny, Columns: 1       [double*]
        - phi       Solution of the linar system. Rows: nx*ny, Columns: 1       [double*]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - phi       Solution of the linar system. Rows: nx*ny, Columns: 1       [double*]
    */

    int n = nx * ny;
    double* AA = (double*) malloc(n * n * sizeof(double*));
    // Check if memory was allocated
    if(AA) {
        std::fill_n(AA, n*n, 0);
        assembleMatrix(nx, ny, A, AA);
        // LUP method
        int* perm = (int*) malloc(n * sizeof(int*));
        printf("\tLUP factor...\n");
        int permSign = factorLU(AA, perm, n, TOL);
        if(permSign == 0) {
            printf("\tError. Matrix is singular. Setting solution to zero\n\n");
            std::fill_n(phi, n, 0);
        } else {
            printf("\tLUP solving...\n\n");
            solveLUP(AA, b, phi, perm, n);
        }
    } else {
        printf("\tError. Could not allocate memory to assemble matrix. Setting solution to zero\n\n");
        std::fill_n(phi, n, 0);
    }
    free(AA);
}

void assembleMatrix(const int nx, const int ny, const double* A, double** AA) {
    /*
    assembleMatrix: given the matrix A (is in vector form), the function assembles the matrix AA (square matrix) for LUP method.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - nx        Number of nodes in the X axis                               [const int]
        - ny        Number of nodes in the Y axis                               [const int]
        - A         Linear system matrix. Rows: nx*ny, Columns: 5               [const double*]
        - AA        Matrix to be assembled (in square matrix form)              [double**]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: temperature map of phi. Prior to running the program,
        - AA        Assembled matrix (in square matrix form)                    [double**]
    */
    printf("\tAssembling matrix...\n");
    if(AA) {
        // // Lower row
        for(int i = 0; i < nx; i++) {
            int node = i;
            AA[node][node+nx] = -A[5*node+3]; // aN coefficient
            AA[node][node] = A[5*node+4];    // aP coefficient
        }
        // // Central rows
        for(int j = 1; j < ny-1; j++) {
            // First node
            int node = j * nx;
            AA[node][node+1] = -A[5*node+2];  // aE coefficient
            AA[node][node] = A[5*node+4];    // aP coefficient
            // Central nodes
            for(int i = 1; i < nx-1; i++) {
                node = j * nx + i;
                AA[node][node-nx] = -A[5*node];   // aS coefficient
                AA[node][node-1]  = -A[5*node+1]; // aW coefficient
                AA[node][node+1]  = -A[5*node+2]; // aE coefficient
                AA[node][node+nx] = -A[5*node+3]; // aN coefficient
                AA[node][node]    = A[5*node+4]; // aP coefficient
            }
            // Last node
            node = j * nx + nx - 1;
            AA[node][node-1] = -A[5*node+1];  // aW coefficient
            AA[node][node]   = A[5*node+4];  // aP coefficient
        }
        // // Upper row
        for(int i = 0; i < nx; i++) {
            int node = (ny - 1) * nx + i;
            AA[node][node-nx] = -A[5*node];   // aS coefficient
            AA[node][node]    = A[5*node+4]; // aP coefficient
        }
    } else {
        printf("\tError. Matrix AA not allocated\n");
    }
}

void assembleMatrix(const int nx, const int ny, const double* A, double* AA) {
    /*
    assembleMatrix: given the matrix A (is in vector form), the function assembles the matrix AA (square matrix) for LUP method.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - nx        Number of nodes in the X axis                               [const int]
        - ny        Number of nodes in the Y axis                               [const int]
        - A         Linear system matrix. Rows: nx*ny, Columns: 5               [const double*]
        - AA        Matrix to be assembled (in vector form)                     [double**]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: temperature map of phi. Prior to running the program,
        - AA        Assembled matrix (in vector form)                           [double**]
    */
    printf("\tAssembling matrix...\n");
    if(AA) {
        int n = nx * ny;
        // // Lower row
        for(int i = 0; i < nx; i++) {
            AA[i*n+i]    = A[5*i+4];    // aP coefficient
            AA[i*n+i+nx] = -A[5*i+3];    // aN coefficient
        }
        // // Central rows
        for(int j = 1; j < ny-1; j++) {
            // First node
            int node = j * nx;
            AA[node*n+node] = A[5*node+4];    // aP coefficient
            AA[node*n+node+1] = -A[5*node+2];  // aE coefficient
            // Central nodes
            for(int i = 1; i < nx-1; i++) {
                node = j * nx + i;
                AA[node*n+node-nx] = -A[5*node];   // aS coefficient
                AA[node*n+node-1]  = -A[5*node+1]; // aW coefficient
                AA[node*n+node+1]  = -A[5*node+2]; // aE coefficient
                AA[node*n+node+nx] = -A[5*node+3]; // aN coefficient
                AA[node*n+node]    = A[5*node+4]; // aP coefficient
            }
            // Last node
            node = j * nx + nx - 1;
            AA[node*n+node-1] = -A[5*node+1];  // aW coefficient
            AA[node*n+node]   = A[5*node+4];  // aP coefficient
        }
        // // Upper row
        for(int i = 0; i < nx; i++) {
            int node = (ny - 1) * nx + i;
            AA[node*n+node-nx] = -A[5*node];   // aS coefficient
            AA[node*n+node]    = A[5*node+4]; // aP coefficient
        }
    } else {
        printf("\tError. Matrix AA not allocated\n");
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PRINT RESULTS FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void printToFile(const Mesh m, const double* phi, const char* filename, const int precision) {
    /*
    printToFile: prints the solution of the linear system to a file. When the file is open, the function creates a new file or overwrites the previous
    one if it already exists.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - m                 Mesh object                                                 [const Mesh]
        - phi               Solution of the linar system. Rows: nx*ny, Columns: 1       [const double*]
        - filename          Name of the file where the solution is to be written        [const char*]
        - precision         Number of decimal places to write doubles                   [const int]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: the file where the solution has been written.
    */
    printf("Printing the solution to file...\n");
    std::ofstream file;
    file.open(filename);
    if(file.is_open()) {
        printf("\tWriting to file...\n");
        file << std::setprecision(precision) << std::fixed;
        for(int i = 0; i < m.getNX(); i++) {
            for(int j = 0; j < m.getNY(); j++)
                file << m.atNodeX(i) << " " << m.atNodeY(j) << " " << phi[j*m.getNX()+i] << std::endl;
            file << std::endl;
        }
    } else
        printf("\tCould not open file\n");
    file.close();
    printf("\tClosing file...\n\n");
}

void plotSolution(const char* filename) {
    /*
    plotSolution: plots phi vs (x,y) in a temperature map using gnuplot.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - filename          Name of the file where the solution has been previously written     [const char*]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: temperature map of phi. Prior to running the program,
        - the command "export DISPLAY=:0 gnuplot" must be executed (only needed before the first execution)
        - Xming server must be running
    */

    // Plotting command, uses filename parameter
    char plotCommand[30+strlen(filename)];
    sprintf(plotCommand, "plot '%s' with image", filename);
    // Commands sent to gnuplot
    const char* GnuCommands[] = {"set xrange [0:1]", "set yrange [0:1]", "set size ratio 1", "set palette rgb 33,13,10", plotCommand};

    FILE *gnupipe = NULL;
    gnupipe = popen("gnuplot -persistent", "w");
    for(int i = 0; i < 5; i++)
        fprintf(gnupipe, "%s\n", GnuCommands[i]);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CHECK SOLUTION FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void verification(const Mesh m, const double* phi) {
    printf("Verificating solution...\n");

}

double computeSolutionDifference(const int nx, const int ny, const double* A, const double* b, const double* phi) {
    /*
    computeSolutionDifference: checks the solution of the linear system resulting from the 2D convection-diffusion equations in a cartesian mesh. If _A is
    the linear system matrix, the function returns the infinity norm of _A * phi - b.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - nx        Number of nodes in the X axis                               [const int]
        - ny        Number of nodes in the Y axis                               [const int]
        - A         Linear system matrix. Rows: nx*ny, Columns: 5               [double*]
        - b         Vector of indenpendent terms. Rows: nx*ny, Columns: 1       [double*]
        - phi       Solution of the linar system. Rows: nx*ny, Columns: 1       [double*]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - maxDiff   Infinity norm of the difference A * phi - b                 [double]
    */
    double maxDiff = -1;
    // Lower row nodes
    for(int i = 0; i < nx; i++) {
        double diff = b[i] + A[5*i+3] * phi[i+nx] - A[5*i+4] * phi[i];
        maxDiff = std::max(maxDiff, std::abs(diff));
    }
    // Central row nodes
    for(int j = 1; j < ny-1; j++) {
        for(int i = 0; i < nx; i++) {
            int node = j * nx + i;
            double diff = b[node] + (A[5*node]*phi[node-nx]) + (A[5*node+1]*phi[node-1]) + (A[5*node+2]*phi[node+1]) + (A[5*node+3]*phi[node+nx]) - (A[5*node+4]*phi[node]);
            maxDiff = std::max(maxDiff, std::abs(diff));
        }
    }
    // Upper row nodes
    for(int i = 0; i < nx; i++) {
        int node = (ny-1) * nx + i;
        double diff = b[node] + (A[5*node] * phi[node-nx]) - (A[5*node+4] * phi[node]);
        maxDiff = std::max(maxDiff, std::abs(diff));
    }
    return maxDiff;
}
