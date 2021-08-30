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

#define TOL 1e-12
#define MAXIT 1000000

#define SCHEME_UDS 0
#define SCHEME_CDS 1
#define SCHEME_EDS 2
#define SCHEME_HYBRID 3
#define SCHEME_POWERLAW 4

// Computation of internal nodes discretization coefficients
double schemeUDS(const double P);
double schemeCDS(const double P);
double schemeHybrid(const double P);
double schemePowerlaw(const double P);
double schemeEDS(const double P);

void computeSteadyStateDiscretizationCoefficientsInternalNodes(const Mesh m, const double rho, const double gamma,
double (*vx)(double,double), double (*vy)(double, double), double (*sourceP)(double, double),
double (*sourceC)(double, double), double (*scheme)(double), double* A, double* b);

void computeSteadyStateDiscretizationCoefficientsInternalNodes(const Mesh m, const double rho, const double gamma,
double (*vx)(double,double), double (*vy)(double, double), double (*sourceP)(double, double),
double (*sourceC)(double, double), int scheme, double* A, double* b);

// Diagonal case functions
int solveDiagonalCase(const double L, const double lz, const int N, const double rho, const double gamma, const double phi_low, const double phi_high);
void computeDiscCoefsBoundaryNodesDiagonal(const Mesh m, const double phi_low, const double phi_high, double* A, double* b);
double vxDiagonal(const double, const double);
double vyDiagonal(const double, const double);
double sourcePDiagonal(const double, const double);
double sourceCDiagonal(const double, const double);

// Smith-Hutton case functions
int solveSmithHuttonCase(const double L, const double lz, const int N, const double rho, const double gamma);
void computeDiscretizationCoefficientsBoundaryNodesSmithHuttonCase(const Mesh m, double* A, double* b);
double vxSmithHutton(const double x, const double y);
double vySmithHutton(const double x, const double y);
double sourcePSmithHutton(const double, const double);
double sourceCSmithHutton(const double, const double);

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
void plotSolution(const Mesh m, const char* filename);
void printVelocityField(const Mesh m, double (*vx)(double, double), double (*vy)(double, double), const char* filename_mod, const char* filename_vec, const int precision);

// Check solution functions
double checkLinearSystemSolution(const int nx, const int ny, const double* A, const double* b, const double* phi);

int main(int argc, char* argv[]) {

    // Physical data
    double L = 1;   // Domain size in x and y axis  [m]
    double lz = 1;  // Domain size in z axis        [m]

    // Numerical data
    int N = 200;    // Nodes for uniform discretization

    // Thermophysical properties
    const double rho = 1000;        // Density                  [kg/m^3]
    const double gamma = 1e12*rho;   // Diffusion coefficient

    // solveDiagonalCase(L, lz, N, rho, gamma, 0, 1);

    solveSmithHuttonCase(L, lz, N, rho, gamma);

    return 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// COMPUTATION OF INTERNAL NODES DISCRETIZATION COEFFICIENTS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double schemeUDS(const double P) {
    return 1;
}

double schemeCDS(const double P) {
    return 1 - 0.5 * std::abs(P);
}

double schemeHybrid(const double P) {
    return std::max(0.0, 1 - 0.5 * std::abs(P));
}

double schemePowerlaw(const double P) {
    return std::max(0.0, std::pow(1 - 0.1*std::abs(P), 5));
}

double schemeEDS(const double P) {
    return std::abs(P)/(std::exp(std::abs(P)) - 1);
}

void computeSteadyStateDiscretizationCoefficientsInternalNodes(const Mesh m, const double rho, const double gamma,
double (*vx)(double,double), double (*vy)(double, double), double (*sourceP)(double, double),
double (*sourceC)(double, double), double (*scheme)(double), double* A, double* b) {
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
    printf("Computing internal nodes discretization coefficients...\n\n");

    // Initialize matrix of discretization coefficients (A) and vector of independent terms (b) to zero
    std::fill_n(A, 5*m.getNX()*m.getNY(), 0);
    std::fill_n(b, m.getNX()*m.getNY(), 0);

    // Internal nodes
    for(int j = 1; j < m.getNY()-1; j++) {
        for(int i = 1; i < m.getNX()-1; i++) {
            int node = j * m.getNX() + i;
            double x = m.atNodeX(i);
            double y = m.atNodeY(j);
            // South node
            double D = gamma * m.atSurfY(i) / m.atDistY(j-1);
            double F = rho * (*vy)(x, m.atFaceY(j)) * m.atSurfY(i);
            double P = F / D;
            A[5*node] = D * (*scheme)(P) + std::max(F, 0.0);
            // West node
            D = gamma * m.atSurfX(j) / m.atDistX(i-1);
            F = rho * (*vx)(m.atFaceX(i), y) * m.atSurfX(j);
            P = F / D;
            A[5*node+1] = D * (*scheme)(P) + std::max(F, 0.0);
            // East node
            D = gamma * m.atSurfX(j) / m.atDistX(i);
            F = rho * (*vx)(m.atFaceX(i+1), y) * m.atSurfX(j);
            P = F / D;
            A[5*node+2] = D * (*scheme)(P) + std::max(-F, 0.0);
            // North node
            D = gamma * m.atSurfY(i) / m.atDistY(j);
            F = rho * (*vy)(x, m.atFaceY(j+1)) * m.atSurfY(i);
            P = F / D;
            A[5*node+3] = D * (*scheme)(P) + std::max(-F, 0.0);
            // Central node
            A[5*node+4] = A[5*node] + A[5*node+1] + A[5*node+2] + A[5*node+3] - (*sourceP)(x,y) * m.atVol(i,j);
            // Independent term
            b[node] = (*sourceC)(x,y) * m.atVol(i,j);
        }
    }
}

void computeSteadyStateDiscretizationCoefficientsInternalNodes(const Mesh m, const double rho, const double gamma,
double (*vx)(double,double), double (*vy)(double, double), double (*sourceP)(double, double),
double (*sourceC)(double, double), int schemeCode, double* A, double* b) {



    double (*scheme) (const double);
    switch(schemeCode) {
        // UDS scheme
        case SCHEME_UDS:
            scheme = &schemeUDS;
            break;
        // CDS scheme
        case SCHEME_CDS:
            scheme = &schemeCDS;
            break;
        // EDS scheme
        case SCHEME_EDS:
            scheme = &schemeEDS;
            break;
        // Hybrid scheme
        case SCHEME_HYBRID:
            scheme = &schemeHybrid;
            break;
        // Powerlaw scheme
        default:    // Powerlaw scheme
            scheme = &schemePowerlaw;
    }

    // Initialize matrix of discretization coefficients (A) and vector of independent terms (b) to zero
    std::fill_n(A, 5*m.getNX()*m.getNY(), 0);
    std::fill_n(b, m.getNX()*m.getNY(), 0);

    printf("Computing internal nodes discretization coefficients...\n\n");
    // Internal nodes
    for(int j = 1; j < m.getNY()-1; j++) {
        for(int i = 1; i < m.getNX()-1; i++) {
            int node = j * m.getNX() + i;
            double x = m.atNodeX(i);
            double y = m.atNodeY(j);
            // South node
            double D = gamma * m.atSurfY(i) / m.atDistY(j-1);
            double F = rho * (*vy)(x, m.atFaceY(j)) * m.atSurfY(i);
            double P = F / D;
            A[5*node] = D * (*scheme)(P) + std::max(F, 0.0);
            // West node
            D = gamma * m.atSurfX(j) / m.atDistX(i-1);
            F = rho * (*vx)(m.atFaceX(i), y) * m.atSurfX(j);
            P = F / D;
            A[5*node+1] = D * (*scheme)(P) + std::max(F, 0.0);
            // East node
            D = gamma * m.atSurfX(j) / m.atDistX(i);
            F = rho * (*vx)(m.atFaceX(i+1), y) * m.atSurfX(j);
            P = F / D;
            A[5*node+2] = D * (*scheme)(P) + std::max(-F, 0.0);
            // North node
            D = gamma * m.atSurfY(i) / m.atDistY(j);
            F = rho * (*vy)(x, m.atFaceY(j+1)) * m.atSurfY(i);
            P = F / D;
            A[5*node+3] = D * (*scheme)(P) + std::max(-F, 0.0);
            // Central node
            A[5*node+4] = A[5*node] + A[5*node+1] + A[5*node+2] + A[5*node+3] - (*sourceP)(x,y) * m.atVol(i,j);
            // Independent term
            b[node] = (*sourceC)(x,y) * m.atVol(i,j);
        }
    }

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DIAGONAL CASE FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int solveDiagonalCase(const double L, const double lz, const int N, const double rho, const double gamma, const double phi_low, const double phi_high) {
    // Build mesh
    printf("Building uniform cartesian mesh...\n\n");
    Mesh m;
    m.buildUniformMesh(0, 0, L, L, lz, N, N);
    if(!m.isBuilt()) {
        printf("\tError. Could not build mesh\n");
        return -2;
    }

    // Allocate memory
    double* A = (double*) malloc(5 * m.getNX() * m.getNY() * sizeof(double*));
    if(!A) {
        printf("Error. Could not allocate memory for the linear system matrix.\n");
        return -2;
    }

    double* b = (double*) malloc(m.getNX() * m.getNY() * sizeof(double*));
    if(!b) {
        printf("Error. Could not allocate memory for the linear system vector.\n");
        return -2;
    }

    // Compute discretization coefficients
    computeSteadyStateDiscretizationCoefficientsInternalNodes(m, rho, gamma, vxDiagonal, vyDiagonal, sourcePDiagonal, sourceCDiagonal, SCHEME_UDS, A, b);
    computeDiscCoefsBoundaryNodesDiagonal(m, phi_low, phi_high, A, b);

    // Allocate memory for the linear system solution
    double* phi = (double*) malloc(m.getNX() * m.getNY() * sizeof(double*));
    if(!phi) {
        printf("Error. Could not allocate memory for the linear system solution vector.\n");
        return -2;
    }
    std::fill_n(phi, m.getNX()*m.getNY(), 1);
    solveSystem(m.getNX(), m.getNY(), A, b, phi, 0);

    // const char* filename = "output/diagonal.dat";

    char filename[100];
    sprintf(filename, "output/diagonal_N%d_Pe%.1e.dat", m.getNX(), rho/gamma);

    printToFile(m, phi, filename, 5);
    plotSolution(m, filename);

    // Free memory allocated
    free(A);
    free(b);
    free(phi);
    return 1;
}

void computeDiscCoefsBoundaryNodesDiagonal(const Mesh m, const double phi_low, const double phi_high, double* A, double* b) {
    /*
    computeDiscCoefsBoundaryNodesDiagonal: computes the discretization coefficients for the boundary nodes in the diagonal case
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - m                 Mesh object                                                                     [const Mesh]
        - phi_low           Minimum value of phi on the boundary                                            [const double]
        - phi_high          Maximum value of phi on the boundary                                            [const double]
        - A                 Linear system matrix set to zero. Rows: nx*ny, Columns: 5                       [double*]
        - b                 Vector of indenpendent terms set to zero. Rows: nx*ny, Columns: 1               [double*]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - A                 Linear system matrix. Rows: nx*ny, Columns: 5                                   [double*]
        - b                 Vector of indenpendent terms. Rows: nx*ny, Columns: 1                           [double*]
    */
    printf("Computing boundary nodes discretization coefficients for the diagonal case...\n\n");
    // Lower row nodes: 0 <= i <= nx-1; j=0
    for(int i = 0; i < m.getNX(); i++) {
        A[5*i+4] = 1;
        b[i] = phi_low;
    }
    // Right column nodes: i = nx-1; 1 <= j <= ny-1
    for(int j = 1; j < m.getNY(); j++) {
        int node = (j + 1) * m.getNX() - 1;
        A[5*node+4] = 1;
        b[node] = phi_low;
    }
    // Left column nodes: i = 0; 1 <= j <= ny-1
    for(int j = 1; j < m.getNY(); j++) {
        int node = j * m.getNX();
        A[5*node+4] = 1;
        b[node] = phi_high;
    }
    // Upper row nodes: 1 <= i <= nx-2, j = ny-1
    for(int i = 1; i < m.getNX()-1; i++) {
        int node = (m.getNY() - 1) * m.getNX() + i;
        A[5*node+4] = 1;
        b[node] = phi_high;
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
    return 1*cos(0.25*M_PI);
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
    return 1*sin(0.25*M_PI);
}

double sourcePDiagonal(const double x, const double y) {
    /*
    sourcePDiagonal: computes the SP coefficient of the source term for the diagonal case at the point (x,y)
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

double sourceCDiagonal(const double x, const double y) {
    /*
    sourcePDiagonal: computes the SC coefficient of the source term for the diagonal case at the point (x,y)
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - x         x coordinate    [const double]
        - y         y coordinate    [const double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - source    SC coefficient of the source term for the diagonal case at the point (x,y)  [double]
    */
    return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SMITH-HUTTON CASE FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int solveSmithHuttonCase(const double L, const double lz, const int N, const double rho, const double gamma) {

    // Build mesh
    const int nx = (N % 2 == 1 ? N : N+1);
    const int ny = 0.5*(nx + 1);

    printf("Building uniform cartesian mesh...\n\n");
    Mesh m;
    m.buildUniformMesh(-L, 0, 2*L, L, lz, nx, ny);
    if(!m.isBuilt()) {
        printf("\tError. Could not build mesh\n");
        return -2;
    }

    // Allocate memory
    double* A = (double*) malloc(5 * m.getNX() * m.getNY() * sizeof(double*));
    if(!A) {
        printf("Error. Could not allocate memory for the linear system matrix.\n");
        return -2;
    }

    double* b = (double*) malloc(m.getNX() * m.getNY() * sizeof(double*));
    if(!b) {
        printf("Error. Could not allocate memory for the linear system vector.\n");
        return -2;
    }

    // Compute discretization coefficients
    computeSteadyStateDiscretizationCoefficientsInternalNodes(m, rho, gamma, vxSmithHutton, vySmithHutton, sourcePSmithHutton, sourceCSmithHutton, SCHEME_POWERLAW, A, b);
    computeDiscretizationCoefficientsBoundaryNodesSmithHuttonCase(m, A, b);

    // Allocate memory for the linear system solution
    double* phi = (double*) malloc(m.getNX() * m.getNY() * sizeof(double*));
    if(!phi) {
        printf("Error. Could not allocate memory for the linear system solution vector.\n");
        return -2;
    }
    std::fill_n(phi, m.getNX()*m.getNY(), 1);
    solveSystem(m.getNX(), m.getNY(), A, b, phi, 0);

    // const char* filename = "output/smith_hutton.dat";
    // printToFile(m, phi, filename, 5);
    // plotSolution(m, filename);

    char filename[250];
    sprintf(filename, "../gnuplot/input/case_smith_hutton/smith_hutton_N%d_rg%.1e.dat", m.getNX(), rho/gamma);

    printToFile(m, phi, filename, 5);
    // plotSolution(m, filename);

    // Free memory allocated
    free(A);
    free(b);
    free(phi);
    return 1;
}

void computeDiscretizationCoefficientsBoundaryNodesSmithHuttonCase(const Mesh m, double* A, double* b) {
    /*
    */
    printf("Computing boundary nodes discretization coefficients for Smith-Hutton case...\n\n");
    // Lower row nodes: 0 <= i <= nx-1; j=0
    // Lower boundary: [-L, 0]
    for(int i = 0; i < m.getNX() && m.atNodeX(i) <= 0; i++) {
        A[5*i+4] = 1;
        b[i] = 1 + tanh(10*(2*m.atNodeX(i)+1));
    }

    // Lower boundary: (0,L]
    for(int i = m.getNX()-1; i >= 0 && m.atNodeX(i) > 0; --i) {
        A[5*i+3] = 1;
        A[5*i+4] = 1;
    }

    // Left and right boudnaries
    for(int j = 1; j < m.getNY(); j++) {
        // Left boundary
        int node = j*m.getNX();
        A[5*node+4] = 1;
        b[node] = 1 - tanh(10);
        // Right boundary
        node = (j+1)*m.getNX() - 1;
        A[5*node+4] = 1;
        b[node] = 1 - tanh(10);
    }

    // Top boundary
    for(int i = 1; i < m.getNX()-1; i++) {
        int node = (m.getNY()-1)*m.getNX() + i;
        A[5*node+4] = 1;
        b[node] = 1 - tanh(10);
    }
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

double sourcePSmithHutton(const double, const double) {
    /*
    sourcePSmithHutton: computes the SP coefficient of the source term for the Smith-Hutton case at the point (x,y)
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

double sourceCSmithHutton(const double, const double) {
    /*
    sourceCSmithHutton: computes the SC coefficient of the source term for the Smith-Hutton case at the point (x,y)
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - x         x coordinate    [const double]
        - y         y coordinate    [const double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - source    SC coefficient of the source term for the Smith-Hutton case at the point (x,y)  [double]
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
    printf("Printing the solution to file '%s'...\n", filename);
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

void plotSolution(const Mesh m, const char* filename) {
    /*
    plotSolution: plots phi vs (x,y) in a temperature map using gnuplot.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - m                 Mesh of the problem                                                 [const Mesh]
        - filename          Name of the file where the solution has been previously written     [const char*]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: temperature map of phi. Prior to running the program,
        - the command "export DISPLAY=:0 gnuplot" must be executed (only needed before the first execution)
        - Xming server must be running
    */
    // Set xrange command
    char xrange[50];
    sprintf(xrange, "set xrange [%.5f : %.5f]", m.atNodeX(0), m.atNodeX(m.getNX()-1));

    // Set yrange command
    char yrange[50];
    sprintf(yrange, "set yrange [%.5f : %.5f]", m.atNodeY(0), m.atNodeY(m.getNY()-1));

    // Set size ratio command
    char sizeratio[50];
    sprintf(sizeratio, "set size ratio %.5f", (m.atNodeY(m.getNY()-1)-m.atNodeY(0))/(m.atNodeX(m.getNX()-1) - m.atNodeX(0)));

    // Plotting command, uses filename parameter
    char plotCommand[30+strlen(filename)];
    sprintf(plotCommand, "plot '%s' with image", filename);

    // Commands sent to gnuplot
    const char* GnuCommands[] = {xrange, yrange, sizeratio, "set palette rgb 33,13,10", plotCommand};

    // Send commands to gnuplot
    FILE *gnupipe = NULL;
    gnupipe = popen("gnuplot -persistent", "w");
    printf("Plotting the solution in gnuplot...\n");
    for(int i = 0; i < 5; i++) {
        printf("\t%s\n", GnuCommands[i]);
        fprintf(gnupipe, "%s\n", GnuCommands[i]);
    }
    printf("\n");
}


void printVelocityField(const Mesh m, double (*vx)(double, double), double (*vy)(double, double), const char* filename_mod, const char* filename_vec, const int precision) {
    /*
    printVelocityField: prints the velocity field to a file. filename_mod contains the norm of the velocity field at each point. filename_vec contains
    the norm and the components of the velocity field at each point.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - m                 Mesh of the problem                                                                         [const Mesh]
        - filename          Name of the file where the solution has been previously written                             [const char*]
        - *vx               Function that gives the velocity in the x axis provided the (x,y) coordinates               [returns double, needs (double,double)]
        - *vy               Function that gives the velocity in the y axis provided the (x,y) coordinates               [returns double, needs (double,double)]
        - filename_mod      File where the norm of the velocity field at each point is printed.                         [const char*]
        - filename_vecd     File where the norm and the components of the velocity field at each point are printed.     [const char*]
        - precision         Number of decimal places to which doubles are printed.                                      [const int]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - filename_mod      File where the norm of the velocity field at each point is printed.                         [const char*]
        - filename_vecd     File where the norm and the components of the velocity field at each point are printed.     [const char*]
    */
    // Print the (x,y) coordinates and the norm of the velocity field at (x,y)
    printf("Printing the velocity field norm...\n");
    std::ofstream file;
    file.open(filename_mod);
    if(file.is_open()) {
        printf("\tWriting to file '%s'...\n", filename_mod);
        file << std::setprecision(precision) << std::fixed;
        for(int j = 0; j < m.getNY(); j++) {
            for(int i = 0; i < m.getNX(); i++) {
                double x = m.atNodeX(i);
                double y = m.atNodeY(j);
                double vel = std::sqrt(vx(x,y)*vx(x,y) + vy(x,y)*vy(x,y));
                file << x << " " << y << " " << vel << std::endl;
            }
            file << std::endl;
        }
        file << std::endl;
    } else
        printf("\tCould not open file '%s'...\n", filename_mod);
    printf("\tClosing file...\n\n");
    file.close();

    // Print the (x,y) coordinates, the norm and the components of the velocity field at (x,y)
    printf("Printing the velocity field norm and components...\n");
    file.open(filename_vec);
    if(file.is_open()) {
        printf("\tWriting to file '%s'...\n", filename_vec);
        file << std::setprecision(precision) << std::fixed;
        for(int j = 0; j < m.getNY(); j++) {
            for(int i = 0; i < m.getNX(); i++) {
                double x = m.atNodeX(i);
                double y = m.atNodeY(j);
                double velx = vx(x,y);
                double vely = vy(x,y);
                double vel = std::sqrt(velx*velx + vely*vely);
                file << x << " " << y << " " << vel << " " << velx << " " << vely << std::endl;
            }
            file << std::endl;
        }
    } else
        printf("\tCould not open file '%s'...\n", filename_vec);
    printf("\tClosing file...\n\n");
    file.close();
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CHECK SOLUTION FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double checkLinearSystemSolution(const int nx, const int ny, const double* A, const double* b, const double* phi) {
    /*
    checkLinearSystemSolution: checks the solution of the linear system resulting from the 2D convection-diffusion equations in a cartesian mesh. If _A is
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

struct ProblemInput {
    int type;
    double L;
    double lx;
    double ly;
    double lz;
    double x0;
    double y0;
    double rho;
    double gamma;
    double phi_low;
    double phi_high;
    int N;
    const char* filename;
};

int parseInputFile(const char* filename, ProblemInput &p) {

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

    if(caseDiagonal.compare(word) == 0) {
        p.type = 1;
        double L;
        iss >> L;
        p.lx = L;
        p.ly = L;
        p.filename = "output/outputDiagonal.dat";
    } else if(caseSmithHutton.compare(word) == 0) {
        p.type = 2;
        double L;
        iss >> L;
        p.lx = L;
        p.ly = L;
        p.filename = "output/outputSmithHutton.dat";
    } else {
        printf("Error. Case provided (%s) is not valid.\n", word.c_str());
        return -1;
    }

    // Read lz
    file >> p.lz;
    // Read x0
    file >> p.x0;
    // Read y0
    file >> p.y0;

    // Read density
    file >> p.rho;
    if(p.rho <= 0) {
        printf("Error. The density provided (%.5e) must be positive.\n", p.rho);
        return -1;
    }

    // Read diffusion coefficient
    file >> p.gamma;
    if(p.gamma <= 0) {
        printf("Error. The diffusion coefficient provided (%.5e) must be positive.\n", p.gamma);
        return -1;
    }

    // Read phi_low and phi_high if it is the diagonal case
    if(p.type == 1) {
        file >> p.phi_low;
        file >> p.phi_high;
    }

    // Read N
    file >> p.N;
    if(p.N < 2) {
        printf("Error. Number of control volumes provided (%d) must be at least 2.\n", p.N);
        return -1;
    }

    if(p.type == 2 && p.N % 2 == 0) {
        printf("Warning. Number of control volumes provided (%d) for the Smith-Hutton case has been increased in one.\n", p.N);
        p.N++;
    }

    file.close();
    return 0;
}
