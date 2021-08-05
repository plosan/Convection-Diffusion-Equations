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

double vxDiagonal(double x, double y) {
    // return V0*cos(ALPHA);
    return cos(0.25*M_PI);
}

double vyDiagonal(double x, double y) {
    // return V0*sin(ALPHA);
    return sin(0.25*M_PI);
}

double sourceDiagonal(double x, double y) {
    return 0;
}

void computeSteadyStateDiscretizationCoefficientsInternalNodes(const Mesh m, const double* prop,
double (*vx)(double,double), double (*vy)(double, double), double (*source)(double, double), double* A, double* b, const int scheme);

void computeDiscretizationCoefficientsBoundaryNodesDiagonalCase(const Mesh m, const double* phi_boundary, double* A, double* b);

void computeDiscretizationCoefficientsInternalNodes(const Mesh m, const double* properties,
double (*vx)(double,double), double (*vy)(double, double), double* A, double* b, const int scheme);

void checkSystemMatrix(const unsigned int nx, const unsigned int ny, const double tol, const double* A);

void solveSystem(const unsigned int nx, const unsigned int ny, const double tol, const unsigned int maxIt, const double* A, const double* b, double* phi);

void printToFile(const Mesh m, const double* phi, const char* filename, const int precision);

void plotSolution(const char* filename);

void assembleMatrix(const unsigned int nx, const unsigned int ny, const double* A, double** A_mat);


double checkSystemSolution(const unsigned int nx, const unsigned int ny, const double* A, const double* b, const double* phi);

void verification(const double lx, const double ly, const double lz, const unsigned int nx, const unsigned int ny, const double* nodeX, const double* nodeY,
const double* distX, const double* distY, const double* faceX, const double* faceY, const double* surfX, const double* surfY, const double* vol,
const double* phi_boundary, const double* v, const double rho, const double gamma, double* A, double* b, double* phi, const int scheme);



int main(int arg, char* argv[]) {

    // // Physical data
    // Sizes
    double x0 = 0;  // Lower left corner x coordinate for rectangular domain    [m]
    double y0 = 0;  // Lower left corner y coordinate for rectangular domain    [m]
    double L = 1;   // Domain size in x and y axis                              [m]
    double lz = 1;  // Domain size in z axis                                    [m]

    // Boundary conditions
    const double phi_low = 10;      // Minimum value for phi
    const double phi_high = 20;     // Maximum value for phi

    // Thermophysical properties for water at 20 ºC
    const double lambda = 0.5861;       // Thermal conductivity                         [W/(k*m)]
    const double cv = 4183;             // Specific heat at constant volume (pressure)  [J/(kg*K)]
    const double rho = 998.2;           // Density                                      [kg/m^3]
    const double gamma = lambda / cv;   // Diffusion coefficient
    double* prop = (double*) malloc(2 * sizeof(double*));
    prop[0] = rho;
    prop[1] = gamma;

    // // Numerical data
    unsigned int N = 75;
    unsigned int nx = N;      // Number of nodes in x axis
    unsigned int ny = N;      // Number of nodes in y axis
    const double phi0 = 1;      // Initial value to fill phi vector for linear system resolution
    const double tol = 1e-15;   // Tolerance to stop iteration
    const int maxIt = 1e6;      // Maximum number of iterations

    double* phi_boundary = (double*) malloc(2 * sizeof(double*));
    phi_boundary[0] = phi_low;
    phi_boundary[1] = phi_high;


    Mesh m;
    printf("Building uniform cartesian mesh...\n\n");
    m.buildUniformMesh(x0, y0, L, L, lz, nx, ny);
    if(!m.isBuilt())
        return -2;

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
    computeDiscretizationCoefficientsBoundaryNodesDiagonalCase(m, phi_boundary, A, b);

    double* phi = (double*) malloc(nx * ny * sizeof(double*));
    std::fill_n(phi, nx*ny, phi0);

    printf("Computing boundary nodes discretization coefficients for the diagonal case...\n\n");
    solveSystem(nx, ny, tol, maxIt, A, b, phi);

    const char* filename = "output/output.dat";
    printToFile(m, phi, filename, 5);
    plotSolution(filename);

    double checkSol = checkSystemSolution(nx, ny, A, b, phi);
    printf("checkSol : %.5e\n\n\n", checkSol);

    // Free memory allocated
    free(prop);
    free(phi_boundary);
    free(A);
    free(b);
    free(phi);

    return 1;
}

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
    for(unsigned int j = 1; j < m.getNY()-1; j++) {
        for(unsigned int i = 1; i < m.getNX()-1; i++) {
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
}

void computeDiscretizationCoefficientsBoundaryNodesDiagonalCase(const Mesh m, const double* phi_boundary, double* A, double* b) {
    /*
    computeDiscretizationCoefficientsBoundaryNodesDiagonalCase: computes the discretization coefficients for the boundary nodes in the diagonal case
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
    printf("Computing boundary nodes discretization coefficients for the diagonal case...\n\n");
    // Lower row nodes: 0 <= i <= nx-1; j=0
    for(unsigned int i = 0; i < m.getNX(); i++) {
        A[5*i+4] = 1;
        b[i] = phi_boundary[0];
    }
    // Right column nodes: i = nx-1; 1 <= j <= ny-1
    for(unsigned int j = 1; j < m.getNY(); j++) {
        int node = (j + 1) * m.getNX() - 1;
        A[5*node+4] = 1;
        b[node] = phi_boundary[0];
    }
    // Left column nodes: i = 0; 1 <= j <= ny-1
    for(unsigned int j = 1; j < m.getNY(); j++) {
        int node = j * m.getNX();
        A[5*node+4] = 1;
        b[node] = phi_boundary[1];
    }
    // Upper row nodes: 1 <= i <= nx-2, j = ny-1
    for(unsigned int i = 1; i < m.getNX()-1; i++) {
        int node = (m.getNY() - 1) * m.getNX() + i;
        A[5*node+4] = 1;
        b[node] = phi_boundary[1];
    }
}


void checkSystemMatrix(const unsigned int nx, const unsigned int ny, const double tol, const double* A) {
    /*
    checkSystemMatrix: checks
        - whether the system matrix has a row full of zeros (sufficient condition for incompatible system but not necessary)
        - whether some coefficient aP (5th column) has a zero (prevent division by 0)
    In either case the function prints a message.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - nx: discretization nodes in X axis    [const unsigned int]
        - ny: discretization nodes in Y axis    [const unsigned int]
        - tol: tolerance                        [const double]
        - A: linear system matrix               [double*]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: none
    */
    printf("Checking linear system matrix...\n");
    // Check whether the system matrix has a row full of zeros (sufficient condition for incompatible system but not necessary)
    bool foundNullRow = false;  // Bool variable: true => a null row has been found; false => no null row has been found
    unsigned int node = 0;      // Node (row) being checked
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

void solveSystem(const unsigned int nx, const unsigned int ny, const double tol, const unsigned int maxIt, const double* A, const double* b, double* phi) {
    /*
    solveSystem: solves the linear system resulting from a 2D convection-diffusion problem in a domain discretized with a cartesian mesh using
    Gauss-Seidel algorithm. It has two criterion to stop the iteration:
        - Let phi* and phi be two consecutive vectors of the sequence produced by Gauss-Seidel algorithm. If the infinity norm of phi-phi* is less
        than tol, then the algorithm stops.
        - If the number of iterations reach maxIt, the algorithm stops.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - nx        Number of nodes in the X axis                               [const double]
        - ny        Number of nodes in the Y axis                               [const double]
        - tol       Tolerance to stop iteration                                 [const double]
        - maxIt     Maximum number of iterations                                [const unsigned int]
        - A         Linear system matrix. Rows: nx*ny, Columns: 5               [double*]
        - b         Vector of indenpendent terms. Rows: nx*ny, Columns: 1       [double*]
        - phi       Solution of the linar system. Rows: nx*ny, Columns: 1       [double*]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - A                 Linear system matrix. Rows: nx*ny, Columns: 5                                   [double*]
        - b                 Vector of indenpendent terms. Rows: nx*ny, Columns: 1                           [double*]
    */
    printf("Solving linear system...\n");
    unsigned int it = 0;        // Current iteration
    bool convergence = false;   // Boolean variable to tell whether there is convergence or not. False: no convergence, True: convergence
    // Gauss-Seidel iteration
    while(it < maxIt && !convergence) {
        double maxDiff = -1;    // Infinity norm of the difference phi-phi*
        // Lower row nodes
        for(unsigned int i = 0; i < nx; i++) {
            int node = i;                                                       // Node whose phi is being computed
            double aux = phi[node];                                             // Previous value of phi[node]
            phi[node] = (b[node] + A[5*node+3] * phi[node+nx]) / A[5*node+4];   // Compute new value
            maxDiff = std::max(maxDiff, std::abs(aux - phi[node]));             // Update infinity norm
        }
        // Central rows nodes
        for(unsigned int j = 1; j < ny-1; j++) {
            for(unsigned int i = 0; i < nx; i++) {
                int node = j * nx + i;                                          // Node whose phi is being computed
                double aux = phi[node];                                         // Previous value of phi[node]
                phi[node] = (b[node] + A[5*node] * phi[node-nx] + A[5*node+1] * phi[node-1] + A[5*node+2] * phi[node+1] + A[5*node+3] * phi[node+nx]) / A[5*node+4];    // Compute new value
                maxDiff = std::max(maxDiff, std::abs(aux - phi[node]));         // Update infinity norm
            }
        }
        // Upper row nodes
        for(unsigned int i = 0; i < nx; i++) {
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
        for(unsigned int i = 0; i < m.getNX(); i++) {
            for(unsigned int j = 0; j < m.getNY(); j++)
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

void assembleMatrix(const unsigned int nx, const unsigned int ny, const double* A, double** A_mat) {
    printf("Assembling matrix to check linear system solution...\n\n");
    if(A_mat) {
        // // Lower row
        for(unsigned int i = 0; i < nx; i++) {
            int node = i;
            A_mat[node][node+nx] = A[5*node+3]; // aN coefficient
            A_mat[node][node] = A[5*node+4];    // aP coefficient
        }
        // // Central rows
        for(unsigned int j = 1; j < ny-1; j++) {
            // First node
            int node = j * nx;
            A_mat[node][node+1] = A[5*node+2];  // aE coefficient
            A_mat[node][node] = A[5*node+4];    // aP coefficient
            // Central nodes
            for(unsigned int i = 1; i < nx-1; i++) {
                node = j * nx + i;
                A_mat[node][node-nx] = A[5*node];   // aS coefficient
                A_mat[node][node-1]  = A[5*node+1]; // aW coefficient
                A_mat[node][node+1]  = A[5*node+2]; // aE coefficient
                A_mat[node][node+nx] = A[5*node+3]; // aN coefficient
                A_mat[node][node]    = A[5*node+4]; // aP coefficient
            }
            // Last node
            node = j * nx + nx - 1;
            A_mat[node][node-1] = A[5*node+1];  // aW coefficient
            A_mat[node][node]   = A[5*node+4];  // aP coefficient
        }
        // // Upper row
        for(unsigned int i = 0; i < nx; i++) {
            int node = (ny - 1) * nx + i;
            A_mat[node][node-nx] = A[5*node];   // aS coefficient
            A_mat[node][node]    = A[5*node+4]; // aP coefficient
        }
    } else {
        printf("A_mat not allocated\n");
    }
}

double checkSystemSolution(const unsigned int nx, const unsigned int ny, const double* A, const double* b, const double* phi) {
    /*
    checkSystemSolution: checks the solution of the linear system resulting from the 2D convection-diffusion equations in a cartesian mesh. If _A is
    the linear system matrix, the function returns the infinity norm of _A * phi - b.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - nx        Number of nodes in the X axis                               [const double]
        - ny        Number of nodes in the Y axis                               [const double]
        - A         Linear system matrix. Rows: nx*ny, Columns: 5               [double*]
        - b         Vector of indenpendent terms. Rows: nx*ny, Columns: 1       [double*]
        - phi       Solution of the linar system. Rows: nx*ny, Columns: 1       [double*]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - maxDiff   Infinity norm of the difference A * phi - b                 [double]
    */
    double maxDiff = -1;
    // Lower row nodes
    for(unsigned int i = 0; i < nx; i++) {
        double diff = b[i] - A[5*i+3] * phi[i+nx] - A[5*i+4] * phi[i];
        maxDiff = std::max(maxDiff, std::abs(diff));
    }
    // Central row nodes
    for(unsigned int j = 1; j < ny-1; j++) {
        for(unsigned int i = 0; i < nx; i++) {
            int node = j * nx + i;
            double diff = b[node] - (A[5*node]*phi[node-nx]) - (A[5*node+1]*phi[node-1]) - (A[5*node+2]*phi[node+1]) - (A[5*node+3]*phi[node+nx]) - (A[5*node+4]*phi[node]);
            maxDiff = std::max(maxDiff, std::abs(diff));
        }
    }
    // Upper row nodes
    for(unsigned int i = 0; i < nx; i++) {
        int node = (ny-1) * nx + i;
        double diff = b[node] - (A[5*node] * phi[node-nx]) - (A[5*node+4] * phi[node]);
        maxDiff = std::max(maxDiff, std::abs(diff));
    }
    return maxDiff;
}

void verification(const double lx, const double ly, const double lz, const unsigned int nx, const unsigned int ny, const double* nodeX, const double* nodeY,
const double* distX, const double* distY, const double* faceX, const double* faceY, const double* surfX, const double* surfY, const double* vol,
const double* phi_boundary, const double* v, const double rho, const double gamma, double* A, double* b, double* phi, const int scheme) {
    printf("Verificating solution...\n");
    double stepX = lx / (nx - 1);
    double stepY = ly / (ny - 1);
    double maxDiff = -1;
    for(unsigned int i = 1; i < nx-1; i++) {
        for(unsigned int j = 1; j < ny-1; j++) {
            int node = j * nx + i;
            double phi_xx = (phi[node+1] - 2*phi[node] + phi[node-1]) / (stepX * stepX);
            double phi_yy = (phi[node+nx] - 2*phi[node] + phi[node-nx]) / (stepY * stepY);
            double RHS = gamma * (phi_xx + phi_yy);

            double phi_x = (phi[node+1] - phi[node-1]) / (2*stepX);
            double phi_y = (phi[node+nx] - phi[node-nx]) / (2*stepY);
            double LHS = rho * (v[0] * phi_x + v[1] * phi_y);

            maxDiff = std::max(maxDiff, std::abs(RHS - LHS));
        }
    }
    printf("\tMaximum difference: %.5e\n\n", maxDiff);
}

void computeDiscretizationCoefficientsDiagonalCase(const unsigned int nx, const unsigned int ny, const double* nodeX, const double* nodeY,
const double* distX, const double* distY, const double* faceX, const double* faceY, const double* surfX, const double* surfY, const double* vol,
const double* phi_boundary, const double* v, const double rho, const double gamma, double* A, double* b, const int scheme) {

    printf("Computing discretization coefficients for the diagonal case...\n\n");
    // Initialize matrix of discretization coefficients (A) and vector of independent terms (b) to zero
    std::fill_n(A, 5*nx*ny, 0);
    std::fill_n(b, nx*ny, 0);

    // Internal nodes, 1 <= i <= nx-2; 1 <= j <= ny-2
    for(unsigned int j = 1; j < ny-1; j++) {
        for(unsigned int i = 1; i < nx-1; i++) {
            int node = j * nx + i;
            printf("%10s : %d\n", "Node", node);
            // South node
            double massFlow = -rho * v[1] * surfY[i];
            double c = (massFlow != 0 ? (massFlow + std::abs(massFlow))/(2 * massFlow) : 0);
            A[5*node] = gamma*surfY[i]/distY[j-1] + massFlow*c;
            // West node
            massFlow = -rho * v[0] * surfX[j];
            c = (massFlow != 0 ? (massFlow + std::abs(massFlow))/(2 * massFlow) : 0);
            A[5*node+1] = gamma*surfX[j]/distX[i-1] + massFlow*c;
            // East node
            massFlow = rho * v[0] * surfX[j];
            c = (massFlow != 0 ? (massFlow - std::abs(massFlow))/(2 * massFlow) : 0);
            A[5*node+2] = gamma*surfX[j]/distX[i] + massFlow*c;
            // North node
            massFlow = rho * v[1] * surfY[i];
            c = (massFlow != 0 ? (massFlow - std::abs(massFlow))/(2 * massFlow) : 0);
            A[5*node+3] = gamma*surfY[i]/distY[j] + massFlow*c;
            // Central node
            A[5*node+4] = A[5*node] + A[5*node+1] + A[5*node+2] + A[5*node+3];
        }
    }

    // Lower row nodes: 0 <= i <= nx-1; j=0
    for(unsigned int i = 0; i < nx; i++) {
        A[5*i+4] = 1;
        b[i] = phi_boundary[0];
    }

    // Right column nodes: i = nx-1; 1 <= j <= ny-1
    for(unsigned int j = 1; j < ny; j++) {
        int node = (j + 1) * nx - 1;
        A[5*node+4] = 1;
        b[node] = phi_boundary[0];
    }

    // Left column nodes: i = 0; 1 <= j <= ny-1
    for(unsigned int j = 1; j < ny; j++) {
        int node = j * nx;
        A[5*node+4] = 1;
        b[node] = phi_boundary[1];
    }

    // Upper row nodes: 1 <= i <= nx-2, j = ny-1
    for(unsigned int i = 1; i < nx-1; i++) {
        int node = (ny - 1) * nx + i;
        A[5*node+4] = 1;
        b[node] = phi_boundary[1];
    }
}
