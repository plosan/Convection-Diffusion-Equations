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
#include "mesh.h"

void printToFile(const unsigned int nx, const unsigned int ny, const double* nodeX, const double* nodeY, const double* phi, std::string fileName, const int precision) {
    std::ofstream file;
    file.open(fileName);
    if(file.is_open()) {
        printf("Writing to file...\n");
        file << std::setprecision(precision) << std::fixed;
        for(unsigned int i = 0; i < nx; i++) {
            for(unsigned int j = 0; j < ny; j++)
                file << nodeX[i] << " " << nodeY[j] << " " << phi[j*nx+i] << std::endl;
            file << std::endl;
        }
    } else {
        printf("\tCould not open file\n");
    }
    file.close();
    printf("\n");
}

void checkSystemSolution(const unsigned int nx, const unsigned int ny, double** A_mat, const double* b, const double* phi);

void assembleMatrix(const unsigned int nx, const unsigned int ny, const double* A, double** A_mat);

void checkSystemMatrix(const unsigned int nx, const unsigned int ny, const double tol, const double* A);

void computeDiscretizationCoefficientsDiagonalCase(const unsigned int nx, const unsigned int ny, const double* nodeX, const double* nodeY,
    const double* distX, const double* distY, const double* faceX, const double* faceY, const double* surfX, const double* surfY, const double* vol,
    const double* phi_boundary, const double* v, const double rho, const double gamma, double* A, double* b, const int scheme);

void solveSystem(const unsigned int nx, const unsigned int ny, const double tol, const unsigned int maxIt, const double* A, const double* b, double* phi) {
    printf("Solving linear system...\n");
    unsigned int it = 0;
    bool convergence = false;
    while(it < maxIt && !convergence) {
        double maxDiff = -1;
        // Lower row nodes
        for(unsigned int i = 0; i < nx; i++) {
            int node = i;
            double aux = phi[node];
            phi[node] = (b[node] + A[5*node+3] * phi[node+nx]) / A[5*node+4];
            maxDiff = std::max(maxDiff, std::abs(aux - phi[node]));
        }
        // Central rows nodes
        for(unsigned int j = 1; j < ny-1; j++) {
            for(unsigned int i = 0; i < nx; i++) {
                int node = j * nx + i;
                double aux = phi[node];
                phi[node] = (b[node] + A[5*node] * phi[node-nx] + A[5*node+1] * phi[node-1] + A[5*node+2] * phi[node+1] + A[5*node+3] * phi[node+nx]) / A[5*node+4];
                maxDiff = std::max(maxDiff, std::abs(aux - phi[node]));
            }
        }
        // Upper row nodes
        for(unsigned int i = 0; i < nx; i++) {
            int node = (ny - 1) * nx + i;
            double aux = phi[node];
            phi[node] = (b[node] + A[5*node] * phi[node-nx]) / A[5*node+4];
            maxDiff = std::max(maxDiff, std::abs(aux - phi[node]));
        }
        // Convergence condition
        convergence = (maxDiff < tol);
        // Increase iteration counter
        it++;
    }
    printf("\tIterations: %d\n\n", it);
}

void verification(const double lx, const double ly, const double lz, const unsigned int nx, const unsigned int ny, const double* nodeX, const double* nodeY,
const double* distX, const double* distY, const double* faceX, const double* faceY, const double* surfX, const double* surfY, const double* vol,
const double* phi_boundary, const double* v, const double rho, const double gamma, double* A, double* b, double* phi, const int scheme);

int main(int arg, char* argv[]) {

    double x0 = 0;
    double y0 = 0;
    unsigned int nx = 100;
    unsigned int ny = 100;
    double lx = 1;
    double ly = 1;
    double lz = 1;

    double* nodeX = (double*) malloc(nx * sizeof(double*));
    double* nodeY = (double*) malloc(ny * sizeof(double*));

    double* distX = (double*) malloc((nx - 1) * sizeof(double*));
    double* distY = (double*) malloc((ny - 1) * sizeof(double*));

    double* faceX = (double*) malloc((nx + 1) * sizeof(double*));
    double* faceY = (double*) malloc((ny + 1) * sizeof(double*));

    double* surfX = (double*) malloc(ny * sizeof(double*));
    double* surfY = (double*) malloc(nx * sizeof(double*));

    double* vol = (double*) malloc(nx * ny * sizeof(double*));

    compute2DUniformRectangularMesh(x0, y0, nx, ny, lx, ly, lz, nodeX, nodeY, distX, distY, faceX, faceY, surfX, surfY, vol);

    // printMeshInfo(x0, y0, nx, ny, lx, ly, lz, nodeX, nodeY, distX, distY, faceX, faceY, surfX, surfY, vol);


    double* A = (double*) malloc(5 * nx * ny * sizeof(double*));
    double* b = (double*) malloc(nx * ny * sizeof(double*));

    const double phi_low = 10;
    const double phi_high = 20;
    double* phi_boundary = (double*) malloc(2 * sizeof(double*));
    phi_boundary[0] = phi_low;
    phi_boundary[1] = phi_high;


    const double v0 = 1;
    const double alpha = 0.25 * M_PI;
    double* v = (double*) malloc(2 * sizeof(double*));
    v[0] = v0 * cos(alpha);
    v[1] = v0 * sin(alpha);

    // Thermophysical properties for water at 20 ºC
    const double lambda = 0.5861;       // Thermal conductivity                         [W/(k*m)]
    const double cv = 4183;             // Specific heat at constant volume (pressure)  [J/(kg*K)]
    const double rho = 998.2;           // Density                                      [kg/m^3]
    const double gamma = lambda / cv;   // Diffusion coefficient

    int scheme = 0;
    computeDiscretizationCoefficientsDiagonalCase(nx, ny, nodeX, nodeY, distX, distY, faceX, faceY, surfX, surfY, vol, phi_boundary, v, rho, gamma, A, b, scheme);

    const double tol = 1e-15;
    checkSystemMatrix(nx, ny, tol, A);

    const double phi0 = 1;
    double* phi = (double*) malloc(nx * ny * sizeof(double*));
    std::fill_n(phi, nx*ny, phi0);

    const int maxIt = 1e6;
    solveSystem(nx, ny, tol, maxIt, A, b, phi);

    std::string fileName("output/output.dat");
    printToFile(nx, ny, nodeX, nodeY, phi, fileName, 5);

    // printf("phi = \n");
    // printReversedRowMatrix(phi, ny, nx);

    double** A_mat = (double**) malloc(nx * ny * sizeof(double*));
    for(unsigned int row = 0; row < nx * ny; row++) {
        A_mat[row] = (double*) malloc(nx * ny * sizeof(double*));
    }

    assembleMatrix(nx, ny, A, A_mat);

    checkSystemSolution(nx, ny, A_mat, b, phi);

    for(unsigned int row = 0; row < nx * ny; row++) {
        free(A_mat[row]);
    }
    free(A_mat);

    printf("Verificating solution...\n");
    verification(lx, ly, lz, nx, ny, nodeX, nodeY, distX, distY, faceX, faceY, surfX, surfY, vol, phi_boundary, v, rho, gamma, A, b, phi, scheme);

    // Free memory allocated
    free(nodeX);
    free(nodeY);
    free(distX);
    free(distY);
    free(faceX);
    free(faceY);
    free(surfX);
    free(surfY);
    free(vol);

    free(A);
    free(b);
    free(phi_boundary);
    free(v);
    free(phi);

    return 0;
}

void checkSystemSolution(const unsigned int nx, const unsigned int ny, double** A_mat, const double* b, const double* phi) {
    printf("Checking system solution...\n");
    double maxDiff = -1;
    for(unsigned int row = 0; row < nx*ny; row++) {
        double sum = 0;
        for(unsigned int col = 0; col < nx*ny; col++)
            sum += A_mat[row][col] * phi[col];
        maxDiff = std::max(maxDiff, std::abs(sum - b[row]));
    }
    printf("\tMaximum difference: %.5e\n\n", maxDiff);
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

void verification(const double lx, const double ly, const double lz, const unsigned int nx, const unsigned int ny, const double* nodeX, const double* nodeY,
const double* distX, const double* distY, const double* faceX, const double* faceY, const double* surfX, const double* surfY, const double* vol,
const double* phi_boundary, const double* v, const double rho, const double gamma, double* A, double* b, double* phi, const int scheme) {

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

// void solveSystem(const unsigned int nx, const unsigned int ny, const double tol, const double* A, const double* b, double* phi) {
//     bool convergence = false;
//     while(!convergence) {
//         double maxDiff = -1;
//         // Lower row
//         double aux = phi[node];
//
//         int node = 0;
//         phi[node] = b[node]/A[5*node];
//         double diff = std::abs(aux - phi[node]);
//         maxDiff = (diff > maxDiff ? diff : maxDiff);
//
//         for(int i = 1; i < nx; i++) {
//             node = i;
//             phi[node] = (b[node] - A[5*node])
//         }
//
//         // Central rows
//
//         // Upper row
//     }
// }



void computeDiscretizationCoefficientsDiagonalCase(const unsigned int nx, const unsigned int ny, const double* nodeX, const double* nodeY, const double* distX, const double* distY, const double* faceX, const double* faceY, const double* surfX, const double* surfY, const double* vol,
const double* phi_boundary, const double* v, const double rho, const double gamma, double* A, double* b, const int scheme) {

    printf("Computing discretization coefficients for the diagonal case...\n");
    // Initialize matrix of discretization coefficients (A) and vector of independent terms (b) to zero
    std::fill_n(A, 5*nx*ny, 0);
    std::fill_n(b, nx*ny, 0);

    // Internal nodes, 1 <= i <= nx-2; 1 <= j <= ny-2
    for(unsigned int i = 1; i < nx-1; i++) {
        for(unsigned int j = 1; j < ny-1; j++) {
            int node = j * nx + i;
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

    printf("\n");
}
