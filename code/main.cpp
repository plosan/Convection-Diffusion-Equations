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

void smithHuttonCase();

void computeDiscretizationCoefficientsDiagonalCase(const unsigned int nx, const unsigned int ny, const double* nodeX, const double* nodeY,
    const double* distX, const double* distY, const double* faceX, const double* faceY, const double* surfX, const double* surfY, const double* vol,
    const double* phi_boundary, const double* v, const double rho, const double gamma, double* A, double* b, const int scheme);

void solveSystem(const unsigned int nx, const unsigned int ny, const double tol, const unsigned int maxIt, const double* A, const double* b, double* phi, unsigned int &it) {
    it = 0;
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

    // // // Check if there is a row full of zeros in the matrix
    // unsigned int node = 0;
    // bool null = false;
    // const double tol = 1e-12;
    // while(node < nx*ny && !null) {
    //     bool found = false;
    //     int row = 0;
    //     while(row < 5 && not found) {
    //         if(std::abs(A[5*node+row]) > tol)
    //             found = true;
    //         row++;
    //     }
    //     null = !found;
    //     node++;
    // }
    // if(null) {
    //     printf("found = %d\n", null);
    //     printf("count = %d\n", node-1);
    // }
    //
    // // // Check if there is a zero in the 5th column
    // node = 0;
    // null = false;
    // while(node < nx*ny && !null) {
    //     null = (std::abs(A[5*node+4]) < tol);
    //     node++;
    // }
    // if(null) {
    //     printf("Null: %d\n", null);
    //     printf("node : %d\n", node);
    // }



    // printf("A = \n");
    // printMatrix(A, nx*ny, 5);
    //
    // printf("b = \n");
    // printMatrix(b, nx*ny, 1);

    const double phi0 = 1;
    double* phi = (double*) malloc(nx * ny * sizeof(double*));

    std::fill_n(phi, nx*ny, phi0);

    unsigned int it = 0;
    solveSystem(nx, ny, 1e-12, 500, A, b, phi, it);

    printf("phi = \n");
    printReversedRowMatrix(phi, ny, nx);

    std::ofstream file;
    file.open("output/output.dat");
    if(file.is_open()) {
        file << std::setprecision(5) << std::fixed;
        for(unsigned int i = 0; i < nx; i++) {
            for(unsigned int j = 0; j < ny; j++)
                file << nodeX[i] << " " << nodeY[j] << " " << phi[j*nx+i] << std::endl;
            file << std::endl;
        }
    }
    file.close();

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

    printf("maxDiff : %.5f\n", maxDiff);

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

}
