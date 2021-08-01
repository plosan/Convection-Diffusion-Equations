#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <iostream>
#include <random>
#include <string>

#include "matrix.h"
#include "mesh.h"

void smithHuttonCase();
void computeDiscretizationCoefficients();
void solveSystem(const unsigned int nx, const unsigned int ny, const double tol, const double* A, const double* b, double* phi);

int main(int arg, char* argv[]) {

    double x0 = 0;
    double y0 = 0;
    unsigned int nx = 13;
    unsigned int ny = 9;
    double lx = 12;
    double ly = 8;
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

    printMeshInfo(x0, y0, nx, ny, lx, ly, lz, nodeX, nodeY, distX, distY, faceX, faceY, surfX, surfY, vol);



    double* A = (double*) malloc(5 * nx * ny * sizeof(double*));
    double* b = (double*) malloc(nx * ny * sizeof(double*));

    const double phi_low = 10;
    const double phi_high = 20;

    const double v0 = 1;
    const double alpha = 0.25 * M_PI;
    const double vx = v0 * cos(alpha);
    const double vy = v0 * sin(alpha);

    // Thermophysical properties for water at 20 ºC
    const double lambda = 0.5861;       // Thermal conductivity                         [W/(k*m)]
    const double cv = 4183;             // Specific heat at constant volume (pressure)  [J/(kg*K)]
    const double rho = 998.2;           // Density                                      [kg/m^3]
    const double gamma = lambda / cv;   // Diffusion coefficient

    // Discretization coefficients for internal nodes
    std::fill_n(A, 5*nx*ny, 0);
    std::fill_n(b, nx*ny, 0);
    for(unsigned int i = 1; i < nx-1; i++) {
        for(unsigned int j = 1; j < ny-1; j++) {
            int node = j * nx + i;
            // South node
            double massFlow = -rho * vy * surfY[i];
            double As = (massFlow + std::abs(massFlow)) / (2 * massFlow);
            A[5*node] = (gamma * surfY[i]) / distY[j-1] + massFlow * As;
            // West node
            massFlow = -rho * vx * surfX[j];
            double Aw = (massFlow + std::abs(massFlow)) / (2 * massFlow);
            A[5*node+1] = (gamma * surfX[j]) / distX[i-1] + massFlow * Aw;
            // East node
            massFlow = rho * vx * surfX[j];
            double Ae = (massFlow - std::abs(massFlow)) / (2 * massFlow);
            A[5*node+2] = (gamma * surfX[j]) / distX[i] + massFlow * Ae;
            // North node
            massFlow = rho * vy * surfY[i];
            double An = (massFlow - std::abs(massFlow)) / (2 * massFlow);
            A[5*node+3] = (gamma * surfY[i]) / distY[j] + massFlow * An;
            // Central node
            A[5*node+4] = A[5*node] + A[5*node+1] + A[5*node+2] + A[5*node+3];
        }
    }

    // Discretization coefficients for lower row nodes
    for(unsigned int i = 0; i < nx; i++) {
        A[5*i+4] = 1;
        b[i] = phi_low;
    }

    // Discretization coefficients for right column nodes
    for(unsigned int j = 1; j < ny; j++) {
        int node = (j + 1) * nx - 1;
        A[5*node+4] = 1;
        b[node] = phi_low;
    }

    // Discretization coefficients for left column nodes
    for(unsigned int j = 1; j < ny; j++) {
        int node = j * nx;
        A[5*node+4] = 1;
        b[node] = phi_high;
    }

    // Discretization coefficients for right column nodes
    for(unsigned int i = 1; i < nx-1; i++) {
        int node = (ny - 1) * nx + i;
        A[5*node+4] = 1;
        b[node] = phi_high;
    }

    const double phi0 = 1;
    double* phi = (double*) malloc(nx * ny * sizeof(double*));
    // printReversedRowMatrix(phi, ny, nx);

    std::fill_n(phi, nx*ny, phi0);
    // printReversedRowMatrix(phi, ny, nx);



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
    free(phi);

    return 0;
}


void solveSystem(const unsigned int nx, const unsigned int ny, const double tol, const double* A, const double* b, double* phi) {
    bool convergence = false;
    while(!convergence) {
        double maxDiff = -1;
        // Lower row
        double aux = phi[node];

        int node = 0;
        phi[node] = b[node]/A[5*node];
        double diff = std::abs(aux - phi[node]);
        maxDiff = (diff > maxDiff ? diff : maxDiff);

        for(int i = 1; i < nx; i++) {
            node = i;
            phi[node] = (b[node] - A[5*node])
        }

        // Central rows

        // Upper row
    }
}























adas
