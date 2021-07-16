#include <iostream>
#include <cstdlib>
#include <ctime>
#include <random>
#include <string>
#include "matrix.h"

void compute2DNodesPositionUniform(unsigned int nx, unsigned int ny, double lx, double ly, double* nodeX, double* nodeY);

int main(int arg, char* argv[]) {

    // srand((unsigned) time(0));
    //
    // int rows = 8;
    // int cols = 4;
    //
    // double* mat = (double*) malloc(rows * cols * sizeof(double*));
    // getRandomMatrix(mat, rows, cols, -50, 50);
    // printMatrix(mat, rows, cols);
    // free(mat);
    //
    // double* vec = (double*) malloc(rows*sizeof(double*));
    // getRandomMatrix(vec, rows, 1, -50, 50);
    // printMatrix(vec, rows, 1);
    //
    // printf("%.5f\n", infNorm(vec, rows));
    //
    //
    // double p = 2.45678;
    // printf("%.5f\n", pNorm(vec, rows, p));

    unsigned int nx = 5;
    unsigned int ny = 6;
    double lx = 1;
    double ly = 1;

    double* nodeX = (double*) malloc(nx * sizeof(double*));
    double* nodeY = (double*) malloc(ny * sizeof(double*));

    compute2DNodesPositionUniform(nx, ny, lx, ly, nodeX, nodeY);

    printf("Nodes X\n");
    printMatrix(nodeX, 1, nx);

    printf("Nodes Y\n");
    printMatrix(nodeY, ny, 1);


    return 0;
}

void compute2DUniformMesh(unsigned int nx, unsigned int ny, double lx, double ly, double lz, double* nodeX, double* nodeY, double* faceX, double* faceY, double* surfX, double* surfY, double* vol) {
    /*
    compute2DUniformMesh: computes the geometry of a uniform cartesian discretization for a 2D rectangular domain
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - nx: discretization nodes in X axis    [unsigned int]
        - ny: discretization nodes in Y axis    [unsigned int]
        - lx: domain length in X axis           [double]
        - ly: domain length in Y axis           [double]
        - lz: domain length in Z axis           [double]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs:
        - nodeX: discretization nodes position in X axis                    [double* - dimension nx - units of lx]
        - nodeY: discretization nodes position in Y axis                    [double* - dimension ny - units of ly]
        - faceX: position in X axis of faces perpendicular to the X axis.   [double* - dimension nx+1 - units of lx]
        - faceY: position in Y axis of faces perpendicular to the Y axis.   [double* - dimension ny+1 - units of ly]
        - surfX: surface of faces perpendicular to the X axis.              [double* - dimension ny - units of ly*lz]
        - surfY: surface of faces perpendicular to the Y axis.              [double* - dimension nx - units of lx*lz]
        - vol: volume of control volumes.                                   [double* - dimension nx*ny - units of lx*ly*lz]
    */

    // Nodes position and faces position in X axis
    nodeX[0] = 0;
    faceX[0] = 0;
    faceX[nx] = lx;
    double stepX = lx / (nx - 1);
    for(unsigned int i = 1; i < nx; ++i) {
        nodeX[i] = i*stepX;
        faceX[i] = nodeX[i] - 0.5 * stepX;
    }

    // Nodes position and faces position in Y axis
    nodeY[0] = 0;
    faceY[0] = 0;
    faceY[ny] = ly;
    double stepY = ly / (ny - 1);
    for(unsigned int j = 1; j < ny; ++j) {
        nodeY[j] = j*stepY;
        faceY[j]= nodeY[j] - 0.5*stepY;
    }

    // Faces surface in X axis
    surfX[0] = stepY / 2;
    surfX[ny-1] = stepY / 2;
    for(unsigned int j = 1; j < ny - 1; ++j)
        surfX[j] = stepY;

    // Faces surface in Y axis
    surfY[0] = stepX / 2;
    surfY[nx-1] = stepX / 2;
    for(unsigned int i = 1; i < nx - 1; ++i)
        surfY[i] = stepX;

    // Volumes
    for(unsigned int i = 0; i < nx; ++i)
        for(unsigned int j = 0; j < ny; ++j)
            vol[i*nx+j] = surfX[j] * surfY[i];
}

void compute2DNodesPositionUniform(unsigned int nx, unsigned int ny, double lx, double ly, double* nodeX, double* nodeY) {
    double stepX = lx / (nx - 1);
    for(unsigned int i = 0; i < nx; ++i)
        nodeX[i] = i*stepX;

    double stepY = ly / (ny - 1);
    for(unsigned int j = 0; j < ny; ++j)
        nodeY[j] = j*stepY;
}

void compute2DFacesUniform(unsigned int nx, unsigned int ny, double lx, double ly, double* faceX, double* faceY) {

}
