#include <iostream>
#include <cstdlib>
#include <ctime>
#include <random>
#include <string>
#include <cstring>
#include "matrix.h"
#include "mesh.h"

// void getAdjacencyList(int* list, int );

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

    unsigned int nx = 7;
    unsigned int ny = 5;
    double lx = 1;
    double ly = 1;
    double lz = 1;

    double* nodeX = (double*) malloc(nx * sizeof(double*));
    double* nodeY = (double*) malloc(ny * sizeof(double*));
    double* faceX = (double*) malloc((nx + 1) * sizeof(double*));
    double* faceY = (double*) malloc((ny + 1) * sizeof(double*));
    double* surfX = (double*) malloc(ny * sizeof(double*));
    double* surfY = (double*) malloc(nx * sizeof(double*));
    double* vol   = (double*) malloc(nx * ny * sizeof(double*));
    // compute2DNodesPositionUniform(nx, ny, lx, ly, nodeX, nodeY);

    printf("Here\n");

    compute2DUniformMesh(nx, ny, lx, ly, lz, nodeX, nodeY, faceX, faceY, surfX, surfY, vol);

    printMeshInfo(nx, ny, lx, ly, lz, nodeX, nodeY, faceX, faceY, surfX, surfY, vol);

    free(nodeX);
    free(nodeY);
    free(faceX);
    free(faceY);
    free(surfX);
    free(surfY);
    free(vol);

    int n = 6;
    int* x = (int*) malloc(n*sizeof(int*));

    printf("x = \n");
    printMatrix(x, 1, n);

    memset(x, -1, n*sizeof(int*));

    printf("x = \n");
    printMatrix(x, 1, n);



    return 0;
}
