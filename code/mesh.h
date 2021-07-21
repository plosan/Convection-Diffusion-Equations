#ifndef CDE_MESH_H_
#define CDE_MESH_H_

#include <iostream>
#include "matrix.h"

void compute2DUniformMesh(const unsigned int nx, const unsigned int ny, const double lx, const double ly, const double lz,
    double* nodeX, double* nodeY, double* faceX, double* faceY, double* surfX, double* surfY, double* vol);

void printMeshInfo(const unsigned int nx, const unsigned int ny, const double lx, const double ly, const double lz,
    const double* nodeX, const double* nodeY, const double* faceX, const double* faceY, const double* surfX, const double* surfY, const double* vol);

#endif // CDE_MESH_H_
