#include "mesh.h"


void compute2DUniformMesh(const unsigned int nx, const unsigned int ny, const double lx, const double ly, const double lz, double* nodeX, double* nodeY,
    double* faceX, double* faceY, double* surfX, double* surfY, double* vol) {
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

    nodeX[0] = 0;
    double stepX = lx / (nx - 1);
    for(unsigned int i = 1; i < nx; i++)
        nodeX[i] = i * stepX;

    nodeY[0] = 0;
    double stepY = ly / (ny - 1);
    for(unsigned int j = 1; j < ny; j++)
        nodeY[j] = j * stepY;

    faceX[0] = 0;
    for(unsigned int i = 1; i < nx; i++)
        faceX[i] = nodeX[i] - 0.5 * stepX;
    faceX[nx] = lx;

    faceY[0] = 0;
    for(unsigned int j = 1; j < ny; j++)
        faceY[j] = nodeY[j] - 0.5 * stepY;
    faceY[ny] = ly;

    surfX[0] = 0.5 * stepY * lz;
    for(unsigned int j = 1; j < ny - 1; j++)
        surfX[j] = stepY * lz;
    surfX[ny-1] = 0.5 * stepY * lz;

    surfY[0] = 0.5 * stepX * lz;
    for(unsigned int i = 1; i < nx - 1; i++)
        surfY[i] = stepX * lz;
    surfY[nx-1] = 0.5 * stepX * lz;

    for(unsigned int i = 0; i < nx; i++)
        for(unsigned int j = 0; j < ny; j++)
            vol[j * nx + i] = surfX[j] * surfY[i] * lz;

}


void printMeshInfo(const unsigned int nx, const unsigned int ny, const double lx, const double ly, const double lz, const double* nodeX, const double* nodeY, const double* faceX, const double* faceY, const double* surfX, const double* surfY, const double* vol) {
    /*
    printMeshInfo: prints the information of a 2D mesh
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs:
        - nx: number of discretization nodes in X axis                      [const unsigned int]
        - ny: number of discretization nodes in Y axis                      [const unsigned int]
        - lx: domain length in X axis                                       [const double]
        - ly: domain length in Y axis                                       [const double]
        - lz: domain length in Z axis                                       [const double]
        - nodeX: discretization nodes position in X axis                    [double* - dimension nx - units of lx]
        - nodeY: discretization nodes position in Y axis                    [double* - dimension ny - units of ly]
        - faceX: position in X axis of faces perpendicular to the X axis.   [double* - dimension nx+1 - units of lx]
        - faceY: position in Y axis of faces perpendicular to the Y axis.   [double* - dimension ny+1 - units of ly]
        - surfX: surface of faces perpendicular to the X axis.              [double* - dimension ny - units of ly*lz]
        - surfY: surface of faces perpendicular to the Y axis.              [double* - dimension nx - units of lx*lz]
        - vol: volume of control volumes.                                   [double* - dimension nx*ny - units of lx*ly*lz]
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: none
    */

    // Basic mesh information: nx, ny, lx, ly, lzs
    printf("Basic mesh information:\n");
    printf("%6s : %d\n", "nx", nx);
    printf("%6s : %d\n", "ny", ny);
    printf("%6s : %.2f\n", "lx", lx);
    printf("%6s : %.2f\n", "ly", ly);
    printf("%6s : %.2f\n\n", "lz", lz);
    // Nodes position in X axis
    printf("Nodes X\n");
    printMatrix(nodeX, 1, nx);
    // Nodes position in Y axis
    printf("Nodes Y\n");
    printReversedRowMatrix(nodeY, ny, 1);
    // Position of the faces perpendicular to the X axis
    printf("Faces X\n");
    printMatrix(faceX, 1, nx + 1);
    // Position of the faces perpendicular to the Y axis
    printf("Faces Y\n");
    printReversedRowMatrix(faceY, ny + 1, 1);
    // Surface of the faces perpendicular to the X axis
    printf("Surf X\n");
    printReversedRowMatrix(surfX, ny, 1);
    // Surface of the faces perpendicular to the Y axis
    printf("Surf Y\n");
    printMatrix(surfY, 1, nx);
    // CVs volumes
    printf("Volumes:\n");
    printReversedRowMatrix(vol, ny, nx);
}
