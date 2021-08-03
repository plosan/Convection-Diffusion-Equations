#include "meshC.h"

void compute2DUniformRectangularMesh(const double x0, const double y0, const unsigned int nx, const unsigned int ny, const double lx, const double ly,
const double lz, double* nodeX, double* nodeY, double* distX, double* distY, double* faceX, double* faceY, double* surfX, double* surfY, double* vol) {

    double stepX = lx / (nx - 1);
    for(unsigned int i = 0; i < nx; i++)
        nodeX[i] = x0 + i * stepX;

    double stepY = ly / (ny - 1);
    for(unsigned int j = 0; j < ny; j++)
        nodeY[j] = y0 + j * stepY;

    for(unsigned int i = 0; i < nx-1; i++)
        distX[i] = stepX;

    for(unsigned int j = 0; j < ny-1; j++)
        distY[j] = stepY;

    faceX[0] = x0;
    faceX[nx] = x0 + lx;
    for(unsigned int i = 0; i < nx-1; i++)
        faceX[i+1] = nodeX[i] + 0.5 * stepX;

    faceY[0] = y0;
    faceY[ny] = y0 + ly;
    for(unsigned int j = 0; j < ny-1; j++)
        faceY[j+1] = nodeY[j] + 0.5 * stepY;

    surfX[0] = 0.5 * stepY * lz;
    surfX[ny-1] = surfX[0];
    for(unsigned int j = 1; j < ny-1; j++)
        surfX[j] = stepY * lz;

    surfY[0] = 0.5 * stepX * lz;
    surfY[nx-1] = 0.5 * stepX * lz;
    for(unsigned int i = 1; i < nx-1; i++)
        surfY[i] = stepX * lz;

    for(unsigned int i = 0; i < nx; i++)
        for(unsigned int j = 0; j < ny; j++)
            vol[j*nx+i] = surfX[j] * surfY[i] / lz;

}

void compute2DUniformMesh(const unsigned int nx, const unsigned int ny, const double lx, const double ly, const double lz, double* nodeX, double* nodeY, double* faceX, double* faceY, double* surfX, double* surfY, double* vol) {
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

void computeAdjacencyList(int* list, const unsigned int nx, const unsigned int ny) {
    if(list) {
        /*
        Corner nodes: 0, nx-1, nx*(ny-1), nx*ny-1
        Lower row nodes: 1 to nx-2
        Upper row nodes: nx*(ny-1)+1 to nx*ny-2
        Left column nodes: 0, nx, 2*nx, ..., nx*(ny-2)
        Right column nodes: nx-1, 2*nx-1, 3*nx-1, ..., (ny-2)*nx-1
        A node is identified by its (i,j) coordinates, i referring to the column and j to the row. Since all matrices are treated as vectors,
        the index of node (i,j) in the vector is j*nx+i.
        Given the (i,j) node, its neighbours are:
        - East node: (i+1,j) --> j*nx+i+1
        - West node: (i-1,j) --> j*nx+i-1
        - North node: (i,j+1) --> (j+1)*nx+i
        - South node: (i,j-1) --> (j-1)*nx+i
        */
        const unsigned int n = 5 * nx * ny;
        memset(list, -1, n*sizeof(int*));

    }
}

void printMeshInfo(const double x0, const double y0, const unsigned int nx, const unsigned int ny, const double lx, const double ly, const double lz,
    const double* nodeX, const double* nodeY, const double* distX, const double* distY, const double* faceX, const double* faceY,
    const double* surfX, const double* surfY, const double* vol) {
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
    printf("%6s : %.5f\n", "x0", x0);
    printf("%6s : %.5f\n", "y0", y0);
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
    // Distances between nodes in X axis
    printf("Distance X\n");
    printMatrix(distX, 1, nx-1);
    // Distances between nodes in Y axis
    printf("Distance Y\n");
    printReversedRowMatrix(distY, ny-1, 1);
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
