#include "Mesh.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructors
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Mesh::Mesh(void) {
    /*
    Mesh: constructor. Sets built to false. Sets x0, y0, lx, ly, lz, nx, ny to zero.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs: none
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: none
    */
    built = false;
    x0 = 0;
    y0 = 0;
    lx = 0;
    ly = 0;
    lz = 0;
    nx = 0;
    ny = 0;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Destructor
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Mesh::~Mesh(void) {
    /*
    Mesh: destructor. Sets built to false. Sets x0, y0, lx, ly, lz, nx, ny to zero. Frees pointers.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs: none
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: none
    */
    // Reset variables
    built = false;
    x0 = 0;
    y0 = 0;
    lx = 0;
    ly = 0;
    lz = 0;
    nx = 0;
    ny = 0;
    // Free pointers
    if(nodeX)
        free(nodeX);
    if(nodeY)
        free(nodeY);
    if(distX)
        free(distX);
    if(distY)
        free(distY);
    if(faceX)
        free(faceX);
    if(faceY)
        free(faceY);
    if(surfX)
        free(surfX);
    if(surfY)
        free(surfY);
    if(vol)
        free(vol);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GETTERS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Mesh::isBuilt(void) const {
    return built;
}

double Mesh::getX0(void) const {
    return x0;
}

double Mesh::getY0(void) const {
    return y0;
}

double Mesh::getLX(void) const {
    return lx;
}

double Mesh::getLY(void) const {
    return ly;
}

double Mesh::getLZ(void) const {
    return lz;
}

unsigned int Mesh::getNX(void) const {
    return nx;
}

unsigned int Mesh::getNY(void) const {
    return ny;
}

double* Mesh::getNodeX(void) const {
    return nodeX;
}

double* Mesh::getNodeY(void) const {
    return nodeY;
}

double* Mesh::getDistX(void) const {
    return distX;
}

double* Mesh::getDistY(void) const {
    return distY;
}

double* Mesh::getFaceX(void) const {
    return faceX;
}

double* Mesh::getFaceY(void) const {
    return faceY;
}

double* Mesh::getSurfX(void) const {
    return surfX;
}

double* Mesh::getSurfY(void) const {
    return surfY;
}

double* Mesh::getVol(void) const {
    return vol;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SAFE ACCESS TO DATA
// These functions check whether the mesh is built and the argument/s give/n are within the range. If the argument/s is/are in the correct range,
// the function returns the specified value. Otherwise it returns a value which should not cause issues, for instance, if a distance is requested,
// it returns 1 instead of 0.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Mesh::satNodeX(unsigned int i) const {
    if(built) { // Mesh is built
        if(i < nx)      // Safe range: 0 <= i < nx
            return nodeX[i];
        else            // Unsuitable argument
            printf("\tError accessing nodeX. Argument provided (%d) is not within the range 0 <= i < %d\n", i, nx);
    } else      // Mesh is not built
        printf("\tError accessing nodeX. Mesh is not built\n");
    return 0;
}

double Mesh::satNodeY(unsigned int j) const {
    if(built) { // Mesh is built
        if(j < ny)  // Safe range: 0 <= j < ny
            return nodeY[j];
        else        // Unsuitable argument
            printf("\tError accessing nodeY. Argument provided (%d) is not within the range 0 <= j < %d\n", j, ny);
    } else      // Mesh is not built
        printf("\tError accessing nodeY. Mesh is not built\n");
    return 0;
}

double Mesh::satDistX(unsigned int i) const {
    if(built) { // Mesh is built
        if(i < nx-1)    // Safe range: 0 <= i < nx-1
            return distX[i];
        else            // Unsuitable argument
            printf("\tError accessing distX. Argument provided (%d) is not within the range 0 <= i < %d\n", i, nx-1);
    } else      // Mesh is not built
        printf("\tError accessing distX. Mesh is not built\n");
    return 1;
}

double Mesh::satDistY(unsigned int j) const {
    if(built) { // Mesh is built
        if(j < ny-1)    // Safe range: 0 <= j < ny-1
            return distY[j];
        else            // Unsuitable argument
            printf("\tError accessing distY. Argument provided (%d) is not within the range 0 <= j < %d\n", j, ny-1);
    } else      // Mesh is not built
        printf("\tError accessing distY. Mesh is not built\n");
    return 1;
}

double Mesh::satFaceX(unsigned int i) const {
    if(built) { // Mesh is built
        if(i < nx+1)    // Safe range: 0 <= i < nx+1
            return faceX[i];
        else            // Unsuitable argument
            printf("\tError accessing faceX. Argument provided (%d) is not within the range 0 <= i < %d\n", i, nx+1);
    } else      // Mesh is not built
        printf("\tError accessing faceX. Mesh is not built\n");
    return 1;
}

double Mesh::satFaceY(unsigned int j) const {
    if(built) { // Mesh is built
        if(j < ny+1)    // Safe range: 0 <= j < ny+1
            return faceY[j];
        else            // Unsuitable argument
            printf("\tError accessing faceY. Argument provided (%d) is not within the range 0 <= j < %d\n", j, ny+1);
    } else      // Mesh is not built
        printf("\tError accessing faceY. Mesh is not built\n");
    return 1;
}

double Mesh::satSurfX(unsigned int j) const {
    if(built) { // Mesh is built
        if(j < ny)  // Safe range: 0 <= j < ny
            return surfX[j];
        else        // Unsuitable argument
            printf("\tError accessing surfX. Argument provided (%d) is not within the range 0 <= j < %d\n", j, ny);
    } else      // Mesh is not built
        printf("\tError accessing surfX. Mesh is not built\n");
    return 1;
}

double Mesh::satSurfY(unsigned int i) const {
    if(built) { // Mesh is built
        if(i < nx)  // Safe range: 0 <= i < nx
            return surfY[i];
        else        // Unsuitable argument
            printf("\tError accessing surfY. Argument provided (%d) is not within the range 0 <= i < %d\n", i, nx);
    } else      // Mesh is not built
        printf("\tError accessing surfY. Mesh is not built\n");
    return 1;
}

double Mesh::satVol(unsigned int i) const {
    if(built) { // Mesh is built
        if(i < nx*ny)   // Safe range: 0 <= i < nx*ny
            return vol[i];
        else            // Unsuitable argument
            printf("\tError accessing vol. Argument provided (%d) is not within the range 0 <= i < %d\n", i, nx*ny);
    } else      // Mesh is not built
        printf("\tError accessing vol. Mesh is not built\n");
    return 1;
}

double Mesh::satVol(unsigned int i, unsigned int j) const {
    if(built) { // Mesh is built
        if(i < nx && j < ny)    // Safe range: 0 <= i < nx; 0 <= j < ny
            return vol[j*nx+i];
        else                    // Unsuitable arguments
            printf("\tError accessing vol. Arguments provided (i = %d, j = %d) are not within the ranges 0 <= i < %d, 0 <= j < %d\n", i, j, nx, ny);
    } else      // Mesh is not built
        printf("\tError accessing vol. Mesh is not built\n");
    return 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UNSAFE ACCESS TO DATA
// These functions do not check if the argument/s give/n are within the range.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Mesh::atNodeX(int i) const {
    return nodeX[i];
}

double Mesh::atNodeY(int j) const {
    return nodeY[j];
}

double Mesh::atDistX(int i) const {
    return distX[i];
}

double Mesh::atDistY(int j) const {
    return distY[j];
}

double Mesh::atFaceX(int i) const {
    return faceX[i];
}

double Mesh::atFaceY(int j) const {
    return faceY[j];
}

double Mesh::atSurfX(int j) const {
    return surfX[j];
}

double Mesh::atSurfY(int i) const {
    return surfY[i];
}

double Mesh::atVol(int i) const {
    return vol[i];
}

double Mesh::atVol(int i, int j) const {
    return vol[j*nx+i];
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// BUILD MESH
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Mesh::buildUniformMesh(double _x0, double _y0, double _lx, double _ly, double _lz, unsigned int _nx, unsigned int _ny) {
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

    // nodeX[0] = 0;
    // double stepX = lx / (nx - 1);
    // for(unsigned int i = 1; i < nx; i++)
    //     nodeX[i] = i * stepX;
    //
    // nodeY[0] = 0;
    // double stepY = ly / (ny - 1);
    // for(unsigned int j = 1; j < ny; j++)
    //     nodeY[j] = j * stepY;
    //
    // faceX[0] = 0;
    // for(unsigned int i = 1; i < nx; i++)
    //     faceX[i] = nodeX[i] - 0.5 * stepX;
    // faceX[nx] = lx;
    //
    // faceY[0] = 0;
    // for(unsigned int j = 1; j < ny; j++)
    //     faceY[j] = nodeY[j] - 0.5 * stepY;
    // faceY[ny] = ly;
    //
    // surfX[0] = 0.5 * stepY * lz;
    // for(unsigned int j = 1; j < ny - 1; j++)
    //     surfX[j] = stepY * lz;
    // surfX[ny-1] = 0.5 * stepY * lz;
    //
    // surfY[0] = 0.5 * stepX * lz;
    // for(unsigned int i = 1; i < nx - 1; i++)
    //     surfY[i] = stepX * lz;
    // surfY[nx-1] = 0.5 * stepX * lz;
    //
    // for(unsigned int i = 0; i < nx; i++)
    //     for(unsigned int j = 0; j < ny; j++)
    //         vol[j * nx + i] = surfX[j] * surfY[i] * lz;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OTHER FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Mesh::printMeshInfo(void) const {
    /*
    printMeshInfo: prints the information of a 2D mesh
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs: node
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
