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
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UNSAFE ACCESS TO DATA
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


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
