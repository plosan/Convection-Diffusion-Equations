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

// Mesh::~Mesh(void) {
//     /*
//     Mesh: destructor. Only executes if the mesh is not built. Sets built to false. Sets x0, y0, lx, ly, lz, nx, ny to zero. Frees pointers.
//     --------------------------------------------------------------------------------------------------------------------------------------------------
//     Inputs: none
//     --------------------------------------------------------------------------------------------------------------------------------------------------
//     Outputs: none
//     */
//     printf("Deleting mesh object...\n");
//     // Set member variables (double and unsigned int) to zero
//     x0 = 0;
//     y0 = 0;
//     lx = 0;
//     ly = 0;
//     lz = 0;
//     nx = 0;
//     ny = 0;
//     // Free memory
//     if(nodeX)
//         free(nodeX);
//     if(nodeY)
//         free(nodeY);
//     if(distX)
//         free(distX);
//     if(distY)
//         free(distY);
//     if(distNFX)
//         free(distNFX);
//     if(distNFY)
//         free(distNFY);
//     if(faceX)
//         free(faceX);
//     if(faceY)
//         free(faceY);
//     if(surfX)
//         free(surfX);
//     if(surfY)
//         free(surfY);
//     if(vol)
//         free(vol);
//     // Set built to false
//     built = false;
// }

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

double* Mesh::getDistNFX(void) const {
    return distNFX;
}

double* Mesh::getDistNFY(void) const {
    return distNFY;
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
// SAFE ACCESS TO POINTERS
// These functions check whether the mesh is built and the argument/s give/n are within the range. If the argument/s is/are in the correct range,
// the function returns the specified value. Otherwise it returns a value which should not cause issues, for instance, if a distance is requested,
// it returns 1 instead of 0.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Mesh::satNodeX(unsigned int i) const {
    if(built) { // Mesh is built
        if(i < nx)      // Safe range: 0 <= i < nx
            return nodeX[i];
        else            // Unsuitable argument
            printf("\tError accessing nodeX. Argument provided (%d) is not within the range 0 <= i < %d.\n", i, nx);
    } else      // Mesh is not built
        printf("\tError accessing nodeX. Mesh is not built.\n");
    return 0;
}

double Mesh::satNodeY(unsigned int j) const {
    if(built) { // Mesh is built
        if(j < ny)  // Safe range: 0 <= j < ny
            return nodeY[j];
        else        // Unsuitable argument
            printf("\tError accessing nodeY. Argument provided (%d) is not within the range 0 <= j < %d.\n", j, ny);
    } else      // Mesh is not built
        printf("\tError accessing nodeY. Mesh is not built.\n");
    return 0;
}

double Mesh::satDistX(unsigned int i) const {
    if(built) { // Mesh is built
        if(i < nx-1)    // Safe range: 0 <= i < nx-1
            return distX[i];
        else            // Unsuitable argument
            printf("\tError accessing distX. Argument provided (%d) is not within the range 0 <= i < %d.\n", i, nx-1);
    } else      // Mesh is not built
        printf("\tError accessing distX. Mesh is not built.\n");
    return 1;
}

double Mesh::satDistY(unsigned int j) const {
    if(built) { // Mesh is built
        if(j < ny-1)    // Safe range: 0 <= j < ny-1
            return distY[j];
        else            // Unsuitable argument
            printf("\tError accessing distY. Argument provided (%d) is not within the range 0 <= j < %d.\n", j, ny-1);
    } else      // Mesh is not built
        printf("\tError accessing distY. Mesh is not built.\n");
    return 1;
}

double Mesh::satDistNFX(unsigned int i) const {
    if(built) { // Mesh is built
        if(i < 2*nx)    // Safe range: 0 <= i < 2*nx
            return distNFX[i];
        else            // Unsuitable argument
            printf("\tError accessing distNFX. Argument provided (%d) is not within the range 0 <= i < %d.\n", i, 2*nx);
    } else      // Mesh is not built
        printf("\tError accessing distNFX. Mesh is not built.\n");
    return 1;
}

double Mesh::satDistNFY(unsigned int j) const {
    if(built) { // Mesh is built
        if(j < 2*ny)    // Safe range: 0 <= j < 2*ny
            return distNFY[j];
        else            // Unsuitable argument
            printf("\tError accessing distNFY. Argument provided (%d) is not within the range 0 <= j < %d.\n", j, 2*ny);
    } else      // Mesh is not built
        printf("\tError accessing distNFY. Mesh is not built.\n");
    return 1;
}

double Mesh::satFaceX(unsigned int i) const {
    if(built) { // Mesh is built
        if(i < nx+1)    // Safe range: 0 <= i < nx+1
            return faceX[i];
        else            // Unsuitable argument
            printf("\tError accessing faceX. Argument provided (%d) is not within the range 0 <= i < %d.\n", i, nx+1);
    } else      // Mesh is not built
        printf("\tError accessing faceX. Mesh is not built.\n");
    return 1;
}

double Mesh::satFaceY(unsigned int j) const {
    if(built) { // Mesh is built
        if(j < ny+1)    // Safe range: 0 <= j < ny+1
            return faceY[j];
        else            // Unsuitable argument
            printf("\tError accessing faceY. Argument provided (%d) is not within the range 0 <= j < %d.\n", j, ny+1);
    } else      // Mesh is not built
        printf("\tError accessing faceY. Mesh is not built.\n");
    return 1;
}

double Mesh::satSurfX(unsigned int j) const {
    if(built) { // Mesh is built
        if(j < ny)  // Safe range: 0 <= j < ny
            return surfX[j];
        else        // Unsuitable argument
            printf("\tError accessing surfX. Argument provided (%d) is not within the range 0 <= j < %d.\n", j, ny);
    } else      // Mesh is not built
        printf("\tError accessing surfX. Mesh is not built.\n");
    return 1;
}

double Mesh::satSurfY(unsigned int i) const {
    if(built) { // Mesh is built
        if(i < nx)  // Safe range: 0 <= i < nx
            return surfY[i];
        else        // Unsuitable argument
            printf("\tError accessing surfY. Argument provided (%d) is not within the range 0 <= i < %d.\n", i, nx);
    } else      // Mesh is not built
        printf("\tError accessing surfY. Mesh is not built.\n");
    return 1;
}

double Mesh::satVol(unsigned int i) const {
    if(built) { // Mesh is built
        if(i < nx*ny)   // Safe range: 0 <= i < nx*ny
            return vol[i];
        else            // Unsuitable argument
            printf("\tError accessing vol. Argument provided (%d) is not within the range 0 <= i < %d.\n", i, nx*ny);
    } else      // Mesh is not built
        printf("\tError accessing vol. Mesh is not built.\n");
    return 1;
}

double Mesh::satVol(unsigned int i, unsigned int j) const {
    if(built) { // Mesh is built
        if(i < nx && j < ny)    // Safe range: 0 <= i < nx; 0 <= j < ny
            return vol[j*nx+i];
        else                    // Unsuitable arguments
            printf("\tError accessing vol. Arguments provided (i = %d, j = %d) are not within the ranges 0 <= i < %d, 0 <= j < %d.\n", i, j, nx, ny);
    } else      // Mesh is not built
        printf("\tError accessing vol. Mesh is not built.\n");
    return 1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// UNSAFE ACCESS TO POINTERS
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

double Mesh::atDistNFX(int i) const {
    return distNFX[i];
}

double Mesh::atDistNFY(int j) const {
    return distNFY[j];
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
int Mesh::buildUniformMesh(double _x0, double _y0, double _lx, double _ly, double _lz, unsigned int _nx, unsigned int _ny) {
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
    if(!built) {
        // Check if the provided arguments are correct
        if(_lx < 0) {
            printf("\tError building uniform mesh. The length provided _lx (%.5f) must be positive.\n", _lx);
            return -1;
        }
        if(_ly < 0) {
            printf("\tError building uniform mesh. The length provided _ly (%.5f) must be positive.\n", _ly);
            return -1;
        }
        if(_lz < 0) {
            printf("\tError building uniform mesh. The length provided _lz (%.5f) must be positive.\n", _lz);
            return -1;
        }
        if(_nx <= 1) {
            printf("\tError building uniform mesh. The number of nodes provided _nx (%d) must be strictly greater than 1.\n", _nx);
            return -1;
        }
        if(_ny <= 1) {
            printf("\tError building uniform mesh. The number of nodes provided _nY (%d) must be strictly greater than 1.\n", _ny);
            return -1;
        }
        // Provided arguments are correct. Now set member variables
        x0 = _x0;
        y0 = _y0;
        lx = _lx;
        ly = _ly;
        lz = _lz;
        nx = _nx;
        ny = _ny;
        // Allocating memory for pointers

        // nodeX
        double stepX = lx / (nx - 1);
        nodeX = (double*) malloc(nx * sizeof(double*));
        if(nodeX) {
            for(unsigned int i = 0; i < nx; i++)
                nodeX[i] = x0 + i * stepX;
        } else {
            printf("\tError building uniform mesh. Could not allocate memory for nodeX.\n");
            return -2;
        }

        // nodeY
        double stepY = ly / (ny - 1);
        nodeY = (double*) malloc(ny * sizeof(double*));
        if(nodeY) {
            for(unsigned int j = 0; j < ny; j++)
                nodeY[j] = y0 + j * stepY;
        } else {
            printf("\tError building uniform mesh. Could not allocate memory for nodeY.\n");
            return -2;
        }

        // distX
        distX = (double*) malloc((nx-1) * sizeof(double*));
        if(distX) {
            for(unsigned int i = 0; i < nx-1; i++)
                distX[i] = stepX;
        } else {
            printf("\tError building uniform mesh. Could not allocate memory for distX.\n");
            return -2;
        }

        // distY
        distY = (double*) malloc((ny-1) * sizeof(double*));
        if(distY) {
            for(unsigned int j = 0; j < ny-1; j++)
                distY[j] = stepY;
        } else {
            printf("\tError building uniform mesh. Could not allocate memory for distY.\n");
            return -2;
        }

        // faceX
        faceX = (double*) malloc((nx+1) * sizeof(double*));
        if(faceX) {
            faceX[0] = x0;
            faceX[nx] = x0 + lx;
            for(unsigned int i = 0; i < nx-1; i++)
                faceX[i+1] = nodeX[i] + 0.5 * stepX;
        } else {
            printf("\tError building uniform mesh. Could not allocate memory for faceX.\n");
            return -2;
        }

        // faceY
        faceY = (double*) malloc((ny+1) * sizeof(double*));
        if(faceY) {
            faceY[0] = y0;
            faceY[ny] = y0 + ly;
            for(unsigned int j = 0; j < ny-1; j++)
                faceY[j+1] = nodeY[j] + 0.5 * stepY;
        } else {
            printf("\tError building uniform mesh. Could not allocate memory for faceY.\n");
            return -2;
        }

        // distNFX
        distNFX = (double*) malloc(2 * nx * sizeof(double*));
        if(distNFX) {
            for(unsigned int i = 0; i < nx; i++) {
                distNFX[2*i] = nodeX[i] - faceX[i];
                distNFX[2*i+1] = faceX[i+1] - nodeX[i];
            }
        } else {
            printf("\tError building uniform mesh. Could not allocate memory for distNFX.\n");
            return -2;
        }

        // distNFY
        distNFY = (double*) malloc(2 * ny * sizeof(double*));
        if(distNFY) {
            for(unsigned int j = 0; j < ny; j++) {
                distNFY[2*j] = nodeY[j] - faceY[j];
                distNFY[2*j+1] = faceY[j+1] - nodeY[j];
            }
        } else {
            printf("\tError building uniform mesh. Could not allocate memory for distNYX.\n");
            return -2;
        }

        // surfX
        surfX = (double*) malloc(ny * sizeof(double*));
        if(surfX) {
            surfX[0] = 0.5 * stepY * lz;
            surfX[ny-1] = surfX[0];
            for(unsigned int j = 1; j < ny-1; j++)
                surfX[j] = stepY * lz;
        } else {
            printf("\tError building uniform mesh. Could not allocate memory for surfX.\n");
            return -2;
        }

        // surfY
        surfY = (double*) malloc(nx * sizeof(double*));
        if(surfY) {
            surfY[0] = 0.5 * stepX * lz;
            surfY[nx-1] = 0.5 * stepX * lz;
            for(unsigned int i = 1; i < nx-1; i++)
                surfY[i] = stepX * lz;
        } else {
            printf("\tError building uniform mesh. Could not allocate memory for surfY.\n");
            return -2;
        }

        // vol
        vol = (double*) malloc(nx * ny * sizeof(double*));
        if(vol) {
            for(unsigned int i = 0; i < nx; i++)
                for(unsigned int j = 0; j < ny; j++)
                    vol[j*nx+i] = surfX[j] * surfY[i] / lz;
        } else {
            printf("\tError building uniform mesh. Could not allocate memory for vol.\n");
            return -2;
        }
        // The uniform mesh has been built successfully
        built = true;
        return 1;
    } else
        printf("\tError building uniform mesh. Mesh already constructed.\n");
    return -1;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RESET MESH
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Mesh::resetMesh(void) {
    /*
    resetMesh: sets all member variables to false/zero and frees the memory previously allocated for pointers.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs: none
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: none
    */
    printf("Reset mesh...\n");
    // Set member variables (double and unsigned int) to zero
    x0 = 0;
    y0 = 0;
    lx = 0;
    ly = 0;
    lz = 0;
    nx = 0;
    ny = 0;
    // Free memory
    if(nodeX)
        free(nodeX);
    if(nodeY)
        free(nodeY);
    if(distX)
        free(distX);
    if(distY)
        free(distY);
    if(distNFX)
        free(distNFX);
    if(distNFY)
        free(distNFY);
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
    // Set built to false
    built = false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OTHER FUNCTIONS
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Mesh::printInfo(void) const {
    /*
    printMeshInfo: prints the mesh information if it is built. Otherwise prints an error message.
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Inputs: node
    --------------------------------------------------------------------------------------------------------------------------------------------------
    Outputs: none
    */
    if(built) {
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
        // Distance between adjacent nodes and faces in the X axis
        printf("Node-Face distance X\n");
        printMatrix(distNFX, 1, 2*nx);
        // Distance between adjacent nodes and faces in the Y axis
        printf("Node-Face distance Y\n");
        printReversedRowMatrix(distNFY, 2*ny, 1);
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
    } else
        printf("Error printing mesh information. Mesh is not built.\n");
}
