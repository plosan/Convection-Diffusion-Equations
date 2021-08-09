#ifndef MESH_H_
#define MESH_H_

#include <cstdlib>
#include <cstring>
#include <iostream>

#include "matrix.h"

class Mesh {

private:
    bool built;         // Gives whether the mesh object is built (true) or not (false)
    double x0;          // X coordinate of the lower left corner
    double y0;          // Y coordinate of the lower left corner
    double lx;          // Domain size in the X axis
    double ly;          // Domain size in the Y axis
    double lz;          // Domain size in the Z axis
    int nx;             // Number of nodes in the X axis
    int ny;             // Number of nodes in the Y axis
    double* nodeX;      // Nodes position in the X axis. Size: nx
    double* nodeY;      // Nodes position in the Y axis. Size: ny
    double* distX;      // Distance between nodes in the X axis. Size: nx - 1
    double* distY;      // Distance between nodes in the Y axis. Size: ny - 1
    double* distNFX;    // Distance between nodes and faces in the X axis. Size: 2*nx
    double* distNFY;    // Distance between nodes and faces in the Y axis. Size: 2*ny
    double* faceX;      // Faces position in the X axis. Size: nx + 1
    double* faceY;      // Faces position in the Y axis. Size: ny + 1
    double* surfX;      // Surface of the faces perpendicular to the X axis. Size: ny
    double* surfY;      // Surface of the faces perpendicular to the Y axis. Size: nx
    double* vol;        // Volume of each control volume. Size: nx * ny

public:

    // Constructor
    Mesh();

    // Destructor
    // ~Mesh();

    // Setters
    // void setInitialCoordinates(double _x0, double _y0);

    // Getters
    bool isBuilt(void) const;           // Returns built
    double getX0(void) const;           // Returns x0
    double getY0(void) const;           // Returns y0
    double getLX(void) const;           // Returns lx
    double getLY(void) const;           // Returns ly
    double getLZ(void) const;           // Returns lz
    int getNX(void) const;     // Returns nx
    int getNY(void) const;     // Returns ny
    double* getNodeX(void) const;       // Returns nodeX
    double* getNodeY(void) const;       // Returns nodeY
    double* getDistX(void) const;       // Returns distX
    double* getDistY(void) const;       // Returns distY
    double* getDistNFX(void) const;     // Returns distNFX
    double* getDistNFY(void) const;     // Returns distNFY
    double* getFaceX(void) const;       // Returns faceX
    double* getFaceY(void) const;       // Returns faceY
    double* getSurfX(void) const;       // Returns surfX
    double* getSurfY(void) const;       // Returns surfY
    double* getVol(void) const;         // Returns vol

    // Safe access to pointers
    double satNodeX(int) const;                // Returns nodeX[i], where i is the argument passed
    double satNodeY(int) const;                // Returns nodeY[j], where j is the argument passed
    double satDistX(int) const;                // Returns distX[i], where i is the argument passed
    double satDistY(int) const;                // Returns distY[j], where j is the argument passed
    double satDistNFX(int) const;              // Returns distNFX[i], where i is the argument passed
    double satDistNFY(int) const;              // Returns distNFY[j], where j is the argument passed
    double satFaceX(int) const;                // Returns faceX[i], where i is the argument passed
    double satFaceY(int) const;                // Returns faceY[j], where j is the argument passed
    double satSurfX(int) const;                // Returns surfX[i], where i is the argument passed
    double satSurfY(int) const;                // Returns surfY[j], where j is the argument passed
    double satVol(int) const;                  // Returns vol[i], where i is the argument passed
    double satVol(int, int) const;    // Returns vol[j*nx+i], where (i,j) are the arguments passed


    // Unsafe access to pointers
    double atNodeX(int) const;
    double atNodeY(int) const;
    double atDistX(int) const;
    double atDistY(int) const;
    double atDistNFX(int) const;
    double atDistNFY(int) const;
    double atFaceX(int) const;
    double atFaceY(int) const;
    double atSurfX(int) const;
    double atSurfY(int) const;
    double atVol(int) const;
    double atVol(int, int) const;

    // Build mesh
    int buildUniformMesh(double _x0, double _y0, double _lx, double _ly, double _lz, int _nx, int _ny);

    // Reset mesh
    void resetMesh(void);

    // Other functions
    void printInfo(void) const;

};

#endif // MESH_H_
