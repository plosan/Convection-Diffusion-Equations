#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <cstring>

#include "Mesh.h"
#include "matrix.h"

int n;
int* mat;

int at(int i) {
    return mat[i];
}

int add(int x, int y) {
    return x + y;
}

int prod(int x, int y) {
    return x*y;
}

void printFunction(int x, int y, int (*f)(int, int), int (*g)(int,int)) {
    int z = (*f)(x,y);
    int p = (*g)(x,y);
    printf("%d + %d = %d\n", x, y, z);
    printf("%d * %d = %d\n", x, y, p);
}


int main(void) {

    // printFunction(9,6,add);

    // Mesh m;
    // m.printInfo();
    // m.buildUniformMesh(0, 0, 1, 1, 1, 5, 5);
    // m.printInfo();

    printFunction(9, 9, add, prod);

    return 0;
}
