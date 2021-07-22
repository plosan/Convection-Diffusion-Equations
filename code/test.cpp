#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <cstring>
#include "matrix.h"


int main(void) {

    unsigned int nx = 7;
    unsigned int ny = 5;
    int* mat = (int*) malloc(nx * ny * sizeof(int*));

    if(mat) {
        int counter = 0;
        for(int j = 0; j < ny; j++)
            for(int i = 0; i < nx; i++)
                mat[j*nx+i] = counter++;

        printf("%10s : %d\n", "nx", nx);
        printf("%10s : %d\n\n", "ny", ny);

        printf("mat = \n");
        printReversedRowMatrix(mat, ny, nx);

        free(mat);
    }
}
