#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <cstring>
// #include "matrix.h"

int n;
int* mat;

int at(int i) {
    return mat[i];
}

int add(int x, int y) {
    return x + y;
}

void printFunction(int x, int y, int (*f)(int, int)) {
    int z = (*f)(x,y);
    printf("%d + %d = %d\n", x, y, z);
}


int main(void) {

    // printFunction(9,6,add);

    n = 5;
    mat = (int*) malloc(n * sizeof(int*));
    for(int i = 0; i < n; i++)
        mat[i] = i;

    for(int i = 0; i < n; i++)
        printf("%5d", mat[i]);
    printf("\n");

    for(int i = 0; i < n; i++)
        printf("mat[%d] : %d\n", i, at(i));
    printf("\n");

    return 0;
}
