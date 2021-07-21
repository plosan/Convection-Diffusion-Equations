#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <cstring>
#include "matrix.h"


int main(void) {

    int n = 6;
    int* x = (int*) malloc(n*sizeof(int*));

    printf("x = \n");
    printMatrix(x, 1, n);

    memset(x, -1, n*sizeof(int*));

    printf("x = \n");
    printMatrix(x, 1, n);
}
