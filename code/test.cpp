#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <cfloat>
#include <cstring>
// #include "matrix.h"


int add(int x, int y) {
    return x + y;
}

void printFunction(int x, int y, int (*f)(int, int)) {
    int z = (*f)(x,y);
    printf("%d + %d = %d\n", x, y, z);
}


int main(void) {

    printFunction(9,6,add);

}
