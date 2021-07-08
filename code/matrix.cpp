#include "matrix.h"

void printMatrix(double* mat, unsigned int rows, unsigned int cols) {
    if(mat) {
        for(unsigned int i = 0; i < rows; ++i) {
            for(unsigned int j = 0; j < cols; ++j)
                printf("%10.2f", mat[i*cols+j]);
            printf("\n");
        }
        printf("\n");
    }
}

void getRandomMatrix(double* mat, unsigned int rows, unsigned int cols, int min_range, int max_range) {
    if(mat) {
        for(unsigned int i = 0; i < rows; ++i)
            for(unsigned int j = 0; j < cols; ++j)
                mat[i*cols+j] = rand()%(max_range - min_range + 1) + min_range;
    }
}

void getSDDMatrix(double* mat, unsigned int rows, int min_range, int max_range) {
    if(mat) {
        for(unsigned int i = 0; i < rows; ++i) {
            int sum = 0;
            for(unsigned int j = 0; j < rows; ++j) {
                mat[i*rows+j] = rand()%(max_range - min_range + 1) + min_range;
                sum += abs(mat[i*rows+j]);
            }
            mat[i*(rows+1)] = (rand()%2 == 0 ? 1 : -1)*sum;
        }
    }
}

void gaussSeidel(double* mat, double* vec, double* sol, unsigned int rows, double& error) {
    if(mat && vec && sol) {

    } else {
        error = -1;
    }
}

double pNorm(double* vec, unsigned int size, double p) {
    if(!vec)
        return -1;

    if(p < 1)
        return -1;

    double norm = 0;
    for(unsigned int i = 0; i < size; ++i)
        norm += pow(abs(vec[i]), p);
    return pow(norm, 1/p);
}

double infNorm(double* vec, unsigned int size) {
    if(!vec)
        return -1;

    double norm = -1;
    for(unsigned int i = 0; i < size; ++i)
        norm = (norm > abs(vec[i]) ? norm : abs(vec[i]));
    return norm;
}
