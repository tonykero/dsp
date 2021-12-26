#pragma once
#include <fftw3.h>
#include <stdio.h>
#include <math.h>

void print_complex(fftw_complex* arr, size_t len) {
    for(int i = 0; i < len; i++) {
        printf("%f + %fi\n", arr[i][0], arr[i][1]);
    }
}

void print_mat(fftw_complex* m, size_t cols, size_t rows) {
    for(size_t i = 0; i < rows; i++) {
        for(size_t j = 0; j < cols; j++) {
            printf("%f + %fi \t", m[i*cols + j][0],m[i*cols + j][1]);
        }
        printf("\n");
    }
}

bool equal(fftw_complex* arr1, fftw_complex* arr2, size_t len, double eps = 10e-5) {
    int i = 0;
    while(i < len && fabs(arr1[i][0]-arr2[i][0]) < eps && fabs(arr1[i][1]-arr2[i][1]) < eps) i++;

    if(i != len) {
        printf("err at %d\n",i);
    } else {
        printf("equal\n");
    }
    return i == len;
}