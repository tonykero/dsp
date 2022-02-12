#pragma once
#include <fftw3.h>
#include <stdio.h>
#include <cmath>

void print_complex(fftw_complex* arr, size_t len);
void print_mat(fftw_complex* m, size_t cols, size_t rows);
bool equal(fftw_complex* arr1, fftw_complex* arr2, size_t len, double eps = 10e-5);