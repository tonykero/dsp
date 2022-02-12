#pragma once

#include <fftw3.h>
#include <matplot/matplot.h>

void plot_complex(fftw_complex* arr, size_t len);
void plot_fft(fftw_complex* arr, size_t len);
void plot_fft_signal(fftw_complex* arr, size_t len, double os = 4.0);