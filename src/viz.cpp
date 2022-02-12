#include "viz.hpp"

#define _USE_MATH_DEFINES
#include <math.h>

#include <fftw3.h>

void plot_complex(fftw_complex* arr, size_t len) {
    std::vector<std::vector<double>> y_data(2);
    for(int i = 0; i < len; i++) {
        double re = arr[i][0];
        double im = arr[i][1];
        y_data[0].push_back(re);
        y_data[1].push_back(im);
    }
    matplot::bar(y_data);
}

void plot_fft(fftw_complex* arr, size_t len) {
    std::vector<std::vector<double>> y_data(2);
    for(int i = 0; i < len; i++) {
        double re = arr[i][0];
        double im = arr[i][1];
        double mag = sqrt(re*re + im*im);
        double pha = atan2(im,re);
        y_data[0].push_back(mag);
        y_data[1].push_back(pha);
    }
    matplot::stem(y_data);
}

void plot_fft_signal(fftw_complex* arr, size_t len, double os) {
    std::vector<double> y_data;
    for(int i = 0; i < len*os; i++) {
        double i_ = i / os;
        double y = 0.0;
        for(int j = 0; j < len/2 + 1; j++) {
            double re = arr[j][0];
            double im = arr[j][1];
            double x = 2.0 * M_PI * j * i_ / (double)len;

            double mag = sqrt(re*re + im*im);
            double phi = atan2(im,re);
            y += mag * cos(x + phi);
        }
        y_data.push_back(y / (len/2.0));
    }
    matplot::stem(y_data);
}