#include <cstdio>
#include <fftw3.h>
#include "viz.hpp"


void print_complex(fftw_complex* arr, size_t len) {
    for(int i = 0; i < len; i++) {
        printf("%f + %fi\n", arr[i][0], arr[i][1]);
    }
}


void fftw_dft(fftw_complex* in, fftw_complex* out, size_t len) {
    fftw_plan p = fftw_plan_dft_1d(len, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
}

bool equal(fftw_complex* arr1, fftw_complex* arr2, size_t len, double eps = 10e-5) {
    int i = 0;
    while(i < len && abs(arr1[i][0]-arr2[i][0]) < eps && abs(arr1[i][1]-arr2[i][1]) < eps) i++;

    if(i != len) {
        printf("err at %d\n",i);
    } else {
        printf("equal\n");
    }
    return i == len;
}

int main() {
    const size_t N = 32;
    fftw_complex *in, *out;
    
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    for(int i = 0; i < N; i++) {
        in[i][0] = in[i][1] = 0.0;
    }
    in[1][0] = 5;
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_dft(in,out,N);

    print_complex(out,N);

    //plot_complex(in,N);
    //plot_complex(out,N);
    //plot_fft(out,N);
    plot_fft_signal(out,N);

    equal(out,out,N);
    
    fftw_free(in);
    fftw_free(out);
    return 0;
}