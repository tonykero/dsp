#include <cstdio>
#include <fftw3.h>

#include "utils.h"
#include "viz.hpp"
#include "dft.h"

void fftw_dft(fftw_complex* in, fftw_complex* out, size_t len) {
    fftw_plan p = fftw_plan_dft_1d(len, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
}

int main() {
    const size_t N = 32;
    fftw_complex *in, *out_ref, *out;
    
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    for(int i = 0; i < N; i++) {
        in[i][0] = in[i][1] = 0.0;
    }
    in[1][0] = 5;
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out_ref = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    fftw_dft(in,out_ref,N);
    //naive_dft(in,out,N);
    naive_dft_rad2(in,out,N);

    //plot_complex(in,N);
    //plot_complex(out,N);
    //plot_fft(out,N);
    plot_fft_signal(out_ref,N);
    plot_fft_signal(out,N);

    equal(out,out_ref,N);

    fftw_free(in);
    fftw_free(out);
    fftw_free(out_ref);
    return 0;
}