#pragma once

#include <stdint.h>
#include <fftw3.h>

fftw_complex* cadd(fftw_complex* a, fftw_complex* b);
fftw_complex* csub(fftw_complex* a, fftw_complex* b);
fftw_complex* cmul(fftw_complex* a, fftw_complex* b);

// out = m1_rows * m2_cols
fftw_complex* matmul(fftw_complex* m1, size_t m1_cols, size_t m1_rows, fftw_complex* m2, size_t m2_cols, size_t m2_rows, fftw_complex* out);