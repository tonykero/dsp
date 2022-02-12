#pragma once

#include "cx_math.hpp"

void naive_dft(fftw_complex* in, fftw_complex* out, size_t len);

void naive_dft_rad2(fftw_complex* in, fftw_complex* out, size_t len);
void naive_dft_rad2_it(fftw_complex* in, fftw_complex* out, size_t len);