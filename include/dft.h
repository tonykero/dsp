#pragma once

#include <fftw3.h>

fftw_complex* cadd(fftw_complex* a, fftw_complex* b) {
    (*a)[0] += (*b)[0];
    (*a)[1] += (*b)[1];
    return a;
}

fftw_complex* cmul(fftw_complex* a, fftw_complex* b) {
    double ac = (*a)[0] * (*b)[0];
    double bd = (*a)[1] * (*b)[1];
    double ad = (*a)[0] * (*b)[1];
    double bc = (*a)[1] * (*b)[0];

    (*a)[0] = ac - bd;
    (*a)[1] = ad + bc;

    return a;
}

// out = m1_rows * m2_cols
fftw_complex* matmul(fftw_complex* m1, size_t m1_cols, size_t m1_rows, fftw_complex* m2, size_t m2_cols, size_t m2_rows, fftw_complex* out) {
    if(m1_cols != m2_rows) { return nullptr; }

    for(size_t i = 0; i < m1_rows; i++) {
        for(size_t j = 0; j < m2_cols; j++) {
            fftw_complex acc = {0.0,0.0};
            for(size_t k = 0; k < m1_cols; k++) {
                fftw_complex a_tmp = {m1[i * m1_cols + k][0], m1[i * m1_cols + k][1]};
                cadd(&acc, cmul(&a_tmp,&m2[k * m2_cols + j]));
            }
            out[i * m2_cols + j][0] = acc[0];
            out[i * m2_cols + j][1] = acc[1];
        }
    }

    return out;
}

void naive_dft(fftw_complex* in, fftw_complex* out, size_t len) {
    fftw_complex* dft_mat = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * len * len);

    for(size_t i = 0; i < len; i++) {
        for(size_t j = 0; j < len; j++) {
            double theta = -2.0 * M_PI * i * j / len;
            dft_mat[i * len + j][0] = cos(theta);
            dft_mat[i * len + j][1] = sin(theta);
        }
    }

    matmul(&dft_mat[0], len,len,in,1,len,out);
}