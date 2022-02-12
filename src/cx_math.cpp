#include "cx_math.hpp"

fftw_complex* cadd(fftw_complex* a, fftw_complex* b) {
    (*a)[0] += (*b)[0];
    (*a)[1] += (*b)[1];
    return a;
}
fftw_complex* csub(fftw_complex* a, fftw_complex* b) {
    (*a)[0] -= (*b)[0];
    (*a)[1] -= (*b)[1];
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
    if(m1_cols != m2_rows) { return NULL; }

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