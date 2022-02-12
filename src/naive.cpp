#include "dft.hpp"

#define _USE_MATH_DEFINES
#include <math.h>

#include <fftw3.h>

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
    fftw_free(dft_mat);
}

void naive_dft_rad2_rec(fftw_complex* in, fftw_complex* out, size_t len, size_t stride) {
    if(len == 1) {
        out[0][0] = in[0][0];
        out[0][1] = in[0][1];
    } else {
        size_t half = len/2;
        naive_dft_rad2_rec(in,out,half,2*stride);
        naive_dft_rad2_rec(in+stride,out+half,half,2*stride);
        for(int k = 0; k < half; k++) {
            fftw_complex p1 = {out[k][0], out[k][1]};
            fftw_complex p2 = {out[k][0], out[k][1]};
            
            double theta = -2.0 * M_PI * k / len;
            fftw_complex q = {cos(theta), sin(theta)};
            cmul(&q, out+k+half);
            cadd(&p1, &q);
            csub(&p2, &q);

            out[k][0] = p1[0];
            out[k][1] = p1[1];
            out[k+half][0] = p2[0];
            out[k+half][1] = p2[1];
        }
    }
}

void naive_dft_rad2(fftw_complex* in, fftw_complex* out, size_t len) {
    naive_dft_rad2_rec(in,out,len,1);
}

//https://stackoverflow.com/a/9144870
//TODO: divide & conquer bit reversal
uint32_t reverse(uint32_t x,uint32_t bits)
{
    x = ((x >> 1) & 0x55555555u) | ((x & 0x55555555u) << 1);
    x = ((x >> 2) & 0x33333333u) | ((x & 0x33333333u) << 2);
    x = ((x >> 4) & 0x0f0f0f0fu) | ((x & 0x0f0f0f0fu) << 4);
    x = ((x >> 8) & 0x00ff00ffu) | ((x & 0x00ff00ffu) << 8);
    x = ((x >> 16) & 0xffffu) | ((x & 0xffffu) << 16);
    return x >>(sizeof(uint32_t)*8 - bits);
}

void naive_dft_rad2_it(fftw_complex* in, fftw_complex* out, size_t len) {
    int bits = (int)log2((int)len);
    for(int i = 0; i < len; i++) {
        uint32_t rev = reverse(i,bits);
        out[rev][0] = in[i][0];
        out[rev][1] = in[i][1];
    }

    for(int i = 0; i < bits; i++) {
        int s = 1 << (i+1);
        for(int j = 0; j < len; j += s) {
            for(int k = 0; k < s/2; k++) {
                int b = j+k;
                int half = s/2;
                fftw_complex p1 = {out[b][0], out[b][1]};
                fftw_complex p2 = {out[b][0], out[b][1]};

                double theta = -2.0 * M_PI * k / s;
                fftw_complex q = {cos(theta), sin(theta)};

                cmul(&q, out+b+half);
                cadd(&p1, &q);
                csub(&p2, &q);

                out[b][0] = p1[0];
                out[b][1] = p1[1];
                out[b+half][0] = p2[0];
                out[b+half][1] = p2[1];
            }
        }
    }
}