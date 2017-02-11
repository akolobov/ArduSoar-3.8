//
//  random.cpp
//  
//
//  Created by Iain Guilliard on 6/09/2016.
//
//

//#include <stdio.h>
#include <stdint.h>
#include "random.h"
#include <math.h>




// xorshift128
//
// See "Xorshift RNGs" by George Marsaglia, The Florida State University
//
// Code based on snippet from Wikipedia
// https://en.wikipedia.org/wiki/Xorshift
//
// This is prbably under the wikimedia commons licence
//
// We should use the MIT licenced library of Mutsuo Saito and
// Makoto Matsumoto:
// https://github.com/MersenneTwister-Lab/XSadd

// non zero seed
uint32_t s[] = { 12793,912503,84501,73115 };

uint32_t xorshift128(void) {
    uint32_t t = s[0];
    t ^= t << 11;
    t ^= t >> 8;
    s[0] = s[1]; s[1] = s[2]; s[2] = s[3];
    s[3] ^= s[3] >> 19;
    s[3] ^= t;
    return s[3];
}

// polar form of Box-Muller transform

void gaussian(float* y1, float* y2) {
    float x1, x2, w;
    uint32_t t;
    int i = 0;
    //printf("g0 s[] = [%u, %u, %u, %u]\n",s[0],s[1],s[2],s[3]);
    do {
        t = s[0];
        t ^= s[0] << 11;
        t ^= t >> 8;
        s[0] = s[1]; s[1] = s[2]; s[2] = s[3];
        s[3] ^= s[3] >> 19;
        s[3] ^= t;
        x1 = 2.0 * ((s[3] >> 8) * XSADD_FLOAT_MUL) - 1.0;
        t = s[0];
        t ^= s[0] << 11;
        t ^= t >> 8;
        s[0] = s[1]; s[1] = s[2]; s[2] = s[3];
        s[3] ^= s[3] >> 19;
        s[3] ^= t;
        x2 = 2.0 * ((s[3] >> 8) * XSADD_FLOAT_MUL) - 1.0;
        w = x1 * x1 + x2 * x2;
        i++;
    } while ( w >= 1.0 );
    //printf("g1 s[] = [%u, %u, %u, %u]\n",s[0],s[1],s[2],s[3]);
    //printf("i=%d\n",i);
    w = sqrtf( (-2.0 * logf( w ) ) / w );
    //printf("y1=%f y2=%f\n",x1 * w,x2 * w);
    *y1 = x1 * w;
    *y2 = x2 * w;
}

// cholesky decomposition of 4x4 matix A into L
// A must be positive definite (not tested for)
// L will be set to the lower triangular cholesky decomposition of A
void cholesky44(float (&A)[4][4], float (&L)[4][4]) {
    L[0][0] =  sqrtf(A[0][0]);
    L[1][0] =  1.0 / L[0][0] * (A[1][0]);
    L[1][1] =  sqrtf(A[1][1] - L[1][0] * L[1][0]);
    L[2][0] =  1.0 / L[0][0] * (A[2][0] );
    L[2][1] =  1.0 / L[1][1] * (A[2][1] - L[2][0] * L[1][0]);
    L[2][2] =  sqrtf(A[2][2] - (L[2][0] * L[2][0] + L[2][1] * L[2][1]));
    L[3][0] =  1.0 / L[0][0] * (A[3][0]);
    L[3][1] =  1.0 / L[1][1] * (A[3][1] - L[3][0] * L[1][0]);
    L[3][2] =  1.0 / L[2][2] * (A[3][2] - (L[3][0] * L[2][0] + L[3][1] * L[2][1]));
    L[3][3] =  sqrtf(A[3][3] - (L[3][0] * L[3][0] + L[3][1] * L[3][1] + L[3][2] * L[3][2]));
}

void multivariate_normal(float (&samp)[MAX_GAUSS_SAMPLES][4], float (&mean)[4], float (&cov)[4][4], const int size) {
    
    float m[4];
    float L[4][4];
    cholesky44(cov, L);
    
    for (int j = 0; j < size; j++) {
        gaussian(&m[0],&m[2]);
        gaussian(&m[1],&m[3]);
        samp[j][0] = mean[0] + L[0][0] * m[0];
        samp[j][1] = mean[1] + L[1][0] * m[0] + L[1][1] * m[1];
        samp[j][2] = mean[2] + L[2][0] * m[0] + L[2][1] * m[1] + L[2][2] * m[2];
        samp[j][3] = mean[3] + L[3][0] * m[0] + L[3][1] * m[1] + L[3][2] * m[2] + L[3][3] * m[3];
    }
}

