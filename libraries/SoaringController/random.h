//
//  random.hpp
//  
//
//  Created by Iain Guilliard on 6/09/2016.
//
//

#ifndef random_h
#define random_h

#include <stdint.h>

#define XSADD_FLOAT_MUL (1.0f / 16777216.0f)
#define MAX_GAUSS_SAMPLES 1000
void gaussian(float* y1, float* y2);
uint32_t xorshift128(void);
void cholesky44(float (&A)[4][4], float (&L)[4][4]);
void multivariate_normal(float (&samp)[MAX_GAUSS_SAMPLES][4], float (&mean)[4], float (&cov)[4][4], const int size);

extern uint32_t s0[];
#endif /* random_h */
