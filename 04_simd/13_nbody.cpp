#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>
#include <iostream>
using namespace std;

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }
  for(int i=0; i<N; i++) {
    __m256 xi = _mm256_set1_ps(x[i]);
    __m256 yi = _mm256_set1_ps(y[i]);
    __m256 tempx,tempy;
    for(int j=0; j<N; j++) {
        __m256 xj = _mm256_load_ps(x);
        __m256 yj = _mm256_load_ps(y);
	__m256 rx = _mm256_sub_ps(xi,xj);
	__m256 ry = _mm256_sub_ps(yi,yj);
        __m256 vm = _mm256_load_ps(m);
	__m256 zeros = _mm256_setzero_ps();
	__m256 mask = _mm256_cmp_ps(xi,xj,_CMP_NEQ_OQ);
        __m256 r = _mm256_rsqrt_ps(_mm256_add_ps(_mm256_mul_ps(rx,rx), _mm256_mul_ps(ry, ry)));
	__m256 r3 = _mm256_mul_ps(r,_mm256_mul_ps(r,r));
	r3 = _mm256_blendv_ps(zeros,r3,mask);
	tempx = _mm256_mul_ps(rx,_mm256_mul_ps(vm,r3));
	tempy = _mm256_mul_ps(ry,_mm256_mul_ps(vm,r3));
    }
    tempx = _mm256_hadd_ps(tempx,_mm256_setzero_ps());
    tempx = _mm256_hadd_ps(tempx,_mm256_setzero_ps());
    fx[i] = -(((float*)&tempx)[0]+((float*)&tempx)[4]);
    tempy = _mm256_hadd_ps(tempy,tempy);
    tempy = _mm256_hadd_ps(tempy,tempy);
    fy[i] = -(((float*)&tempy)[0]+((float*)&tempy)[4]);
    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
