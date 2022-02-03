
#ifndef SAMPLEORIENTATION_CUH
#define SAMPLEORIENTATION_CUH
#include "SetPrecision.cuh"
#include<curand_kernel.h>


__device__  void GetPermutation(const real *, int *);
__device__  void ComputeMap1(const real &a, const real& aprime, const real *xyz, real *map1);
__device__  void ComputeMap2(const real &beta,const real *xyz, real *map2);
__device__  void ComputeMap3(const real *xyz, real *map3);
__device__  void ComputeMeasure(const real &t,real & f);
__device__ void GenerateSamples(const int Nsample,unsigned int seedL, int subsq, curandState_t &s1, real *axisAngle);

#endif
