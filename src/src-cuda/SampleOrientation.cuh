
#ifndef SAMPLEORIENTATION_CUH
#define SAMPLEORIENTATION_CUH

#include<curand_kernel.h>


__device__  void GetPermutation(const double *, int *);
__device__  void ComputeMap1(const double &a, const double& aprime, const double *xyz, double *map1);
__device__  void ComputeMap2(const double &beta,const double *xyz, double *map2);
__device__  void ComputeMap3(const double *xyz, double *map3);
__device__  void ComputeMeasure(const double &t,double & f);
__device__ void GenerateSamples(const int Nsample,unsigned int seedL, int subsq, curandState_t &s1, double *axisAngle);

#endif
