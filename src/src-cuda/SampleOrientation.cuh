// grid class that contains dx,dy,Nx,Ny

#ifndef SAMPLEORIENTATION_CUH
#define SAMPLEORIENTATION_CUH


__device__ inline void GetPermutation(const double *, double *);
__device__ inline void ComputeMap1(const double &a, const double& aprime, const double *xyz, double *map1);
__device__ inline void ComputeMap2(const double &beta,const double *xyz, double *map2);
__device__ inline void ComputeMap3(const double *xyz, double *map3);
__device__ inline void ComputeMeasure(const double &t,double & f);
__device__ void GenerateSamples(const int Nsample,unsigned int seedL, int subsq, curandState_t &s1, double *axisAngle);

#endif
