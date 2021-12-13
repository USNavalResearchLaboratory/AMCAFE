// grid class that contains dx,dy,Nx,Ny

#ifndef BASEPLATE_CUH
#define BASEPLATE_CUH

#include "vector"
#include "Grid.cuh"
#include <random>
//#include <math.h>

void GenerateGrainSites(const Grid &, std::vector<double> & )
__global__ void createBasePlateGrains(VoxelsCA *, int *, int *,
                              Grid *,double *,double *, double *,
                              const int &)
__global__ void createBasePlateOrientations(VoxelsCA *, double *)
#endif
