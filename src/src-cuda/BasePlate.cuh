#ifndef BASEPLATE_CUH
#define BASEPLATE_CUH
#include "SetPrecision.cuh"
#include "vector"
#include "VoxelsCA.cuh"
#include "Grid.cuh"
#include <random>
//#include <math.h>


class VoxelsCA;

__global__ void createBasePlateGrains(VoxelsCA* vx, int *gid, int *vs,Grid *g,real *site,real *ext, real *centroids,const int ntot);
__global__ void createBasePlateOrientations(VoxelsCA *vx, real *ctheta, Grid *gg);

void GenerateGrainSites(const Grid &, std::vector<real> & );


#endif
