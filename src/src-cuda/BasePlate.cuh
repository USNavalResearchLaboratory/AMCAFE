#ifndef BASEPLATE_CUH
#define BASEPLATE_CUH

#include "vector"
#include "VoxelsCA.cuh"
#include "Grid.cuh"
#include <random>
//#include <math.h>


class VoxelsCA;

__global__ void createBasePlateGrains(VoxelsCA* vx, int *gid, int *vs,Grid *g,double *site,double *ext, double *centroids,const int ntot);
__global__ void createBasePlateOrientations(VoxelsCA *vx, double *ctheta);

void GenerateGrainSites(const Grid &, std::vector<double> & );


#endif
