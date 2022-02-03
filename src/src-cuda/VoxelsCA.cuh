#ifndef VOXELSCA_CUH
#define VOXELSCA_CUH
#include "SetPrecision.cuh"
#include "Grid.cuh"
#include "TempField.cuh"
#include "BasePlate.cuh"
#include "vector"
#include "string"
#include <random>

// #include <math.h>


class VoxelsCA
{
 public:
  VoxelsCA(Grid &);
  void WriteToHDF1(const std::string &filname, const Grid &, const real *);
  void AddLayer1Macro(VoxelsCA *d_vx,Grid &g,Grid *d_g,real **d_cthptr,
                    real *d_troids, int *d_gid, int *d_vst, int &nbuf2);
  void CleanLayerMacro(VoxelsCA *dvx,int *dgid,real **dcthetaptr, int &ntot);
  void ConvertSolid1Macro(Grid *dg,int *vstate,real *dextents,const int &iswitch,int nThreads,const int &ntot);
  void SetLiquid3Macro(Grid *dg,int *dgid, int *dvstate,real *dtempval, real *dexts,int &nThreads, int &nBlocks);
  void UpdateVoxelsMacro(Grid *dg, Grid &gg,VoxelsCA *dvox,int *dgid, int *dvstate,real *dtempval,real *dexts,
                            real *troids, real *dctheta,int &nThreads, int &nBlocks, int &ntot);
  inline void getNumPowderGrains(const Grid &g,int &numPG)
  {
    real rate = g.lrate* (g.dX[0]*1e6)*(g.dX[1]*1e6)*
      g.nX[0]*g.nX[1]* (g.layerT*1e6);
    std::poisson_distribution<int> ng(rate);
    numPG =0;
    while (numPG == 0){numPG = ng(genlayer);} // end while
  } // end getNumPowderGrains
  int *gID,*ineighID,*ineighptr,*vState;
  real *cTheta,*extents,*centroidOct;
  real vmax;
  int nGrain,seed0,seed1,NzhBP;
  std::default_random_engine genlayer;
}; // end class VoxelCA

__global__ void addLayer1Part1(Grid *g,VoxelsCA *vx,real *xs, real *troids,
                          int *gD, int *vs, int *itmp, int numPG,int ntot);
__global__ void addLayer1Part2(Grid *g,VoxelsCA *vx, real *cTheta,
                          int *gD, int *itmp, int numPG,int ntot);
__global__ void addLayer1Part3(const Grid *g,int *gD, int *vs);
__global__ void copyGlobal(real *x1,real *x0, int n);
__global__ void getSites(Grid *g,VoxelsCA *vx,real *xs, int numPG);
__global__ void cleanLayerPart1(VoxelsCA *dvx,int *dgid,int *gvolflg,int Ntot);
__global__ void cleanLayerPart2(VoxelsCA *dvx,int *gvolflg, int *itmp);
__global__ void cleanLayerPart3(VoxelsCA *dvx,int *dgid,int *gvolflg, int *itmp,
                                real *ctmp, real *ctheta,int Ntot);
__global__ void cleanLayerPart4(VoxelsCA *dvx, int *gvolflg);
__global__ void convertSolid1Part1(Grid *dg,int *vst,real *dexts, int ntot);
__global__ void convertSolid1Part2(Grid *dg,int *vst,int ntot);
__global__ void setLiquid3(Grid *dg,int *dgid,int *dvstate,real *dtempval, real *dexts);
__global__ void reduceGlobalArray(int *ig, int n,int isw);
__device__ void projectPointLine(real *A, real *x0, real *x1, real *xproj);
__device__  void loadS(real S[][3],int sInd[][3]);
__device__ void loadRotMat(float omega, float *ax, float rRot[][3]);
void resizeGlobalArray(real **y, int &n0, int &n1);
void resizeArray(real **y, int &n);
#endif
