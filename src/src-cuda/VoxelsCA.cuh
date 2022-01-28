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
  void WriteToHDF1(const std::string &filname, const Grid &, const double *);
  void AddLayer1Macro(VoxelsCA *d_vx,Grid &g,Grid *d_g,double **d_cthptr,
                    double *d_troids, int *d_gid, int *d_vst, int &nbuf2);
  void CleanLayerMacro(VoxelsCA *dvx,int *dgid,double **dcthetaptr, int &ntot);
  void ConvertSolid1Macro(Grid *dg,int *vstate,double *dextents,const int &iswitch,int nThreads,const int &ntot);
  void SetLiquid3Macro(Grid *dg,int *dgid, int *dvstate,double *dtempval, double *dexts,int &nThreads, int &nBlocks);
  void UpdateVoxelsMacro(Grid *dg, Grid &gg,VoxelsCA *dvox,int *dgid, int *dvstate,double *dtempval,double *dexts,
                            double *troids, double *dctheta,int &nThreads, int &nBlocks, int &ntot);
  inline void getNumPowderGrains(const Grid &g,int &numPG)
  {
    double rate = g.lrate* (g.dX[0]*1e6)*(g.dX[1]*1e6)*
      g.nX[0]*g.nX[1]* (g.layerT*1e6);
    std::poisson_distribution<int> ng(rate);
    numPG =0;
    while (numPG == 0){numPG = ng(genlayer);} // end while
  } // end getNumPowderGrains
  int *gID,*ineighID,*ineighptr,*vState;
  double *cTheta,*extents,*centroidOct;
  double vmax;
  int nGrain,seed0,seed1,NzhBP;
  std::default_random_engine genlayer;
}; // end class VoxelCA

__global__ void addLayer1Part1(Grid *g,VoxelsCA *vx,double *xs, double *troids,
                          int *gD, int *vs, int *itmp, int numPG,int ntot);
__global__ void addLayer1Part2(Grid *g,VoxelsCA *vx, double *cTheta,
                          int *gD, int *itmp, int numPG,int ntot);
__global__ void addLayer1Part3(const Grid *g,int *gD, int *vs);
__global__ void copyGlobal(double *x1,double *x0, int n);
__global__ void getSites(Grid *g,VoxelsCA *vx,double *xs, int numPG);
__global__ void cleanLayerPart1(VoxelsCA *dvx,int *dgid,int *gvolflg,int Ntot);
__global__ void cleanLayerPart2(VoxelsCA *dvx,int *gvolflg, int *itmp);
__global__ void cleanLayerPart3(VoxelsCA *dvx,int *dgid,int *gvolflg, int *itmp,
                                double *ctmp, double *ctheta,int Ntot);
__global__ void cleanLayerPart4(VoxelsCA *dvx, int *gvolflg);
__global__ void convertSolid1Part1(Grid *dg,int *vst,double *dexts, int ntot);
__global__ void convertSolid1Part2(Grid *dg,int *vst,int ntot);
__global__ void setLiquid3(Grid *dg,int *dgid,int *dvstate,double *dtempval, double *dexts);
__global__ void reduceGlobalArray(int *ig, int n,int isw);
__device__ void projectPointLine(double *A, double *x0, double *x1, double *xproj);
__device__  void loadS(double S[][3],int sInd[][3]);
__device__ void loadRotMat(float omega, float *ax, float rRot[][3]);
void resizeGlobalArray(double **y, int &n0, int &n1);
void resizeArray(double **y, int &n);
#endif
