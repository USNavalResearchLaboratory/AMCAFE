// member function definitions for Grid.C
#include "SetPrecision.cuh"
#include "BasePlate.cuh"
#include "Grid.cuh"
#include "SampleOrientation.cuh"
#include "random"
#include <algorithm>
#include <numeric>
#include "iostream"
#include <math.h>
#include <curand_kernel.h>

__global__ void createBasePlateGrains(VoxelsCA *vx, int *gD, int *vs,
                              Grid *gg,real *xSite,real *exts, real *troids,
                              const int Ntot)
{
  int tid = threadIdx.x + blockIdx.x*blockDim.x,ng=vx->nGrain;
  int j1,j2,j3,igid=0,js,stride=blockDim.x*gridDim.x;
  real xloc,yloc,zloc,dsq=1e6,dsqtmp;
  js=tid;
  while (js <Ntot){

    if (js < (gg->Nzhg*gg->nX[0]*gg->nX[1])){
      j3 = js/ (gg->nX[0]*gg->nX[1]); // note that int/int is equivalent to floor
      j2 = (js - gg->nX[0]*gg->nX[1]*j3)/(gg->nX[0]);
      j1 = js - gg->nX[0]*gg->nX[1]*j3 - gg->nX[0]*j2;
      zloc = (j3+.5)*gg->dX[2];
      yloc = (j2+.5)*gg->dX[1];
      xloc = (j1+.5)*gg->dX[0];
      for (int j=0;j<ng;++j){
        dsqtmp = pow(xloc - xSite[3*j],2.)+pow(yloc-xSite[3*j+1],2.)+pow(zloc-xSite[3*j+2],2.);
        if (dsqtmp<dsq){
          dsq=dsqtmp;
          igid=j;
        }
      } // for (int j=...
      gD[js]=igid;
      if (tid < (gg->Nzhg-1)*gg->nX[0]*gg->nX[1]){
        vs[js]=3;
      } else {
        vs[js]=2;
      } // if (tid < gg
      troids[3*js] = xloc;
      troids[3*js+1] = yloc;
      troids[3*js+2] = zloc;
    } // if (js < (gg...
    js+=stride;
  } // whil (js < Ntot
}

__global__ void createBasePlateOrientations(VoxelsCA *vx, real *cTheta,Grid *gg)
{
  int tid = threadIdx.x + blockIdx.x*blockDim.x,ng=vx->nGrain, 
    nsamp=1,subsq=0,jc,stride=blockDim.x*gridDim.x;

  unsigned int seedL = vx->seed0 + tid*1000;
  real axAng[4];
  jc=tid;
  while (jc<ng){
    GenerateSamples(nsamp,seedL,subsq,gg->s1, axAng);
    cTheta[4*jc]=axAng[0];
    cTheta[4*jc+1]=axAng[1];
    cTheta[4*jc+2]=axAng[2];
    cTheta[4*jc+3]=axAng[3];
    jc+=stride;
  }
}

void GenerateGrainSites(const Grid &g, std::vector<real> & Xsite)
{

  real Lx = g.lX[0], Ly = g.lX[1],height = g.Nzhg*g.dX[2];
  real rate = g.mu*(Lx*1e6)*(Ly*1e6)*(height*1e6);
  std::default_random_engine generator;
  std::poisson_distribution<int> ng(rate);
  int num_grain =0;
  while (num_grain == 0){
    num_grain = ng(generator);
  } // end while 
  Xsite.assign(3*num_grain,0.);
  std::uniform_real_distribution<real> xrand(0.0,1.0);
  for (int j=0; j<num_grain;++j){
    Xsite[3*j] = xrand(generator)*Lx;
    Xsite[3*j+1] = xrand(generator)*Ly;
    Xsite[3*j+2] = xrand(generator)*height;
  } // end for   
} // end GenerateGrainSites()

