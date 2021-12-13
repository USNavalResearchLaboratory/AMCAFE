// member function definitions for Grid.C

#include "BasePlate.cuh"
#include "Grid.cuh"
#include "SampleOrientation.cuh"
#include "random"
#include <algorithm>
#include <numeric>
#include "iostream"
#include <math.h>

__global__ void createBasePlateGrains(VoxelsCA *vx, int *gD, int *vs,
                              Grid *gg,double *xSite,double *exts, double *troids,
                              const int &Ntot)
{
  int tid = threadIdx.x + blockIdx.x*blockDim.x,ng=vx->nGrain;
  int j1,j2,j3,igid=0;
  double xloc,yloc,zloc,dsq=1e6,dsqtmp;
  if (tid <Ntot){
    if (tid < (gg->Nzhg*gg->nX[0]*gg.nX[1])){
      j3 = floor(tid/(gg->nX[0]*gg->nX[1]));
      j2 = floor( (tid - gg->nX[0]*gg->nX[1]*j3)/gg->nX[0]);
      j1 = tid - gg->nX[0]*gg->nX[1]*j3 - gg->nX[0]*j2;
      zloc = (j3+.5)*gg->dX[2];
      yloc = (j2+.5)*gg->dX[1];
      xloc = (j1+.5)*gg->dX[0];
      for (int j=0;j<ng;++j){
        dsqtmp = pow(xloc - xSite[3*j-2],2.)+pow(yloc-xSite[3*j-1],2.)+pow(zloc-xSite[3*j],2.);
        if (dsqtmp==fmin(dsqtmp,dsq)){
          dsq=dsqtmp;
          igid=j;
        }
      } // for (int j=...
      gID[tid]=igid;
      if (tid < (gg->Nzhg-1)*gg->nX[0]*gg.nX[1]){
        vs[tid]=3;
      } else {
        vs[tid]=2;
      } // if (tid < gg
      troid[3*tid] = xloc;
      troid[3*tid+1] = yloc;
      troid[3*tid+2] = zloc;
    } // if (tid < (gg...
  } // if (tid < Ntot
}

__global__ void createBasePlateOrientations(VoxelsCA *vx, double *cTheta)
{
  int tid = threadIdx.x + blockIdx.x*blockDim.x,ng=vx->nGrain, &
    nsamp=1,subsq=0;
  curandState_t = s1;
  unsigned int seedL = vx->seed + tid*1000;
  double axAng[4];
  if (tid < ng){
    GenerateSamples(nsamp,seedL,subsq,s1, axAng);
    cTheta[4*tid]=axAng[0];
    cTheta[4*tid+1]=axAng[1];
    cTheta[4*tid+2]=axAng[2];
    cTheta[4*tid+3]=axAng[3];
  }
}

void GenerateGrainSites(const Grid &g, std::vector<double> & Xsite)
{
  double rate = g.mu* (g.dX[0]*1e6)*(g.dX[1]*1e6)* g.nX[0]*g.nX[1]* (height*1e6);
  std::poisson_distribution<int> ng(rate);
  int num_grain =0;
  while (num_grain == 0){
    num_grain = ng(generator);
  } // end while 
  Xsite.assign(3*num_grain,0.);
  std::uniform_real_distribution<double> xrand(0.0,1.0);
  double Lx = g.dX[0] * g.nX[0];
  double Ly = g.dX[1] * g.nX[1];
  double height= g.Nzhg*g.dX[2];
  for (int j=0; j<num_grain;++j){
    Xsite[3*j] = xrand(generator)*Lx;
    Xsite[3*j+1] = xrand(generator)*Ly;
    Xsite[3*j+2] = xrand(generator)*height;
  } // end for   
} // end GenerateGrainSites()

