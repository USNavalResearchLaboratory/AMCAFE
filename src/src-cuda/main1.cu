
#include "Grid.cuh"
#include "VoxelsCA.cuh"
#include "BasePlate.cuh"
#include "TempField.cuh"
#include "iostream"
#include "vector"
#include <math.h>
#include <chrono>
#include <thread>
#include <algorithm>
#include "fstream"

static void HandleError(cudaError_t err) {
    if (err != cudaSuccess) {
        printf("%s \n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char *argv[])
{
  //-----------------------------------------------
  // initialize and create base plate
  // set up all pointers for arrays in class
  double *d_lcoor,*d_lcoor2;
  // voxels
  int *d_gID,*d_ineighID,*d_neighptr,*d_vState,nBlocks,nThreads,
  double *d_cTheta,*d_extents,*d_centroidOct;
  // tempfield
  double *d_Temp,*d_ispvec;
  // initialize class variables
  auto texec1 = std::chrono::high_resolution_clock::now();
  std::string filbaseOut,filout,filLogOut,filParamIn;
  filParamIn = argv[1];
  Grid g(filParamIn);
  Grid *d_g;
  HandleError(cudaMallocManaged((void**)&d_g,sizeof(Grid)));
  HandleError(cudaMemcpy(d_g, &g, sizeof(Grid), cudaMemcpyHostToDevice));
  HandleError(cudaMallocManaged((void**)&d_lcoor,2*g.NpT*sizeof(double)));
  HandleError(cudaMallocManaged((void**)&d_lcoor2,2*g.NpT*sizeof(double)));
  HandleError(cudaMemcpy(d_lcoor,(g.lcoor), 2*g.NpT*sizeof(double), cudaMemcpyHostToDevice));
  HandleError(cudaMemcpy(d_lcoor2,(g.lcoor2), 2*g.NpT*sizeof(double), cudaMemcpyHostToDevice));
  TempField TempF(g);
  TempF.InitializeAnalytic();
  TempField *d_TempF;
  HandleError(cudaMallocManaged((void**)&d_TempF,sizeof(Grid)));
  HandleError(cudaMemcpy(d_TempF, &TempF, sizeof(TempField), cudaMemcpyHostToDevice));
  HandleError(cudaMallocManaged((void**)&d_ispvec,g.NpT*sizeof(int)));
  HandleError(cudaMemcpy(d_ispvec,(TempF.ispvec), g.NpT*sizeof(int), cudaMemcpyHostToDevice));
  std::vector<double> bpSites;
  GenerateGrainSites(g,bpSites);
  double *d_bpSites;
  HandleError(cudaMallocManaged((void**)&d_bpSites,(int)(bpSites.size())*sizeof(double)));
  HandleError(cudaMemcpy(d_bpSites,bpSites, (int)bpSites.size()*sizeof(double), cudaMemcpyHostToDevice));
  VoxelsCA vox(g);
  vox.nGrain = bpSites.size()/3;
  VoxelsCA *d_vox;
  HandleError(cudaMallocManaged((void**)&d_vox,sizeof(VoxelsCA)));
  HandleError(cudaMemcpy(d_vox, &vox, sizeof(VoxelsCA), cudaMemcpyHostToDevice));
  nThreads=1024;
  HandleError(cudaMallocManaged((void**)&d_gID,Ntot*sizeof(int)));
  HandleError(cudaMallocManaged((void**)&d_vState,Ntot*sizeof(int)));
  HandleError(cudaMallocManaged((void**)&d_extents,Ntot*sizeof(double)));
  HandleError(cudaMallocManaged((void**)&d_centroidOct,3*Ntot*sizeof(double)));
  HandleError(cudaMemset(d_gID,Ntot*sizeof(int)));
  HandleError(cudaMemset(d_vState,Ntot*sizeof(int)));
  HandleError(cudaMemset(d_extents,Ntot*sizeof(double)));
  HandleError(cudaMemset(d_centroidOct,3*Ntot*sizeof(double)));
  HandleError(cudaMallocManaged((void**)&d_cTheta,4*vox.nGrain*sizeof(double)));
  int Ntot=g.nX[0]*g.nX[1]*g.nX[2];
  nBlocks=Ntot/nThreads;
  CreateBasePlateGrains<<<nBlocks,nThreads>>>(d_vox,d_gID,d_vState,d_g,d_bpSites,d_extends, &
					    d_centroidOct,Ntot);
  nBlocks=vox.nGrain/nThreads;
  CreateBasePlateOrientations<<<nBlocks,nThreads>>>(d_vox,d_cTheta);
  cudaFree(d_bpSites);
  // end initialize and create baseplate
  //-----------------------------------------------

  //-----------------------------------------------
  // run simulation loop 

  // end run simulation loop 
  //-----------------------------------------------

  //free cpu memory

  //free memory gpu
  cudaFree(d_gID);
  cudaFree(d_vState);
  cudaFree(d_extents);
  cudaFree(d_centroidOct);
  cudaFree(d_cTheta);
  
  return 0;
}
