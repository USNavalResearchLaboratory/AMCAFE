

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



static void HandleError( cudaError_t err,
                         const char *file,
                         int line ) {
    if (err != cudaSuccess) {
        printf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}

#define HandleError( err ) (HandleError( err, __FILE__, __LINE__ ))


int main(int argc, char *argv[])
{
  //-----------------------------------------------
  // initialize and create base plate
  // set up all pointers for arrays in class
  double *d_lcoor,*d_lcoor2;
  // voxels
  int *d_gID,*d_ineighID,*d_neighptr,*d_vState,*d_itmp,nBlocks,
    nThreads,*d_ispvec;
  double *d_cTheta,*d_extents,*d_centroidOct;
  // tempfield
  double *d_Tempvals;
  // initialize class variables
  auto texec1 = std::chrono::high_resolution_clock::now();
  std::string filbaseOut,filout,filLogOut,filParamIn;
  filParamIn = argv[1];
  Grid g(filParamIn);
  Grid *d_g;
  int Ntot=g.nX[0]*g.nX[1]*g.nX[2];
  HandleError(cudaMallocManaged((void**)&d_g,sizeof(Grid)));
  HandleError(cudaMemcpy(d_g, &g, sizeof(Grid), cudaMemcpyHostToDevice));
  HandleError(cudaMallocManaged((void**)&d_lcoor,2*g.NpT*sizeof(double)));
  HandleError(cudaMallocManaged((void**)&d_lcoor2,2*g.NpT*sizeof(double)));
  HandleError(cudaMemcpy(d_lcoor,(g.lcoor), 2*g.NpT*sizeof(double), cudaMemcpyHostToDevice));
  HandleError(cudaMemcpy(d_lcoor2,(g.lcoor2), 2*g.NpT*sizeof(double), cudaMemcpyHostToDevice));
  TempField TempF(g);
  TempF.InitializeAnalytic(g);
  TempField *d_TempF;
  HandleError(cudaMallocManaged((void**)&d_TempF,sizeof(Grid)));
  HandleError(cudaMalloc((void**)&d_Tempvals,Ntot*sizeof(double)));
  HandleError(cudaMemcpy(d_TempF, &TempF, sizeof(TempField), cudaMemcpyHostToDevice));
  HandleError(cudaMallocManaged((void**)&d_ispvec,g.NpT*sizeof(int)));
  HandleError(cudaMemcpy(d_ispvec,(TempF.ispvec), g.NpT*sizeof(int), cudaMemcpyHostToDevice));
  std::vector<double> Sites;
  GenerateGrainSites(g,Sites);
  VoxelsCA vox(g);
  vox.nGrain = Sites.size()/3;
  double *d_Sites;
  HandleError(cudaMallocManaged((void**)&d_Sites,vox.nGrain*3*sizeof(double)));
  HandleError(cudaMemcpy(d_Sites,Sites.data(), vox.nGrain*3*sizeof(double), cudaMemcpyHostToDevice));
  vox.cTheta=(double*)malloc(vox.nGrain*4*sizeof(double));
  VoxelsCA *d_vox;
  HandleError(cudaMallocManaged((void**)&d_vox,sizeof(VoxelsCA)));
  HandleError(cudaMemcpy(d_vox, &vox, sizeof(VoxelsCA), cudaMemcpyHostToDevice));
  nThreads=1024;
  HandleError(cudaMallocManaged((void**)&d_gID,Ntot*sizeof(int)));
  HandleError(cudaMallocManaged((void**)&d_vState,Ntot*sizeof(int)));
  HandleError(cudaMallocManaged((void**)&d_extents,Ntot*sizeof(double)));
  HandleError(cudaMallocManaged((void**)&d_centroidOct,3*Ntot*sizeof(double)));
  HandleError(cudaMemset(d_gID,0,Ntot*sizeof(int)));
  HandleError(cudaMemset(d_vState,0,Ntot*sizeof(int)));
  HandleError(cudaMemset(d_extents,0,Ntot*sizeof(double)));
  HandleError(cudaMemset(d_centroidOct,0,3*Ntot*sizeof(double)));
  HandleError(cudaMallocManaged((void**)&d_cTheta,4*vox.nGrain*sizeof(double)));
  nBlocks=Ntot/nThreads;  
  createBasePlateGrains<<<nBlocks,nThreads>>>(d_vox,d_gID,d_vState,d_g,d_Sites,d_extents, 
					    d_centroidOct,Ntot);
  nThreads=128;
  nBlocks=vox.nGrain/nThreads;
  cudaDeviceSynchronize();
  createBasePlateOrientations<<<nBlocks,nThreads>>>(d_vox,d_cTheta,d_g);
  HandleError( cudaPeekAtLastError() );
  cudaDeviceSynchronize();
  cudaFree(d_Sites);
  Sites.clear();
  Sites.shrink_to_fit();
  // end initialize and create baseplate
  //-----------------------------------------------


  //-----------------------------------------------
  // run simulation loop 
  int indout,nlayerTot,npg,nbuf1,nbuf2;
  nlayerTot=int(ceil( (double)(g.nX[2]-g.Nzhg)/(double)g.nZlayer));

  // addlayer and update laser
  vox.AddLayer1Macro(d_vox,g,d_g,&d_cTheta,d_centroidOct,d_gID,d_vState,nbuf2);
  UpdateLaserGlobal<<<1,1>>>(d_g,d_lcoor,d_lcoor2);

  
  HandleError(cudaMemcpy(&(vox.nGrain), &(d_vox->nGrain), sizeof(int), cudaMemcpyDeviceToHost));
  HandleError(cudaMemcpy(vox.gID, d_gID, Ntot*sizeof(int), cudaMemcpyDeviceToHost));
  HandleError(cudaMemcpy(vox.cTheta, d_cTheta, nbuf2*sizeof(double), cudaMemcpyDeviceToHost));
  filout="tmp";
  vox.WriteToHDF1(filout, g, TempF.TempCurr);



  
  //while (!g.bcheck){
  TempF.tInd = int(round(g.time/TempF.DelT));
  calcTemptInd<<<1,1>>>(d_g,d_TempF);
  nThreads=512; nBlocks=Ntot/nThreads;
  TempF.AnalyticalTempCalcMacro(d_g,d_TempF, d_Tempvals,d_lcoor,d_lcoor2,d_ispvec,
                                        nThreads, nBlocks, Ntot);
  // call global temp analytic
  // call global updatevoxels; 
  UpdateLaserGlobal<<<1,1>>>(d_g,d_lcoor,d_lcoor2);
  HandleError(cudaMemcpy(&g.bcheck,&(d_g->bcheck),sizeof(int),cudaMemcpyDeviceToHost));
  HandleError(cudaMemcpy(&g.inewlayerflg,&(d_g->inewlayerflg),sizeof(int),cudaMemcpyDeviceToHost));


  if (g.inewlayerflg==1){ 
    HandleError(cudaMemcpy(&(vox.nGrain), &(d_vox->nGrain), sizeof(int), cudaMemcpyDeviceToHost));
    vox.CleanLayerMacro(g,d_vox,d_gID,&d_cTheta);
  }    
    // cpu: assign indoutput and if statement 
    // if true, then do memcpy of necessary arrays and write out hdf5


    // right now g.updatetime is here, think about putting it in update laser
  g.UpdateTime2(TempF.DelT);
  UpdateTime2Global<<<1,1>>>(d_g,TempF.DelT);
  if (g.inewlayerflg==1){
    vox.AddLayer1Macro(d_vox,g,d_g,&d_cTheta,d_centroidOct,d_gID,d_vState,nbuf2);
  }



  //}





  // end run simulation loop 
  //-----------------------------------------------

  //free cpu memory

  //free memory gpu

  cudaFree(d_lcoor);
  cudaFree(d_lcoor2);
  //cudaFree(d_neighID);
  //cudaFree(d_neighptr);
  cudaFree(d_Tempvals);
  cudaFree(d_ispvec);
  cudaFree(d_gID);
  cudaFree(d_vState);
  cudaFree(d_extents);
  cudaFree(d_centroidOct);
  cudaFree(d_cTheta);
  cudaFree(d_vox);
  cudaFree(d_g);
  cudaFree(d_TempF);
  
  
  return 0;
}
