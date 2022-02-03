
#include "SetPrecision.cuh"
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
  auto texec1 = std::chrono::high_resolution_clock::now();
  real *d_lcoor,*d_lcoor2;
  // voxels
  int *d_gID,*d_vState,nBlocks,
    nThreads,*d_ispvec;
  real *d_cTheta,*d_extents,*d_centroidOct;
  // tempfield
  real *d_Tempvals;
  // initialize class variables
  std::string filbaseOut,filout,filLogOut,filParamIn;
  filParamIn = argv[1];
  filbaseOut = "CA3D"+filParamIn.substr(0,filParamIn.find("."));
  filLogOut="CA3D"+filParamIn.substr(0,filParamIn.find("."))+".log";
  Grid g(filParamIn);
  Grid *d_g;
  int Ntot=g.nX[0]*g.nX[1]*g.nX[2];
  HandleError(cudaMallocManaged((void**)&d_g,sizeof(Grid)));
  HandleError(cudaMemcpy(d_g, &g, sizeof(Grid), cudaMemcpyHostToDevice));
  HandleError(cudaMallocManaged((void**)&d_lcoor,2*g.NpT*sizeof(real)));
  HandleError(cudaMallocManaged((void**)&d_lcoor2,2*g.NpT*sizeof(real)));
  HandleError(cudaMemcpy(d_lcoor,(g.lcoor), 2*g.NpT*sizeof(real), cudaMemcpyHostToDevice));
  HandleError(cudaMemcpy(d_lcoor2,(g.lcoor2), 2*g.NpT*sizeof(real), cudaMemcpyHostToDevice));
  TempField TempF(g);
  TempF.InitializeAnalytic(g);
  TempField *d_TempF;
  HandleError(cudaMallocManaged((void**)&d_TempF,sizeof(Grid)));
  HandleError(cudaMalloc((void**)&d_Tempvals,Ntot*sizeof(real)));
  HandleError(cudaMemcpy(d_TempF, &TempF, sizeof(TempField), cudaMemcpyHostToDevice));
  HandleError(cudaMallocManaged((void**)&d_ispvec,g.NpT*sizeof(int)));
  HandleError(cudaMemcpy(d_ispvec,(TempF.ispvec), g.NpT*sizeof(int), cudaMemcpyHostToDevice));
  std::vector<real> Sites;
  GenerateGrainSites(g,Sites);
  VoxelsCA vox(g);
  vox.nGrain = Sites.size()/3;
  real *d_Sites;
  HandleError(cudaMallocManaged((void**)&d_Sites,vox.nGrain*3*sizeof(real)));
  HandleError(cudaMemcpy(d_Sites,Sites.data(), vox.nGrain*3*sizeof(real), cudaMemcpyHostToDevice));
  vox.cTheta=(real*)malloc(vox.nGrain*4*sizeof(real));
  VoxelsCA *d_vox;
  HandleError(cudaMallocManaged((void**)&d_vox,sizeof(VoxelsCA)));
  HandleError(cudaMemcpy(d_vox, &vox, sizeof(VoxelsCA), cudaMemcpyHostToDevice));
  nThreads=512;
  HandleError(cudaMallocManaged((void**)&d_gID,Ntot*sizeof(int)));
  HandleError(cudaMallocManaged((void**)&d_vState,Ntot*sizeof(int)));
  HandleError(cudaMallocManaged((void**)&d_extents,Ntot*sizeof(real)));
  HandleError(cudaMallocManaged((void**)&d_centroidOct,3*Ntot*sizeof(real)));
  HandleError(cudaMemset(d_gID,0,Ntot*sizeof(int)));
  HandleError(cudaMemset(d_vState,0,Ntot*sizeof(int)));
  HandleError(cudaMemset(d_extents,0,Ntot*sizeof(real)));
  HandleError(cudaMemset(d_centroidOct,0,3*Ntot*sizeof(real)));
  HandleError(cudaMallocManaged((void**)&d_cTheta,4*vox.nGrain*sizeof(real)));
  nBlocks=Ntot/nThreads;  
  createBasePlateGrains<<<nBlocks,nThreads>>>(d_vox,d_gID,d_vState,d_g,d_Sites,d_extents, 
					    d_centroidOct,Ntot);

  HandleError( cudaPeekAtLastError() );
  nThreads=128;
  nBlocks=vox.nGrain/nThreads;
  createBasePlateOrientations<<<nBlocks,nThreads>>>(d_vox,d_cTheta,d_g);
  HandleError( cudaPeekAtLastError() );
  HandleError(cudaFree(d_Sites));
  Sites.clear();
  Sites.shrink_to_fit();
  // end initialize and create baseplate
  //-----------------------------------------------


  //-----------------------------------------------
  // run simulation loop 
  int indOut,nbuf;
  std::ofstream fplog;
  fplog.open(filLogOut.c_str());
  fplog << "Time index= ,Total clock time passed(s)"<<std::endl;
  // addlayer and update laser
  vox.AddLayer1Macro(d_vox,g,d_g,&d_cTheta,d_centroidOct,d_gID,d_vState,nbuf);
  UpdateLaserGlobal<<<1,1>>>(d_g,d_lcoor,d_lcoor2);
  int cc=0;
  while (!g.bcheck){
    TempF.tInd = int(round(g.time/TempF.DelT));
    calcTemptInd<<<1,1>>>(d_g,d_TempF);
    nThreads=512; nBlocks=Ntot/nThreads;
    // call global temp analytic
    TempF.AnalyticalTempCalcMacro(d_g,d_TempF, d_Tempvals,d_lcoor,d_lcoor2,d_ispvec,
					  nThreads, nBlocks, Ntot);
    // call global updatevoxels;
    vox.UpdateVoxelsMacro(d_g,g,d_vox,d_gID,d_vState,d_Tempvals,d_extents,
    			d_centroidOct,d_cTheta,nThreads,nBlocks,Ntot);
    UpdateLaserGlobal<<<1,1>>>(d_g,d_lcoor,d_lcoor2);
    //HandleError( cudaPeekAtLastError() );
    HandleError(cudaMemcpy(&g.bcheck,&(d_g->bcheck),sizeof(bool),cudaMemcpyDeviceToHost));
    HandleError(cudaMemcpy(&g.inewlayerflg,&(d_g->inewlayerflg),sizeof(int),cudaMemcpyDeviceToHost));
    if (g.inewlayerflg==1){ 
      HandleError(cudaMemcpy(&(vox.nGrain), &(d_vox->nGrain), sizeof(int), cudaMemcpyDeviceToHost));
      cudaDeviceSynchronize();
      vox.CleanLayerMacro(d_vox,d_gID,&d_cTheta,Ntot);
    }    
    indOut = TempF.tInd % g.outint;
    if (indOut==0 || g.bcheck || (g.inewlayerflg==1 && g.outNL==0)){
      // if true, then do memcpy of necessary arrays and write out hdf5
      // (note nGrain is copied over from device in CleanLayerMacro
      nbuf=4*vox.nGrain;
      resizeArray(&(vox.cTheta),nbuf);
      HandleError(cudaMemcpy(vox.gID, d_gID, Ntot*sizeof(int), cudaMemcpyDeviceToHost));
      HandleError(cudaMemcpy(vox.vState, d_vState, Ntot*sizeof(int), cudaMemcpyDeviceToHost));
      HandleError(cudaMemcpy(TempF.TempCurr, d_Tempvals, Ntot*sizeof(real), cudaMemcpyDeviceToHost));
      HandleError(cudaMemcpy(vox.cTheta, d_cTheta, 4*vox.nGrain*sizeof(real), cudaMemcpyDeviceToHost));
      filout = filbaseOut+std::to_string(TempF.tInd);
      vox.WriteToHDF1(filout, g, TempF.TempCurr);
    } // if (indOut==0 ...
    g.UpdateTime2(TempF.DelT);
    UpdateTime2Global<<<1,1>>>(d_g,TempF.DelT);
    HandleError( cudaPeekAtLastError() ); 
    if (g.inewlayerflg==1){
      vox.AddLayer1Macro(d_vox,g,d_g,&d_cTheta,d_centroidOct,d_gID,d_vState,nbuf);
    }
    auto texec2 = std::chrono::high_resolution_clock::now();
    auto delTexec = std::chrono::duration_cast<std::chrono::seconds>( texec2 - texec1 ).count();
    fplog << TempF.tInd<<","<<delTexec<<std::endl;
    std::cout<<"time step and total time (s)"<<TempF.tInd<<","<<delTexec<<std::endl;
    cc+=1;
  } // while (!g.bcheck...
  cudaDeviceSynchronize();
  fplog.close();



  // end run simulation loop 
  //-----------------------------------------------

  //free cpu memory

  //free memory gpu

  HandleError(cudaFree(d_lcoor));
  HandleError(cudaFree(d_lcoor2));
  HandleError(cudaFree(d_Tempvals));
  HandleError(cudaFree(d_ispvec));
  HandleError(cudaFree(d_gID));
  HandleError(cudaFree(d_vState));
  HandleError(cudaFree(d_extents));
  HandleError(cudaFree(d_centroidOct));
  HandleError(cudaFree(d_cTheta));
  HandleError(cudaFree(d_vox));
  HandleError(cudaFree(d_g));
  HandleError(cudaFree(d_TempF));
  
  
  return 0;
}
