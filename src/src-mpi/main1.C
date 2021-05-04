#include "Grid.h"
#include "VoxelsCA.h"
#include "BasePlate.h"
#include "TempField.h"
#include "Partition.h"
#include "SampleOrientation.h"
#include "iostream"
#include "vector"
#include <math.h>
#include <chrono>
#include <thread>
#include <algorithm>
#include "mpi.h"
#include "fstream"

int main(int argc, char *argv[])
{
  /*-----------------------------------------------
    initialization step */
  auto texec1 = std::chrono::high_resolution_clock::now();
  MPI_Init(NULL,NULL);
  int nprocs,myid;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  std::string filbaseOut,filout,filLogOut,filParamIn;
  // schwalbach parameters
  double beamVel,beamPower,wEst,cP,rho,kappa,beamEta,rcut,T0targ;
  filParamIn = argv[1];
  Grid g(filParamIn,myid,nprocs);
  Partition part(g,myid,nprocs);
  //part.PartitionGraph();
  part.PartitionGraph2();
  BasePlate bp(g, part);
  TempField TempF(g,part,bp);
  int Ntot=part.ncellLoc+ part.nGhost;
  TempF.InitializeAnalytic();
  VoxelsCA vox(g,TempF, part);
  vox.InitializeVoxels(bp);
  /*-----------------------------------------------
    execute simulation */
  int indOut,nFils;
  double filSize=(g.nX[0]/1e2*g.nX[1]/1e2*g.nX[2]/1e2*4*(9+3+3+8)+12)/1e3;
  nFils = int(ceil(filSize/1.5 ));
  std::vector<int> filinds,out2(2,0);
  std::vector<double> filtime;
  int nlayerTot,icheck = 1,ichecktmp;
  std::ofstream fplog;
  bool bcheck=0;
  filbaseOut = "CA3D"+filParamIn.substr(0,filParamIn.find("."));
  filLogOut="CA3D"+filParamIn.substr(0,filParamIn.find("."))+".log";
  if (part.myid==0){
    fplog.open(filLogOut.c_str());
    fplog << "Time index= ,Total clock time passed(s)"<<std::endl;
  }
  nlayerTot=int(ceil( (double)(g.nX[2]-bp.Nzh)/(double)g.nZlayer));
  while (!bcheck){
    if (g.inewlayerflg==1){vox.AddLayer1();}
    // update temperature field
    TempF.tInd = int(round(g.time/TempF.DelT));
    g.UpdateLaser();
    TempF.AnalyticTempCurr(g.time,TempF.TempCurr,part.icellidLoc,Ntot);
    bcheck=g.indlayer>nlayerTot;

    // update next step for voxels 
    vox.UpdateVoxels();
    //write out
    indOut = TempF.tInd % g.outint;
    if (g.inewlayerflg==1){vox.CleanLayer();}
    if (indOut==0 || bcheck || (g.inewlayerflg==1 && g.outNL==0)){
      filinds.push_back(TempF.tInd);
      filtime.push_back(g.time);
      filout = filbaseOut+std::to_string(TempF.tInd);
      filout = filbaseOut+std::to_string(TempF.tInd);
      vox.WriteToHDF1(filout);
      MPI_Barrier(MPI_COMM_WORLD);
    } // (indOut==0 ...
    g.UpdateTime2(TempF.DelT);    
    auto texec2 = std::chrono::high_resolution_clock::now();
    auto delTexec = std::chrono::duration_cast<std::chrono::seconds>( texec2 - texec1 ).count();
    if (part.myid==0){std::cout << TempF.tInd<<","<<delTexec<< std::endl;}
    if (part.myid==0){fplog << TempF.tInd<<","<<delTexec<<std::endl;}
    } // while
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
