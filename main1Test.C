/*
This runs the cellular automata model using a given temperature
input (obtained from Moose runs) to simulate microstructure
evolution of AM builds (20190123)

*/


// include files here
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

int main()
{
  /*-----------------------------------------------
    initialization step */
  MPI_Init(NULL,NULL);
  int nprocs,myid;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  
  int NtM,cc1,nDim;
  std::vector<int> nXM,nX,nTTemp(3,0);
  std::vector<double> dX,dXM,LX;
  double tL,mL,c0,Gamma,dP,dL,dtM,muN,layerThickness,T0;
  std::string filbaseTemp,filbaseOut,filout,neighOrder;
  double mu = 2e11; //2e11;//5e9 // rate for baseplate voronoi;
  double heightBase = .1e-3;
  // schwalbach parameters
  int patternID;
  double beamVel,beamSpacing,beamPower,wEst,cP,rho,kappa,beamEta,rcut;
  std::vector<double> beamSTD,LxAll;
  // schwalbach parameters

  nXM = {20,20,20};
  nX = {32,32,32};
  nDim = nX.size();
  LX = {.002,.002,.002};
  dX.assign(nDim,0.0);
  dXM.assign(nDim,0.0);
  for (int j=0;j<nDim;++j){
    dX[j] = LX[j]/double(nX[j]);
    dXM[j] = LX[j]/double(nXM[j]);
  }
  NtM = 50;
  dtM = .05; // must set based on moose results
  tL = 1609; // K
  mL = -10.9; // (K / wt%)
  dL = 3e-9; // (m^2/s)
  Gamma = 1e-7;  // (K m)
  muN = 9e-2; // rate for nucleation
  dP = .48;
  c0 = 4.85; // (wt %)
  filbaseTemp = "/Users/kteferra/Documents/research/projects/AMICME/codes/CA/tempData/tempField0.";
  neighOrder = "second"; // can equal "first", "second", "third"
  rho = 8000.0; // kg /m^3
  cP = 502.0; // J/kg-K)
  kappa = 18.0; // W/(m-K)
  layerThickness = 5e-5; //30e-6; // m
  Grid g(dX,nX,tL,mL,c0,Gamma,dP,dL,muN,rho,cP,kappa,layerThickness,neighOrder,nDim);
  Partition part(g,myid,nprocs);
  part.PartitionGraph();
  BasePlate bp(g,heightBase,mu, part);
  TempField TempF(g,part,bp);
  // initialize appropriate temperature model
  //TempF.InitializeMoose(filbaseTemp,NtM,dtM,nXM,dXM);
  beamVel = 70e-3; // m/s
  beamSTD = {5e-5,5e-5,8e-5}; //  {20e-6,20e-6,20e-6}; // m
  beamPower = 20; //300; // W
  beamEta = 1.0;
  wEst  = pow(8*beamPower/(exp(1.0)*M_PI*rho*cP*(tL-298.0)*beamVel),.5); // see EQ (1) in schwalbach
  beamSpacing = 4.0/3.0*beamSTD[1]; //wEst*.8;
  patternID = 0;
  LxAll = {.002+beamSpacing,.002+beamSpacing,.002+beamSpacing}; //{.01,.01,.01}; // domain of entire piece, used for laser pattern (m)
  nTTemp = {int(floor(LxAll[0]/beamSpacing))+1,int(floor(LxAll[1]/beamSpacing))+1,
	    int(floor(LxAll[2]/layerThickness))+1};
  T0 = 300.0; // initial temperature in (K)
  TempF.InitializeSchwalbach(patternID,beamSTD,beamSpacing,beamVel,beamPower,beamEta,LxAll,T0);
  //TempF.SchwalbachTempCurr();
  //TempF.SchwalbachDDtTemp();
  //TempF.ReadCSVMoose2();
  TempF.Test2();
  //TempF.ComputeDDtTemp();
    
  
  VoxelsCA vox(g,TempF, part);
  //vox.InitializeVoxels(bp);
  vox.InitializeTest1();

  /*
  if (part.myid==0){
    std::cout << bp.Ngrain << std::endl;
    for (int j=0;j<bp.Ngrain;++j){std::cout <<bp.gNucleus[j] << ",";}
  }
  MPI_Barrier(MPI_COMM_WORLD);
  */
  /*-----------------------------------------------
    execute simulation */
  cc1=0;
  int outskip=20,indOut;
  std::vector<int> filinds,out2(3,0),j123(3,0);
  std::vector<double> filtime;
  int icheck = 1,ichecktmp,cc2=0, irep=0;
  filbaseOut = "CA3D";
  out2 = {6,6,3}; // the increment to skip output per direction
  while (icheck!=0){
    cc2+=1;
    j123[2] = floor(TempF.tInd /(nTTemp[0]*nTTemp[1]));
    j123[1] = floor((TempF.tInd - (nTTemp[0]*nTTemp[1])*j123[2])/ nTTemp[0]);
    j123[0] = TempF.tInd - (nTTemp[0]*nTTemp[1])*j123[2] - nTTemp[0]*j123[1];
    indOut = j123[2] % out2[2] + j123[1] % out2[1] + j123[0] % out2[0];
    //if (cc2>1000){icheck=0;}
    icheck=!std::all_of(vox.vState.begin(),vox.vState.end(),[](int n){return n==3;});
    ichecktmp = icheck;
    MPI_Allreduce(&ichecktmp,&icheck,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

    // update next step for voxels (time is updated in vox.ComputeExtents() )
    vox.UpdateVoxels();

    if (cc2%2==0){
      filinds.push_back(cc2);
      filtime.push_back(g.time);
      filout = filbaseOut+std::to_string(cc2);
      vox.WriteToVTU1(filout);
      filout = filbaseOut+".csv";
      vox.WriteCSVData(filout);
      cc1+=1;
      if (cc1 % 20 || icheck==0){
	filout=filbaseOut;
	vox.WriteToPVD(filout,filinds,filtime);
      } // if (cc1
      MPI_Barrier(MPI_COMM_WORLD);
    }
    } // while
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
