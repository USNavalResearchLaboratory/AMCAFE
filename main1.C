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
  auto texec1 = std::chrono::high_resolution_clock::now();

  MPI_Init(NULL,NULL);
  int nprocs,myid,ictrl;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  int NtM,cc1,nDim;
  std::vector<int> nXM,nX;
  std::vector<double> dX,dXM,LX;
  double tL,mL,c0,Gamma,dP,dL,dtM,muN,layerThickness,T0,dTempM,dTempS,rNmax;
  std::string filbaseTemp,filbaseOut,filout,neighOrder,filLogOut;
  double mu; //2e11;//5e9 // rate for baseplate voronoi;
  double heightBase;
  // schwalbach parameters
  int patternID;
  double beamVel,beamPower,wEst,cP,rho,kappa,beamEta,rcut;
  std::vector<double> beamSTD;
  // schwalbach parameters
  ictrl=3;
  nXM = {20,20,20};
  nX = {128,128,64};
  nX = {64,64,32};
  LX = {.002,.002,.001};
  //nX = {128,32,16};
  //LX = {.001,.00025,12.5e-5};  
  nDim = nX.size();
  dX.assign(nDim,0.0);
  dXM.assign(nDim,0.0);
  for (int j=0;j<nDim;++j){
    dX[j] = LX[j]/double(nX[j]);
    dXM[j] = LX[j]/double(nXM[j]);
  }
  NtM = 50;
  dtM = .05; // must set based on moose results
  tL = 1609; // K
  dTempM = 2.5; // K (mean undercooling for nucleation)
  dTempS = 1.0; // K (standard dev undercooling for nucleation)
  rNmax = 7e14*pow(.001/dX[0],3.0); // m^{-3} maximum nucleation density
  mL = -10.9; // (K / wt%)
  dL = 3e-9; // (m^2/s)
  Gamma = 1e-7;  // (K m)
  muN = 9; // 9e-2; // rate for nucleation
  mu = 2e14; // 2e11
  heightBase = 2e-5;
  dP = .48;
  c0 = 4.85; // (wt %)
  filbaseTemp = "/Users/kteferra/Documents/research/projects/AMICME/codes/CA/tempData/tempField0.";
  neighOrder = "first"; // can only equal "first"
  rho = 8000.0; // kg /m^3
  cP = 502.0; // J/kg-K)
  kappa = 18.0; // W/(m-K)
  beamVel = 70e-3; //70e-3 // m/s
  beamSTD = {5e-5,10e-5,12.5e-5}; //  {20e-6,20e-6,20e-6}; // m
  beamPower = 70; //300; // W
  beamEta = 1.0;
  layerThickness = floor(beamSTD[2]/dX[2])*dX[2]; //30e-6; // m (layer thickness to be multiple of dX[2])
  Grid g(dX,nX,tL,mL,c0,Gamma,dP,dL,muN,rho,cP,kappa,layerThickness,neighOrder,dTempM,dTempS,rNmax,nDim);
  Partition part(g,myid,nprocs);
  part.PartitionGraph();
  BasePlate bp(g,heightBase,mu, part);
  TempField TempF(g,part,bp);
  // initialize appropriate temperature model
  //TempF.InitializeMoose(filbaseTemp,NtM,dtM,nXM,dXM);
  wEst  = pow(8*beamPower/(exp(1.0)*M_PI*rho*cP*(tL-298.0)*beamVel),.5); // see EQ (1) in schwalbach
  patternID = 0;
  //LxAll = {LX[0]+1000*beamSpacing,LX[1],LX[2]};
  T0 = 300.0; // initial temperature in (K)
  TempF.InitializeSchwalbach(patternID,beamSTD,beamVel,beamPower,beamEta,LX,T0);
  TempF.SchwalbachTempCurr();
  //TempF.SchwalbachDDtTemp();
  //TempF.ReadCSVMoose2();
  //TempF.Test2();
  //TempF.ComputeDDtTemp();
  VoxelsCA vox(g,TempF, part);
  vox.InitializeVoxels(bp);
  //vox.ExtentsInitialize();
  //vox.InitializeTest1();
  /*-----------------------------------------------
    execute simulation */
  cc1=0;
  int outskip=20,indOut,nTmax=TempF.nTTemp[0]*TempF.nTTemp[1]*TempF.nTTemp[2];
  std::vector<int> filinds,out2(3,0),j123(3,0);
  std::vector<double> filtime;
  int icheck = 1,ichecktmp,cc2=0, irep=0;
  std::ofstream fplog;
  filbaseOut = "CA3D";
  filLogOut="CA3D.log";
  out2 = {1,1,2}; // the increment to skip output per direction
  if (part.myid==0){fplog.open(filLogOut.c_str());}
  while (TempF.tInd<nTmax){
    cc2+=1;
    j123[2] = floor(TempF.tInd /(TempF.nTTemp[0]*TempF.nTTemp[1]));
    j123[1] = floor((TempF.tInd - (TempF.nTTemp[0]*TempF.nTTemp[1])*j123[2])/ TempF.nTTemp[0]);
    j123[0] = TempF.tInd - (TempF.nTTemp[0]*TempF.nTTemp[1])*j123[2] - TempF.nTTemp[0]*j123[1];
    indOut = j123[2] % out2[2] + j123[1] % out2[1] + j123[0] % out2[0];
    if (irep==0){
      irep=1;
      if (indOut==0 || TempF.tInd ==(nTmax-1)){ 
	filinds.push_back(TempF.tInd);
	filtime.push_back(g.time);
	filout = filbaseOut+std::to_string(TempF.tInd);
	vox.WriteToVTU1(filout);
	filout = filbaseOut+".csv";
	vox.WriteCSVData(filout);
	cc1+=1;
	if (cc1 % 20 || TempF.tInd==(nTmax-1)){
	  filout=filbaseOut;
	  vox.WriteToPVD(filout,filinds,filtime);
	} // if (cc1
	MPI_Barrier(MPI_COMM_WORLD);
      }
    }

    // update next step for voxels (time is updated in vox.ComputeExtents() )
    if (ictrl==0){
      vox.UpdateVoxels();
      vox.CheckTimeSkip();
    }
    if (ictrl==1){
      /*
      if (part.myid==0){std::cout << g.time << std::endl;}
      MPI_Barrier(MPI_COMM_WORLD);
      */
      vox.UpdateVoxels2();
      g.UpdateTime2(TempF.DelT);
    }
    if (ictrl==2){
      vox.UpdateVoxels3();
      g.UpdateTime2(TempF.DelT);
    }
    if (ictrl==3){
      vox.UpdateVoxels4();
      g.UpdateTime2(TempF.DelT);
      MPI_Barrier(MPI_COMM_WORLD);

    }
    // update temperature field
    if (TempF.tInd != int(round(g.time/TempF.DelT))){
      if (ictrl==0){vox.ExtentsInitialize();}
      TempF.tInd = int(round(g.time/TempF.DelT));
      TempF.SchwalbachTempCurr();
      //TempF.SchwalbachDDtTemp();
      irep=0;
    } // if (TempF.tInd !=
    /* This is for reading CSV file
    if (int(floor(g.time/TempF.dtM))!=TempF.indexM){
      TempF.ReadCSVMoose2();
      TempF.ComputeDDtTemp();
    }
    */
    auto texec2 = std::chrono::high_resolution_clock::now();
    auto delTexec = std::chrono::duration_cast<std::chrono::seconds>( texec2 - texec1 ).count();
    if (part.myid==0){std::cout << TempF.tInd<<","<< g.time/TempF.DelT<< std::endl;}
    if (part.myid==0){fplog << "Time index= "<<TempF.tInd<<",Total clock time passed(s)= "<<delTexec<<std::endl;}
    } // while
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
