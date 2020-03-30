/*
This runs the cellular automata model using a given temperature
input (obtained from Moose runs) to simulate microstructure
evolution of AM builds (20190123)
*/
// include files here
#include "Grid.h"
#include "VoxelsCAval1.h"
#include "BasePlateval1.h"
#include "TempFieldval1.h"
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
  int nprocs,myid,ictrl;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  int NtM,cc1,nDim;
  std::vector<int> nXM,nX;
  std::vector<double> dX,dXM,LX;
  double tL,tS,mL,c0,Gamma,dP,dL,dtM,muN,layerThickness,T0,dTempM,dTempS,rNmax;
  std::string filbaseTemp,filbaseOut,filout,neighOrder,neighType,filLogOut;
  double mu,bwidth; //2e11;//5e9 // rate for baseplate voronoi;
  double heightBase;
  // schwalbach parameters
  int patternID;
  double beamVel,beamPower,wEst,cP,rho,kappa,beamEta,rcut,T0targ;
  std::vector<double> beamSTD;
  // schwalbach parameters
  bwidth = std::stod(argv[1]);
  rNmax = std::stod(argv[2]);
  ictrl=3;
  if (ictrl==4){
    nX = {32,32,32}; // KT: THIS IS FOR TEST
    LX = {.002,.002,.002}; // KT: THIS IS FOR TEST
  } else {
  nX = {128,128,64};
  LX = {nX[0]*1.875e-6,nX[1]*1.875e-6,nX[2]*1.875e-6};
  //nX = {128,32,16};
  //LX = {.001,.00025,12.5e-5};  
  }
  nXM = {20,20,20};
  nDim = nX.size();
  dX.assign(nDim,0.0);
  dXM.assign(nDim,0.0);
  for (int j=0;j<nDim;++j){
    dX[j] = LX[j]/double(nX[j]);
    dXM[j] = LX[j]/double(nXM[j]);
  }
  NtM = 50;
  dtM = .05; // must set based on moose results
  tL = 1733; // K
  tS = 1693; // K
  dTempM = (tL-tS)*.75; //7.5; // 2.5 // K (mean undercooling for nucleation)
  dTempS = (tL-tS)/3.0; //5.0; // 1.0 // K (standard dev undercooling for nucleation)
  mL = -10.9; // (K / wt%)
  dL = 3e-9; // (m^2/s)
  Gamma = 1e-7;  // (K m)
  muN = 9e-2; // 9 // 9e-2; 
  dP = .48;
  c0 = 4.85; // (wt %)
  neighOrder = "first"; // can only equal "first"
  neighType = "Moore"; // either "Moore" or "VonNeumann"
  rho = 8000.0; // kg /m^3
  cP = 502; // 502.0; // J/kg-K)
  kappa = 18; //18; //18.0; // W/(m-K)
  beamVel = 250e-3;//250e-3;//250e-3; //70e-3 // m/s
  /*
  layerThickness = 30e-6; // floor(beamSTD[2]/dX[2])*dX[2]; // (layer thickness to be multiple of dX[2])
  beamSTD = {5e-5,5e-5,layerThickness*1.5}; // m
  heightBase = dX[2];// layerThickness;//  ; 
  patternID = 1; // see TempField.C for description
  */
  layerThickness = 25e-6;
  beamSTD = {7.5e-5,7.5e-5,7.5e-5};
  heightBase = layerThickness;
  patternID = 1; // see TempField.C for description
  T0targ = 1500;//2000.0; // target peak temperature for one ellipsoid
  beamEta = 1.0;
  Grid g(dX,nX,tL,tS,mL,c0,Gamma,dP,dL,muN,rho,cP,kappa,layerThickness,neighOrder,dTempM,dTempS,rNmax,nDim,neighType,ictrl);
  Partition part(g,myid,nprocs);
  //part.PartitionGraph();
  part.PartitionGraph2();
  mu = 1e4/LX[0]/LX[1]/dX[2];// heightBase;//2e13; // 2e11 , 2e14  // rate for nucleation for baseplate 
  BasePlate bp(g,bwidth,beamSTD[1],heightBase,mu, part);
  TempField TempF(g,part,bp);
  wEst  = pow(8*beamPower/(exp(1.0)*M_PI*rho*cP*(tL-298.0)*beamVel),.5); // see EQ (1) in schwalbach
  T0 = 300.0; // initial temperature in (K)
  int Ntot=part.ncellLoc+ part.nGhost;
  if (ictrl==4){
    // this is a test case scenario
    TempF.InitializeSchwalbach(patternID,beamSTD,beamVel,T0targ,beamEta,LX,T0);
    //TempF.Test2ComputeTemp(1.02*tL,.97*tL,514880,0.0);
    TempF.Test2ComputeTemp(0.97*tL,.97*tL,0.0,0.0);
  } else{
    //TempF.InitializeSchwalbach(patternID,beamSTD,beamVel,T0targ,beamEta,LX,T0);
    //TempF.SchwalbachTempCurr();
    TempF.InitializeAnalytic(patternID,beamSTD,beamVel,LX,T0);
    TempF.AnalyticTempCurr(g.time,TempF.TempCurr,part.icellidLoc,Ntot);
  }
  VoxelsCA vox(g,TempF, part);
  if (ictrl==4){
  //vox.InitializeTest1();
  vox.InitializeTest2();
  } else{
    vox.InitializeVoxels(bp);
  }

  /*-----------------------------------------------
    execute simulation */
  cc1=0;
  int outskip=20,indOut,nTmax=TempF.nTTemp[0]*TempF.nTTemp[1]*TempF.nTTemp[2];
  std::vector<int> filinds,out2(2,0),j123(3,0);
  std::vector<double> filtime;
  int icheck = 1,ichecktmp,cc2=0, irep=0;
  std::ofstream fplog;
  filbaseOut = "CA3D";
  filLogOut="CA3D.log";
  out2 = {1,1}; // the increment to skip output per direction
  if (part.myid==0){
    fplog.open(filLogOut.c_str());
    fplog << "Time index= ,Total clock time passed(s)"<<std::endl;
  }
  while (TempF.tInd<=nTmax){
    icheck=!std::all_of(vox.vState.begin(),vox.vState.end(),[](int n){return n==3;});
    ichecktmp = icheck;
    MPI_Allreduce(&ichecktmp,&icheck,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    cc2+=1;
    j123[2] = floor(TempF.tInd /(TempF.nTTemp[0]*TempF.nTTemp[1]));
    j123[1] = floor((TempF.tInd - (TempF.nTTemp[0]*TempF.nTTemp[1])*j123[2])/ TempF.nTTemp[0]);
    j123[0] = TempF.tInd - (TempF.nTTemp[0]*TempF.nTTemp[1])*j123[2] - TempF.nTTemp[0]*j123[1];
    //indOut = j123[2] % out2[1] + (TempF.nTTemp[0]*j123[1]+j123[0]) % out2[0];
    indOut = (TempF.nTTemp[0]*TempF.nTTemp[1]*j123[2]+TempF.nTTemp[0]*j123[1]+j123[0]) % out2[0];
    if (irep==0){
      irep=1;
      if (indOut==0 || TempF.tInd ==(nTmax-1)){ 
	filinds.push_back(TempF.tInd);
	filtime.push_back(g.time);
	filout = filbaseOut+std::to_string(TempF.tInd);
	vox.WriteToVTU1(filout);
	//vox.WriteCSVData(filout);
	cc1+=1;
	if (cc1 % 20 || TempF.tInd==(nTmax-1)){
	  filout=filbaseOut;
	  vox.WriteToPVD(filout,filinds,filtime);
	} // if (cc1
	MPI_Barrier(MPI_COMM_WORLD);
      }
    }
    // update next step for voxels (time is updated in vox.ComputeExtents() )
    if (ictrl==3){
      if (fmod(TempF.tInd,TempF.nTTemp[0]*TempF.nTTemp[1])==0){
	filout = filbaseOut+"_t"+std::to_string(TempF.tInd)+".csv";
	vox.UpdateLayer(filout); // WriteCSVData1 called in UpdateLayer
      }
      vox.UpdateVoxels8();
      g.UpdateTime2(TempF.DelT);
    }
    if (ictrl==4){
      //std::cout << TempF.tInd<<",00,"<<g.time<<","<<g.time/TempF.DelT<<std::endl;
      vox.UpdateVoxels4();
      g.UpdateTime2(TempF.DelT);
      //TempF.Test2ComputeTemp(1.02*tL,.97*tL,514880,g.time);
      TempF.Test2ComputeTemp(0.97*tL,.97*tL,0.0,g.time);
    }
    // update temperature field
    if (TempF.tInd != int(round(g.time/TempF.DelT))){
      TempF.tInd = int(round(g.time/TempF.DelT));
      if (ictrl!=4){
	//TempF.SchwalbachTempCurr();
	TempF.AnalyticTempCurr(g.time,TempF.TempCurr,part.icellidLoc,Ntot);
      }
      irep=0;
    } // if (TempF.tInd !=
    auto texec2 = std::chrono::high_resolution_clock::now();
    auto delTexec = std::chrono::duration_cast<std::chrono::seconds>( texec2 - texec1 ).count();
    if (part.myid==0){fplog << TempF.tInd<<","<<delTexec<<std::endl;}
    if (part.myid==0){std::cout << TempF.tInd<<","<<delTexec<<std::endl;}
    } // while
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
