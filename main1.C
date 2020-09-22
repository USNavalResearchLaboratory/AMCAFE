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

int main(int argc, char *argv[])
{
  /*-----------------------------------------------
    initialization step */
  auto texec1 = std::chrono::high_resolution_clock::now();
  MPI_Init(NULL,NULL);
  int nprocs,myid;
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  int cc1;
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
  if (g.ictrl==4){
    // this is a test case scenario
    TempF.InitializeSchwalbach(); 
    //TempF.Test2ComputeTemp(1.02*g.tL,.97*g.tL,514880,0.0);
    TempF.Test2ComputeTemp(0.97*g.tL,.97*g.tL,0.0,0.0);
  } else{
    //TempF.InitializeSchwalbach(g); 
    //TempF.SchwalbachTempCurr();
    TempF.InitializeAnalytic();
    TempF.AnalyticTempCurr(g.time,TempF.TempCurr,part.icellidLoc,Ntot);
  }
  VoxelsCA vox(g,TempF, part);
  if (g.ictrl==4){
  //vox.InitializeTest1();
  vox.InitializeTest2();
  } else{
    vox.InitializeVoxels(bp);
  }
  /*-----------------------------------------------
    execute simulation */
  cc1=0;
  int indOut,nTmax=TempF.nTTemp[0]*TempF.nTTemp[1]*TempF.nTTemp[2],
    nFils,iNL;
  double filSize=(g.nX[0]/1e2*g.nX[1]/1e2*g.nX[2]/1e2*4*(9+3+3+8)+12)/1e3;
  nFils = int(ceil(filSize/1.5 ));
  std::vector<int> filinds,out2(2,0),j123(3,0);
  std::vector<double> filtime;
  int icheck = 1,ichecktmp,cc2=0, irep=0;
  std::ofstream fplog;
  filbaseOut = "CA3D"+filParamIn.substr(0,filParamIn.find("."));
  filLogOut="CA3D"+filParamIn.substr(0,filParamIn.find("."))+".log";
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
    indOut = (TempF.nTTemp[0]*TempF.nTTemp[1]*j123[2]+TempF.nTTemp[0]*j123[1]+j123[0]) % g.outint;
    iNL=fmod(TempF.tInd,TempF.nTTemp[0]*TempF.nTTemp[1]);
    if (iNL==0){vox.CleanLayer();}
    if (irep==0){
      irep=1;
      if (indOut==0 || TempF.tInd ==nTmax || (iNL==0 && g.outNL==0)){
	filinds.push_back(TempF.tInd);
	filtime.push_back(g.time);
	filout = filbaseOut+std::to_string(TempF.tInd);
	cc1+=1;
        filout = filbaseOut+"_t"+std::to_string(TempF.tInd)+".csv";
	vox.WriteToHDF1(filout);
	/*
	vox.WriteToVTU1(filout);
	if (cc1 % 20 || TempF.tInd==(nTmax-1)){
	  filout=filbaseOut;
	  vox.WriteToPVD(filout,filinds,filtime);
	} // if (cc1
	*/
	MPI_Barrier(MPI_COMM_WORLD);
      } // (indOut==0 ...
    } // if (irep==0
    if (iNL==0){vox.AddLayer();}
    // update next step for voxels 
    vox.UpdateVoxels8();
    g.UpdateTime2(TempF.DelT);    
    // update temperature field
    if (TempF.tInd != int(round(g.time/TempF.DelT))){
      TempF.tInd = int(round(g.time/TempF.DelT));
      TempF.AnalyticTempCurr(g.time,TempF.TempCurr,part.icellidLoc,Ntot);
      irep=0;
    } // if (TempF.tInd !=
    auto texec2 = std::chrono::high_resolution_clock::now();
    auto delTexec = std::chrono::duration_cast<std::chrono::seconds>( texec2 - texec1 ).count();
    if (part.myid==0){std::cout << TempF.tInd<<","<<delTexec<< std::endl;}
    if (part.myid==0){fplog << TempF.tInd<<","<<delTexec<<std::endl;}
    } // while
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  return 0;
}
