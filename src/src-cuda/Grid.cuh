// grid class that contains dx,dy,Nx,Ny

#ifndef GRID_CUH
#define GRID_CUH
#include "SetPrecision.cuh"
#include "iostream"
#include <math.h>
#include <vector>
#include <curand_kernel.h>
class Grid
{
 public:
  // default constructor
  Grid(std::string & filIn);  
  // any functions added here
  void readInputFile(std::string &filInput);
  void UpdateTime(const real &velo);
  void UpdateTime2(const real &dtIn);
  void UpdateTimeTest1(const real &velo);
  __device__ void UpdateLaser(real *lasercoor,real *lasercoor2);
  __device__ void GetNeighbors(int &jvox, int *ineigh);
  void SkipTime(const real &DelT);
  void ComputeNeighborhood(int &j, std::string & nO,std::vector<int> & nn);
  void ComputeNeighborhoodFirst(int &j, int & ntype, std::vector<int> &nn);
  void inline ComputeNeighborhoodMooreFirst(int &j, std::vector<int> &nn)
  {
    nn.assign(0,0);
    int j3,j2,j1,jst;
    j3 = floor(j /( nX[0]*nX[1]) );
    j2 = floor( (j - nX[0]*nX[1]*j3)/nX[0] );
    j1 = j - nX[0]*nX[1]*j3 - nX[0]*j2;
    std::vector<int> itmp = {-1,0,1};
    for (int i3 =0;i3<3;++i3){
      if ( (j3+itmp[i3]<0) || (j3+itmp[i3]>=nX[2])){continue;}
      for (int i2 =0;i2<3;++i2){
	if ( (j2+itmp[i2]<0) || (j2+itmp[i2]>=nX[1])){continue;}
	for (int i1 =0;i1<3;++i1){
	  if ( (j1+itmp[i1]<0) || (j1+itmp[i1]>=nX[0])){continue;}
	  jst = nX[0]*nX[1]*(j3+itmp[i3])+nX[0]*(j2+itmp[i2])+j1+itmp[i1];
	  if (jst !=j){nn.push_back(jst);}
	} // for (int i1...
      } // for (int i2...
    } // for (int i3...
  } // end inline Compute...
  real time,dX[3],meltparam[4],beamSTD[3],lX[3],offset[3],gsize[2],gth,gth0,*lcoor,*lcoor2,gbox[4],
    Avel,nvel,bmDelT;
  int nDim,tInd,nnodePerCell,ictrl,nZlayer,patternID,outint,outNL,nTsd,Nsd,Ntd,NpT,
    inewscanflg,inewlayerflg,isp,indlayer,ilaserLoc,Nzhg,nX[3],ntype,nlayerTot;
  real bmV,bmP,bhatch; // beam velocity,power, hatch spacing, nucleation rate
  real tL,tS,T0; // liquidus, solidus, room temp (K)
  real rNmax; // nucleation density (m^{-3})
  real layerT; // thickness of each layer
  real mu; // rate for Voronoi tessellation for baseplate
  real bpH; // base plate height
  real lrate; // layer rate for Voronoi powder microstructure 
  bool bcheck=0;
  curandState_t s1;
}; // end class Grid

__global__ void UpdateLaserGlobal(Grid *gg,real *lasercoor,real *lasercoor2);
__global__ void UpdateTime2Global(Grid *dg, const real dt);

#endif
