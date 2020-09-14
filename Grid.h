// grid class that contains dx,dy,Nx,Ny

#ifndef GRID_H
#define GRID_H
#include "iostream"
#include <math.h>
#include <vector>
class Grid
{
 public:
  // default constructor
  Grid(std::string & filIn, int & myid, int & nprocs);  
  // any functions added here
  void readInputFile();
  void UpdateTime(const double &velo);
  void UpdateTime2(const double &dtIn);
  void UpdateTimeTest1(const double &velo);
  void SkipTime(const double &DelT);
  void ComputeNeighborhood(int &j, std::string & nO,std::vector<int> & nn);
  void ComputeNeighborhoodFirst(int &j, std::string & ntype, std::vector<int> &nn);
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
  std::string filInput;
  double deltaXmin,deltaTcheck,dt,time;
  std::vector<double> dX,meltparam,beamSTD,LX,offset;
  std::vector<int> nX;
  int nDim,tInd,nnodePerCell,ictrl,nZlayer,patternID,outint,myid,nprocs,outNL,nTsd;
  double bmV,bmP,bhatch; // beam velocity,power, hatch spacing, nucleation rate
  double tL,tS,T0; // liquidus, solidus, room temp (K)
  double mL; // liquidus slope of (K/wt%)
  double c0; // initial concentration (wt %)
  double Gamma; // Gibbs-Thompson coefficient (K m)
  double kP; // partition coefficient 
  double dL; // diffusion coefficient (m^2/s)
  double muN; // nucleation rate coefficient
  double rho; // material density (kg/m^3)
  double kappa; // material thermal conductivity (W/(m-K))
  double dTempM; // mean undercooling for nucleation (K)
  double dTempS; // standard dev undercooling for nucleation (K)
  double rNmax; // nucleation density (m^{-3})
  double cP,dP; // specific heat capacity (J/(kg-K)), partition coeff
  double layerT; // thickness of each layer
  double ethresh; // 1- cut off probability for Delta t max (see latex notes)
  double deltaThresh; // threshold value of pmf for allowing if voxel can be captured
  double mu; // rate for Voronoi tessellation for baseplate
  double bpH; // base plate height
  double beamEta,T0targ,bmDelT;
  std::string ntype; // type of neighborhood: Moore or VonNeumann
  std::string neighOrder; // order of neighborhood
}; // end class Grid

#endif
