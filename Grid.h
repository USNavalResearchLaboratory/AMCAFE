// grid class that contains dx,dy,Nx,Ny

#ifndef GRID_H
#define GRID_H
#include "iostream"
#include <vector>
class Grid
{
 public:
  // default constructor
  Grid( const std::vector<double> & dxIn, const std::vector<int> & nXIn, 
	const double &tL,
	const double &mL,const double &c0,const double &Gamma, 
	const double &kP,const double &dL, const double & muN,
	const double & rhoIn, const double & cPIn, const double & kappaIn,
	const double &layerThicknessIn, const std::string neighOrderIn,const int &nDimIn);

  // any functions added here
  void UpdateTime(const double &velo);
  void UpdateTime2(const double &dtIn);
  void UpdateTimeTest1(const double &velo);
  void SkipTime(const double &DelT);
  void ComputeNeighborhood(int &j, std::string & nO,std::vector<int> & nn);
  void ComputeNeighborhoodFirst(int &j, std::vector<int> &nn);
  double deltaXmin,deltaTcheck,dt,time;
  std::vector<double> dX;
  std::vector<int> nX;
  int nDim,tInd,nnodePerCell;
  double tL; // liquid temp (K)
  double mL; // liquidus slope of (K/wt%)
  double c0; // initial concentration (wt %)
  double Gamma; // Gibbs-Thompson coefficient (K m)
  double kP; // partition coefficient 
  double dL; // diffusion coefficient (m^2/s)
  double muN; // nucleation rate coefficient
  double rho; // material density (kg/m^3)
  double kappa; // material thermal conductivity (W/(m-K))
  double cP; // specific heat capacity (J/(kg-K))
  double layerT; // thickness of each layer
  double ethresh; // 1- cut off probability for Delta t max (see latex notes)
  double deltaThresh; // threshold value of pmf for allowing if voxel can be captured
  std::string neighOrder; // order of neighborhood
}; // end class Grid

#endif
