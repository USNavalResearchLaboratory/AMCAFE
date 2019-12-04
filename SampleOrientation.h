// grid class that contains dx,dy,Nx,Ny

#ifndef SAMPLEORIENTATION_H
#define SAMPLEORIENTATION_H

#include "vector"
#include <random>
//#include <math.h>

class SampleOrientation
{
 public:
  // default constructor: 
  SampleOrientation();

  // any functions added here
  void GetPermutation(std::vector<double> & xyz, std::vector<int> & P);
  void ComputeMap1(std::vector<double> & xyz,std::vector<double> &map1);
  void ComputeMap2(std::vector<double> & xyz,std::vector<double> &map2);
  void ComputeMap3(std::vector<double> & xyz,std::vector<double> &map3);
  void ComputeMeasure(double & t, double &f);
  void GenerateSamples(const int Nsample, unsigned int seedin, std::vector<double> &axisAngle);
  int Ngrid, NTgrid;
  double r,a,aprime,R1,beta;
  std::vector<double> coeffA;

}; // end class BasePlate

#endif
