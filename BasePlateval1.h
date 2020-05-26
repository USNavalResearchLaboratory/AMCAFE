// grid class that contains dx,dy,Nx,Ny

#ifndef BASEPLATEVAL1_H
#define BASEPLATEVAL1_H

#include "vector"
#include "Grid.h"
#include "Partition.h"
#include <random>
//#include <math.h>

class BasePlate
{
 public:
  // default constructor: 
  BasePlate(const Grid &g, const double & bwIn, const double & bhatch, const double &height, const double &rate, Partition &);

  // any functions added here
  void GenerateNgrain();
  void GenerateVoxel();
  double height;  // height_of_base_plate
  double mu; // rate of Poisson process for determining # of grains
  double bwidth,bInc;
  int Ngrain,Nzh;
  std::vector<double> cTheta,jVals;
  std::vector<int> X,gNucleus; // container for grain
 private:
  const Grid *_xyz;
  std::default_random_engine generator;
  unsigned seed; // random number generator seed
  Partition *_part;
}; // end class BasePlate

#endif
