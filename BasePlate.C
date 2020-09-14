// member function definitions for Grid.C

#include "BasePlate.h"
#include "Grid.h"
#include "Partition.h"
#include "SampleOrientation.h"
#include "random"
#include <algorithm>
#include <numeric>
#include "iostream"
#include <math.h>


#include <mpi.h>


// constructor
BasePlate::BasePlate(const Grid & g, Partition & part)
{
  _xyz = &g;
  _part = &part;
  mu = _xyz->mu;
  height = ceil(_xyz->bpH/_xyz->dX[2])*_xyz->dX[2];
  Nzh = std::min(_xyz->nX[2],int(ceil(_xyz->bpH/_xyz->dX[2])));
  height = Nzh*(_xyz->dX[2]);
  // if you want seed to be different every time program runs
  //seed = std::chrono::system_clock::now().time_since_epoch().count();
  seed = 1234567;
  std::default_random_engine generator(seed); 
} // end constructor
void BasePlate::GenerateNgrain()
{
  double rate = mu* (_xyz->dX[0]*1e6)*(_xyz->dX[1]*1e6)* _xyz->nX[0]*_xyz->nX[1]* (height*1e6);
  std::poisson_distribution<int> ng(rate);

  Ngrain =0;
  while (Ngrain == 0){
    Ngrain = ng(generator);
  } // end while 
} // end GenerateNgrain()
void BasePlate::GenerateVoxel()
{
  // this is a Voronoi tessellation
  int j1,j2,j3;
  GenerateNgrain();
  std::vector<std::vector<double>> Xsite(Ngrain, std::vector<double>(3));
  unsigned int sdloc;
  std::uniform_real_distribution<double> xrand(0.0,1.0);
  double Lx = _xyz->dX[0] * _xyz->nX[0];
  double Ly = _xyz->dX[1] * _xyz->nX[1];
  gNucleus.assign(Ngrain,0);
  for (int j=0; j<Ngrain;++j){
    Xsite[j] = {xrand(generator)*Lx,xrand(generator)*Ly, xrand(generator)*height};
    j3 = floor(Xsite[j][2]/_xyz->dX[2]);
    j2 = floor(Xsite[j][1]/_xyz->dX[1]);
    j1 = floor(Xsite[j][0]/_xyz->dX[0]);
    gNucleus[j] = _xyz->nX[0]*_xyz->nX[1]*j3 + (_xyz->nX[0])*j2 + j1;
  } // end for
  double xloc, yloc,zloc;
  std::vector<double> dsq(Ngrain);
  std::vector<int> jrange(Nzh*_xyz->nX[0]*_xyz->nX[1],0),jfind;
  std::iota(jrange.begin(),jrange.end(),0);
  std::set_intersection(_part->icellidLoc.begin(),_part->icellidLoc.begin()+_part->ncellLoc,
		     jrange.begin(),jrange.end(),std::back_inserter(jVals));
  X.assign(jVals.size(),0);
  jfind.resize(jVals.size(),0);
  for (int j=0;j<jVals.size();++j){jfind[j]= \
      std::distance(_part->icellidLoc.begin(),std::find(_part->icellidLoc.begin(),
			       _part->icellidLoc.begin()+_part->ncellLoc,jVals[j]));}
  for (int j=0; j < jVals.size();++j){
    j3 = floor(jVals[j]/(_xyz->nX[0]*_xyz->nX[1]));
    j2 = floor( (jVals[j] - _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
    j1 = jVals[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
    zloc = (j3+.5)*_xyz->dX[2];
    yloc = (j2+.5)*_xyz->dX[1];
    xloc = (j1+.5)*_xyz->dX[0];
    for (int j3=0; j3<Ngrain;++j3){
      dsq[j3] = std::pow(xloc - Xsite[j3][0],2)+std::pow(yloc-Xsite[j3][1],2)+std::pow(zloc-Xsite[j3][2],2);
    } // end for j3
    X[jfind[j]] = std::distance(dsq.begin(),std::min_element(dsq.begin(),dsq.end()))+1;
  } // for j

  std::default_random_engine g1(531231);
  //  for (int j=0;j<Ngrain;++j){sdloc.push_back(int(g1()/10000 + 1123));}
  sdloc = unsigned(double(g1())/double(g1.max())*pow(2.0,32.0));
  SampleOrientation sa;
  // randomly generate crystallographic orientations
  std::vector<double> aa;
  sa.GenerateSamples(Ngrain,sdloc,aa);
  cTheta.assign(Ngrain*4,0);
  for (int j1=0;j1<Ngrain;++j1){
    cTheta[4*j1] = aa[4*j1];
    cTheta[4*j1+1] =aa[4*j1+1];
    cTheta[4*j1+2] = aa[4*j1+2];
    cTheta[4*j1+3] = aa[4*j1+3];
  } // end for j1
} // end GenerateVoxel()
