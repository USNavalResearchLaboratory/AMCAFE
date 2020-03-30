// member function definitions for Grid.C

#include "BasePlateval1.h"
//#include "Grid.h"
#include "Partition.h"
#include "SampleOrientation.h"
#include "random"
#include <algorithm>
#include <numeric>
#include "iostream"
#include <math.h>

#include <mpi.h>


// constructor
BasePlate::BasePlate(const Grid & g, const double & bwIn, const double & bstdIn, const double & hIn,const double & muIn, Partition & part)
{
  _xyz = &g;
  _part = &part;
  mu = muIn;
  height = ceil(hIn/_xyz->dX[2])*_xyz->dX[2];
  Nzh = ceil(hIn/_xyz->dX[2]);
  height = Nzh*(_xyz->dX[2]);
  bwidth=bwIn;
  bInc = 1.53*bstdIn;
  // if you want seed to be different every time program runs
  //seed = std::chrono::system_clock::now().time_since_epoch().count();
  seed = 1234567;
  std::default_random_engine generator(seed); 
} // end constructor
void BasePlate::GenerateNgrain()
{
  double rate = mu* _xyz->dX[0]*_xyz->dX[1]* _xyz->nX[0]*_xyz->nX[1]* height;
  
  std::poisson_distribution<int> ng(rate);
  Ngrain =0;
  while (Ngrain == 0){
    Ngrain = ng(generator);
  } // end while 
} // end GenerateNgrain()

void BasePlate::GenerateVoxel()
{
  // generates a grain for each voxel but then rotates
  // to get the desired banded structure
  int j1,j2,j3,j,cc;
  Ngrain = _xyz->nX[0]*_xyz->nX[1]*Nzh;
  std::default_random_engine g1(531231);
  //  for (int j=0;j<Ngrain;++j){sdloc.push_back(int(g1()/10000 + 1123));}
  unsigned int sdloc;
  sdloc = unsigned(double(g1())/double(g1.max())*pow(2.0,32.0));
  SampleOrientation sa;
  std::vector<double> aa;
  sa.GenerateSamples(Ngrain,sdloc,aa);
  cTheta.assign(Ngrain*4,0);
  for (int j1=0;j1<Ngrain;++j1){
    cTheta[4*j1] = aa[4*j1];
    cTheta[4*j1+1] =aa[4*j1+1];
    cTheta[4*j1+2] = aa[4*j1+2];
    cTheta[4*j1+3] = aa[4*j1+3];
  } // end for j1
  double Lx = _xyz->dX[0] * _xyz-> nX[0];
  double Ly = _xyz->dX[1] * _xyz->nX[1];
  gNucleus.assign(Ngrain,0);
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
    X[jfind[j]] = jVals[j]+1;
  } // for j                                                                                                                                       
  for (int j3=0;j3<Nzh;++j3){
    for (int j2=0;j2<_xyz->nX[1];++j2){
      for (int j1=0;j1<_xyz->nX[0];++j1){
	j = _xyz->nX[0]*_xyz->nX[1]*j3+_xyz->nX[0]*j2+j1;
	gNucleus[j] = j;
      }
    }
  } // end for
  // generate crystallographic orientations such that within center of melt pool
  // bands, the <001> is pointing upward and outside bands the <101> is pointing
  // upwards
  // strategy: withinn band:find angle from reference [001] to closest local <001>.
  // then apply rotation of that angle about axis perpendicular to these two vectors.
  // outside band: do same thing for local <101>

  std::vector<std::vector<double> > axl,rRot(3,std::vector<double>(3,0));
  std::vector<double> axr(3,0),ax(3,0),axg(3,0),q0(4,0),qr(4,0),q2(4,0),tmp(3,0);
  double eta,y,omega,th2,theta;
  for (int j=0;j<Ngrain;++j){
    omega = cTheta[4*j];
    ax[0]=cTheta[4*j+1];
    ax[1]=cTheta[4*j+2];
    ax[2]=cTheta[4*j+3];
    rRot[0][0] = cos(omega) + pow(ax[0],2.0)*(1-cos(omega));
    rRot[0][1] = ax[0]*ax[1]*(1-cos(omega)) - ax[2]*sin(omega);
    rRot[0][2] = ax[0]*ax[2]*(1-cos(omega)) + ax[1]*sin(omega);
    rRot[1][0] = ax[0]*ax[1]*(1-cos(omega)) + ax[2]*sin(omega);
    rRot[1][1] = cos(omega) + pow(ax[1],2.0)*(1-cos(omega));
    rRot[1][2] = ax[1]*ax[2]*(1-cos(omega)) - ax[0]*sin(omega);
    rRot[2][0] = ax[2]*ax[0]*(1-cos(omega)) - ax[1]*sin(omega);
    rRot[2][1] = ax[2]*ax[1]*(1-cos(omega)) + ax[0]*sin(omega);
    rRot[2][2] = cos(omega) + pow(ax[2],2.0)*(1-cos(omega));            
    j3 = floor(gNucleus[j]/(_xyz->nX[0]*_xyz->nX[1]));
    j2 = floor( (gNucleus[j] -  _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
    //j1 = gNucleus[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
    // find out if within a band
    y = _xyz->dX[1]*j2;
    eta = round( y / bInc);
//    if ( (eta*bInc - bwidth/2 <= y) && (y<= eta*bInc + bwidth/2)){
    if ( std::fabs(eta*bInc - y)<=bwidth/2.0){
      axl.assign(6,std::vector<double>(3,0));
      axl[0] = {1,0,0};
      axl[1] = {-1,0,0};
      axl[2] = {0,1,0};
      axl[3] = {0,-1,0};
      axl[4] = {0,0,1};
      axl[5] = {0,0,-1};
    } else {
      axl.assign(12,std::vector<double>(3,0));
      axl[0] = {1/pow(2,.5),1/pow(2,.5),0};
      axl[1] = {-1/pow(2,.5),1/pow(2,.5),0};
      axl[2] = {1/pow(2,.5),-1/pow(2,.5),0};
      axl[3] = {-1/pow(2,.5),-1/pow(2,.5),0};
      axl[4] = {1/pow(2,.5),0,1/pow(2,.5)};
      axl[5] = {-1/pow(2,.5),0,1/pow(2,.5)};
      axl[6] = {1/pow(2,.5),0,-1/pow(2,.5)};
      axl[7] = {-1/pow(2,.5),0,-1/pow(2,.5)};
      axl[8] = {0,1/pow(2,.5),1/pow(2,.5)};
      axl[9] = {0,-1/pow(2,.5),1/pow(2,.5)};
      axl[10] = {0,1/pow(2,.5),-1/pow(2,.5)};
      axl[11] = {0,-1/pow(2,.5),-1/pow(2,.5)};
    }
    theta = 1000;
    for (int j1=0;j1<axl.size();++j1){
      axg[0] = rRot[0][0]*axl[j1][0]+rRot[1][0]*axl[j1][1]+rRot[2][0]*axl[j1][2];
      axg[1] = rRot[0][1]*axl[j1][0]+rRot[1][1]*axl[j1][1]+rRot[2][1]*axl[j1][2];
      axg[2] = rRot[0][2]*axl[j1][0]+rRot[1][2]*axl[j1][1]+rRot[2][2]*axl[j1][2];
      th2 = acos( axg[2]);
      if (th2<theta){
	theta = th2;
	axr[0] = axg[1]/pow(pow(axg[0],2.0)+pow(axg[1],2.0),.5);
	axr[1] = -axg[0]/pow(pow(axg[0],2.0)+pow(axg[1],2.0),.5);
	axr[2] = 0.0;	  
      }
    }
    // convert axis angle to quat
    q0[0] = cos(omega/2.0);
    q0[1] = -ax[0]*sin(omega/2.0);
    q0[2] = -ax[1]*sin(omega/2.0);
    q0[3] = -ax[2]*sin(omega/2.0);
    qr[0] = cos(theta/2.0);
    qr[1] = axr[0]*sin(theta/2.0);
    qr[2] = axr[1]*sin(theta/2.0);
    qr[3] = axr[2]*sin(theta/2.0);
    // multiple 2 quats and then convert back to axis-angle
    tmp[0] = qr[2]*q0[3]-q0[2]*qr[3];
    tmp[1] = qr[3]*q0[1]-q0[3]*qr[1];
    tmp[2] = qr[1]*q0[2]-q0[1]*qr[2];
    q2[0] = qr[0]*q0[0] - (qr[1]*q0[1]+qr[2]*q0[2]+qr[3]*q0[3]);
    q2[1] =  (q0[0]*qr[1] + qr[0]*q0[1] + tmp[0]);
    q2[2] =  (q0[0]*qr[2] + qr[0]*q0[2] + tmp[1]);
    q2[3] =  (q0[0]*qr[3] + qr[0]*q0[3] + tmp[2]);
    omega = 2*acos(q2[0]);
    th2 = 1.0 / 
      pow( pow(q2[1],2.0)+pow(q2[2],2.0)+pow(q2[3],2.0),.5);
    ax[0] = q2[1]*th2;
    ax[1] = q2[2]*th2;
    ax[2] = q2[3]*th2;    
    cTheta[4*j]=omega;
    cTheta[4*j+1]=ax[0];
    cTheta[4*j+2]=ax[1];
    cTheta[4*j+3]=ax[2];    
  }
} // end GenerateVoxel()



