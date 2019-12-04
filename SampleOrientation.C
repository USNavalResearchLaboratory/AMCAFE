// member function definitions for Grid.C

#include "SampleOrientation.h"
#include "random"
#include <numeric>
#include <algorithm>
#include "iostream"
#include <math.h>

// constructor
SampleOrientation::SampleOrientation()
{
  // if you want seed to be different every time program runs
  //seed = std::chrono::system_clock::now().time_since_epoch().count();
  // seed = 1234567;
  r = 1;
  Ngrid = 100;
  NTgrid = 2*Ngrid+1;
  R1 = pow(3.0*M_PI/4.0,1.0/3.0)*r;
  a  = pow(2.0*M_PI/3.0,1.0/2.0)*R1;
  aprime = pow(4.0*M_PI/3.0,1.0/3.0)*R1;
  coeffA = {-0.500009615, -0.024866061,-0.004549382,0.000511867,
	    -0.001650083, 0.000759335,-0.000204042};
  beta = pow(M_PI/6.0,1.0/2.0);
} // end constructor

void SampleOrientation::GetPermutation(std::vector<double> & xyz, std::vector<int> &P)
{
  double x,y,z;
  x= xyz[0];
  y = xyz[1];
  z = xyz[2];
  // the 6 if statements refer to the 6 pyramids in rosca2014new 
  // (in AM literature folder). Order is C1,C2,C3,C4,C5,C6
  if ( (fabs(x) <= z) && (fabs(y) <= z) ){ 
    P = {1,2,3};
  }
  if ( (z <= -fabs(x)) && (z<= -fabs(y)) ){ 
    P = {-1,2,-3};
  }
  if ( (fabs(x)<=y) && (fabs(z)<=y)){
    P = {1,-3,2};
  }
  if ( (y<=-fabs(x)) && (y<= -fabs(z))){
    P = {1,3,-2};
  }
  if ( (fabs(z)<x) && (fabs(y)<= x)){
    P = {-3,2,1};
  }
  if ((x<=-fabs(z)) && (x<=-fabs(y)) ){
    P = {3,2,-1};
  }
} // end GetPermutation

void SampleOrientation::ComputeMap1(std::vector<double> & xyz, std::vector<double> & map1)
{
  map1.assign(3,0);
  map1[0] = a/aprime*xyz[0];
  map1[1] = a/aprime*xyz[1];
  map1[2] = a/aprime*xyz[2];
} // end ComputeMap1

void SampleOrientation::ComputeMap2(std::vector<double> & xyz, std::vector<double> & map2)
{
  double x,y,z,Tcurv1,Tcurv2;
  x=xyz[0];
  y=xyz[1];
  z=xyz[2];
  if ( fabs(y) <= fabs(x) ){
    if (fabs(x) < 1e-10){
      Tcurv1=0.0;
      Tcurv2=0.0;
    } else {
      Tcurv1=pow(2.0,1.0/2.0)*cos(y*M_PI/12/x) - 1;
      Tcurv2=pow(2.0,1.0/2.0)*sin(y*M_PI/12/x);
      Tcurv1*=pow(2.0,1.0/4.0)*x/(beta*pow( pow(2.0,1.0/2.0)-cos(y*M_PI/12/x),1.0/2.0));
      Tcurv2*=pow(2.0,1.0/4.0)*x/(beta*pow( pow(2.0,1.0/2.0)-cos(y*M_PI/12/x),1.0/2.0));
    }
  }
  if (fabs(x) <=fabs(y)) {
    if (fabs(y) < 1e-10){
      Tcurv1=0.0;
      Tcurv2=0.0;
    } else {
      Tcurv1=pow(2.0,1.0/2.0)*sin(x*M_PI/12/y);
      Tcurv2=pow(2.0,1.0/2.0)*cos(x*M_PI/12/y) -1;
      Tcurv1*=pow(2.0,1.0/4.0)*y/(beta*pow( pow(2.0,1.0/2.0)-cos(x*M_PI/12/y),1.0/2.0));
      Tcurv2*=pow(2.0,1.0/4.0)*y/(beta*pow( pow(2.0,1.0/2.0)-cos(x*M_PI/12/y),1.0/2.0));
    }
  }
  map2 = {Tcurv1,Tcurv2,z};
} // end ComputeMap2

void SampleOrientation::ComputeMap3(std::vector<double> & xyz, std::vector<double> & map3)
{
  double x,y,z,m1,m2,m3;
  x=xyz[0];
  y=xyz[1];
  z=xyz[2];
  if (fabs(z) < 1e-10){
    map3 = {0,0,0};
  } else {
    m1 = pow(1-(pow(x,2)+pow(y,2))*M_PI/(24*pow(z,2)),1.0/2.0)*x;
    m2 = pow(1-(pow(x,2)+pow(y,2))*M_PI/(24*pow(z,2)),1.0/2.0)*y;
    m3 = pow(6/M_PI,1.0/2.0)*z - (pow(x,2)+pow(y,2))*pow(M_PI,1.0/2.0)/pow(24*z,1.0/2.0);
    map3={m1,m2,m3};
  }
} // end ComputeMap3

void SampleOrientation::ComputeMeasure(double &t,double & f)
{
  if ( fabs((fabs(t)-1.0)) <1e-10){
    f=1.0;
  } else {
    if (t >=0){
      f = r/pow(pow(r,2)-pow(t,2),1.0/2.0)*
	pow(3.0/2.0*(acos(t/r)-t/pow(r,2)*pow(pow(r,2)-pow(t,2),1.0/2.0)),1.0/3.0);
    } else {
      f = r/pow(pow(r,2)-pow(t,2),1.0/2.0)*
	pow(3.0/2.0*(M_PI +t/pow(r,2)*pow(pow(r,2)-pow(t,2),1.0/2.0) -acos(t/r)),1.0/3.0);
     
    }
  }
} // end ComputeMeasure

void SampleOrientation::GenerateSamples(const int Nsample,unsigned int seedL, std::vector<double> &axisAngle)
{
  int pind;
  double rhoP,dx,t,omega,fmeas,tnorm;
  std::uniform_real_distribution<double> xrand(0.0,1.0);
  std::vector<double> xyz(3),xyz2(3),xyzb(3),xyzs(3),map1(3),map2(3),map3(3);
  std::vector<int> Perm(3),pid(3);
  axisAngle.assign(Nsample*4,0.0);
  dx = aprime/2.0/Ngrid;
  std::default_random_engine g1(seedL);
  for (int jn=0;jn<Nsample;++jn){
    xyz = {aprime*(xrand(g1)-.5),aprime*(xrand(g1)-.5),
	   aprime*(xrand(g1)-.5)};
    Perm={1,2,3};
    GetPermutation(xyz,Perm);
    for (int j=0;j<3;++j){xyz2[j] = Perm[j]/abs(Perm[j])*xyz[abs(Perm[j])-1];}
    ComputeMap1(xyz2,map1);
    ComputeMap2(map1,map2);
    ComputeMap3(map2,map3);

    for (int j=0;j<3;++j){xyzb[abs(Perm[j])-1] = Perm[j]/abs(Perm[j])*map3[j];}
    pid = {int(round(xyz[0]/dx)),int(round(xyz[1]/dx)),int(round(xyz[2]/dx))};
    pind = *std::max_element(pid.begin(),pid.end());
    rhoP=R1*pind/Ngrid;
    t = 1;
    for (int j=1;j<8;++j){t+=coeffA[j-1]*pow(rhoP,2.0*j);}
    omega = 2.0*acos(t);
    ComputeMeasure(t,fmeas);
    // dont think fmeas always produces normalized axis
    t=pow(xyzb[0],2.0)+pow(xyzb[1],2.0)+pow(xyzb[2],2.0);
    tnorm = pow(pow(xyzb[0],2.0)+pow(xyzb[1],2.0)+pow(xyzb[2],2.0),.5);
    xyzs = {xyzb[0]/tnorm,xyzb[1]/tnorm,xyzb[2]/tnorm};
    axisAngle[4*jn] = omega;
    axisAngle[4*jn+1] = xyzs[0];
    axisAngle[4*jn+2] = xyzs[1];
    axisAngle[4*jn+3] = xyzs[2];
  } // for (int jn =0...)


}

  
