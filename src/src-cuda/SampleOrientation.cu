// member function definitions for Grid.C
#include "SetPrecision.cuh"
#include "SampleOrientation.cuh"
#include "random"
#include <numeric>
#include <algorithm>
#include <math.h>
#include <curand_kernel.h>
#include <math_constants.h>


__device__ void GetPermutation(const double *xyz, int *P)
{
  double x,y,z;
  x= xyz[0];
  y = xyz[1];
  z = xyz[2];
  // the 6 if statements refer to the 6 pyramids in rosca2014new 
  // (in AM literature folder). Order is C1,C2,C3,C4,C5,C6
  if ( (fabs(x) <= z) && (fabs(y) <= z) ){ 
    P[0]= 1; P[1]= 2; P[2]= 3;
  }
  if ( (z <= -fabs(x)) && (z<= -fabs(y)) ){ 
    P[0]= -1; P[1]= 2; P[2]= -3;
  }
  if ( (fabs(x)<=y) && (fabs(z)<=y)){
    P[0]= 1; P[1]= -3; P[2]= 2;
  }
  if ( (y<=-fabs(x)) && (y<= -fabs(z))){
    P[0]= 1; P[1]= 3; P[2]= -2;
  }
  if ( (fabs(z)<x) && (fabs(y)<= x)){
    P[0]= -3; P[1]= 2; P[2]= 1;
  }
  if ((x<=-fabs(z)) && (x<=-fabs(y)) ){
    P[0]= 3; P[1]= 2; P[2]= -1;
  }
} // end GetPermutation

__device__ void ComputeMap1(const double &a, const double& aprime, const double *xyz, double *map1)
{
 
  map1[0] = a/aprime*xyz[0];
  map1[1] = a/aprime*xyz[1];
  map1[2] = a/aprime*xyz[2];
} // end ComputeMap1

__device__ void ComputeMap2(const double &beta,const double *xyz, double *map2)
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
      Tcurv1=pow(2.0,1.0/2.0)*cos(y*CUDART_PI/12/x) - 1;
      Tcurv2=pow(2.0,1.0/2.0)*sin(y*CUDART_PI/12/x);
      Tcurv1*=pow(2.0,1.0/4.0)*x/(beta*pow( pow(2.0,1.0/2.0)-cos(y*CUDART_PI/12/x),1.0/2.0));
      Tcurv2*=pow(2.0,1.0/4.0)*x/(beta*pow( pow(2.0,1.0/2.0)-cos(y*CUDART_PI/12/x),1.0/2.0));
    }
  }
  if (fabs(x) <=fabs(y)) {
    if (fabs(y) < 1e-10){
      Tcurv1=0.0;
      Tcurv2=0.0;
    } else {
      Tcurv1=pow(2.0,1.0/2.0)*sin(x*CUDART_PI/12/y);
      Tcurv2=pow(2.0,1.0/2.0)*cos(x*CUDART_PI/12/y) -1;
      Tcurv1*=pow(2.0,1.0/4.0)*y/(beta*pow( pow(2.0,1.0/2.0)-cos(x*CUDART_PI/12/y),1.0/2.0));
      Tcurv2*=pow(2.0,1.0/4.0)*y/(beta*pow( pow(2.0,1.0/2.0)-cos(x*CUDART_PI/12/y),1.0/2.0));
    }
  }
  map2[0] = Tcurv1; map2[1] = Tcurv2; map2[2] = z;
} // end ComputeMap2

__device__ void ComputeMap3(const double *xyz, double *map3)
{
  double x,y,z,m1,m2,m3;
  x=xyz[0];
  y=xyz[1];
  z=xyz[2];
  if (fabs(z) < 1e-10){
    map3[0] = 0; map3[1] = 0; map3[2] = 0;
  } else {
    m1 = pow(1-(pow(x,2.)+pow(y,2.))*CUDART_PI/(24*pow(z,2.)),1.0/2.0)*x;
    m2 = pow(1-(pow(x,2.)+pow(y,2.))*CUDART_PI/(24*pow(z,2.)),1.0/2.0)*y;
    m3 = pow(6/CUDART_PI,1.0/2.0)*z - (pow(x,2.)+pow(y,2.))*pow(CUDART_PI,1.0/2.0)/(pow(24,1.0/2.0)*z);
    map3[0]=m1; map3[1]=m2; map3[2]=m3;
  }
} // end ComputeMap3

__device__  void ComputeMeasure(const double &t,double & f)
{
  float r=1.;
  if ( fabs((fabs(t)-1.0)) <1e-10){
    f=1.0;
  } else {
    if (t >=0){
      f = r/pow(pow(r,2.)-pow(t,2),1.0/2.0)*
	pow(3.0/2.0*(acos(t/r)-t/pow(r,2)*pow(pow(r,2)-pow(t,2),1.0/2.0)),1.0/3.0);
    } else {
      f = r/pow(pow(r,2.)-pow(t,2.),1.0/2.0)*
	pow(3.0/2.0*(CUDART_PI +t/pow(r,2.)*pow(pow(r,2.)-pow(t,2),1.0/2.0) -acos(t/r)),1.0/3.0);     
    }
  }
} // end ComputeMeasure

__device__ void GenerateSamples(const int Nsample,unsigned int seedL, int subsq, curandState_t &s1, double *axisAngle)
{
  double a,rhoP,t,omega,tnorm,pmax,coeffA[7],beta,aprime,
    xyz0[3],xyz[3],xyz2[3],xyzb[3],map1[3],map2[3],map3[3];
  float r=1.;
  double R1 = pow(3.0*CUDART_PI/4.0,1.0/3.0)*r;
  int Perm[]={1,2,3};
  a  = pow(2.0*CUDART_PI/3.0,1.0/2.0)*R1;
  aprime = pow(4.0*CUDART_PI/3.0,1.0/3.0)*R1;
  coeffA[0]= -0.500009615; coeffA[1]= -0.024866061; coeffA[2]= -0.004549382;
  coeffA[3]= 0.000511867; coeffA[4]= -0.001650083; coeffA[5]=  0.000759335; coeffA[6]= -0.000204042;
  beta = pow(CUDART_PI/6.0,1.0/2.0);
  curand_init(seedL,subsq,0,&s1);
  for (int jn=0;jn<Nsample;++jn){
    xyz0[0]= curand_uniform(&s1)-.5;
    xyz0[1]= curand_uniform(&s1)-.5;
    xyz0[2]= curand_uniform(&s1)-.5;
    xyz[0] = aprime*xyz0[0];
    xyz[1] = aprime*xyz0[1];
    xyz[2] = aprime*xyz0[2];
    GetPermutation(xyz,Perm);
    xyz2[0] = Perm[0]/abs(Perm[0])*xyz[abs(Perm[0])-1];
    xyz2[1] = Perm[1]/abs(Perm[1])*xyz[abs(Perm[1])-1];
    xyz2[2] = Perm[2]/abs(Perm[2])*xyz[abs(Perm[2])-1];
    // line below is equivalent to ComputeMap1
    ComputeMap1(a,aprime,xyz2,map1);
    ComputeMap2(beta,map1,map2);
    ComputeMap3(map2,map3);
    xyzb[abs(Perm[0])-1] = Perm[0]/abs(Perm[0])*map3[0];
    xyzb[abs(Perm[1])-1] = Perm[1]/abs(Perm[1])*map3[1];
    xyzb[abs(Perm[2])-1] = Perm[2]/abs(Perm[2])*map3[2];
    pmax=max(max(fabs(xyz0[0]),fabs(xyz0[1])),fabs(xyz0[2]));
    rhoP=R1*pmax/.5;
    t = 1;
    for (int j=1;j<8;++j){t+=coeffA[j-1]*pow(rhoP,2.0*j);}
    omega = 2.0*acos(t);    
    //t=pow(xyzb[0],2.0)+pow(xyzb[1],2.0)+pow(xyzb[2],2.0);
    tnorm = pow(pow(xyzb[0],2.0)+pow(xyzb[1],2.0)+pow(xyzb[2],2.0),.5);
    axisAngle[4*jn] = omega;
    axisAngle[4*jn+1] = xyzb[0]/tnorm;
    axisAngle[4*jn+2] = xyzb[1]/tnorm;
    axisAngle[4*jn+3] = xyzb[2]/tnorm;
  } // for (int jn =0...)

}

  
