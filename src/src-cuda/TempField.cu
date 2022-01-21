// member functions for TempField

#include "Grid.cuh"
#include "TempField.cuh"
#include "fstream"
// #include "iostream"
#include "numeric"
#include "math.h"
#include <float.h>
#include <algorithm>

__global__ void calcTemptInd(const Grid *g,TempField *temp)
{
  temp->tInd = int(round(g->time/temp->DelT));
}

__global__ void analyticTempCalcPart1(Grid *dg,TempField *dtempF, int *ispv, 
				 double *dlcoor, double *dtempOut, 
				 int Ntot)
{

  //computes temp field based on Schwalbach et al
  int j1,j2,j3,tid=threadIdx.x+blockDim.x*blockIdx.x,js,stride,nX[2],i1;
  double x0,y0,x,y,z,dsq,xi,xp,yp,dirp,zp,rij1[3],a1m[6],dx,LX[2];
  dx=dg->dX[0];
  nX[0]=dg->nX[0];nX[1]=dg->nX[1];
  LX[0]=nX[0]*dx;;LX[1]=nX[1]*dx;
  stride = gridDim.x*blockDim.x;
  xi = dg->tL*(1+FLT_EPSILON);
  for (int j=0;j<6;++j){a1m[j] = dtempF->a1[j];}
  // x,y,z spatial location of source (in grid reference frame)
  x = dlcoor[2*ispv[dg->isp]]-dtempF->offset[0];
  y = dlcoor[2*ispv[dg->isp]+1]-dtempF->offset[1];
  z = dg->ilaserLoc*dx-dtempF->offset[2];
  // initialize temp field
  js=tid;
  x0=dtempF->T0;
  while (js<Ntot){
    dtempOut[js]=x0;
    js+=stride;
  }
  js=tid;
  i1= dg->isp - ((dg->isp)/(dg->Nsd))*(dg->Nsd);
  while (js < Ntot){
    j3 = js/(nX[0]*nX[1]);
    j2 = (js -nX[0]*nX[1]*j3)/nX[0];
    j1 = js - nX[0]*nX[1]*j3 - nX[0]*j2;
    x0 = (double(j1)+.5)*(dx);
    y0 = (double(j2)+.5)*(dx);
    zp = (double(j3)+.5)*(dx);
    if (zp<=z){
      xp=cos(dg->gth)*(x0-LX[0]/2.)+
	sin(dg->gth)*(y0-LX[1]/2.)+LX[0]/2.;
      yp=-sin(dg->gth)*(x0-LX[0]/2.) +
	cos(dg->gth)*(y0-LX[1]/2.)+LX[1]/2.;
      if (i1==0){
	dirp=(dlcoor[2*ispv[dg->isp+1]]-dlcoor[2*ispv[dg->isp]]);
      } else {
	dirp=(dlcoor[2*ispv[dg->isp]]-dlcoor[2*ispv[dg->isp-1]]);
      }
      rij1[0] = xp-x;
      rij1[1] = yp-y;
      rij1[2] = zp-z;
      if (dirp*rij1[0]>0){
	//xp,yp is in front of laser
	dsq = pow(rij1[1]/a1m[1],2.0)+pow(rij1[2]/a1m[2],2.0);
	if (dsq<1.0 && (fabs(rij1[0])<dtempF->bmDX[0]) ){
	  dtempOut[js]=xi;
	} else {
	  dtempOut[js] = dg->tS;
	}
      } else {
	dsq = pow(rij1[0]/a1m[3],2.0)+pow(rij1[1]/a1m[4],2.0)+pow(rij1[2]/a1m[5],2.0);
	if (dsq<1.0){
	  dtempOut[js]=xi;
	} else {
	  dtempOut[js] = dg->tS;
	}
      } // if (dirp*rij1[0]>0...
    } // if (zp<=z
    js+=stride;
  } // while (js<Ntot
  // check if last simulation of scan                   

}

__global__ void analyticTempCalcPart2(Grid *dg, TempField *dtempF,int *ispv,double *dlcoor2,
				      double *dtempOut, int Ntot)
{
  // check if last simulation of scan
  int tid=threadIdx.x+blockDim.x*blockIdx.x, js, stride=blockDim.x*gridDim.x;
  double tmelt=dg->tL,x,y;
  extern __shared__ volatile bool icheck[];
  x=dlcoor2[2*ispv[dg->isp]]-dtempF->offset[0];
  y=dlcoor2[2*ispv[dg->isp]+1]-dtempF->offset[1];
  if (x<dg->gbox[0] || x>dg->gbox[1] || y<dg->gbox[2] || y>dg->gbox[3]){
    icheck[threadIdx.x]=1;
    js=tid;
    while ((js<Ntot) & (icheck[threadIdx.x])){
      icheck[threadIdx.x]= dtempOut[js]<tmelt;      
      js+=stride;
    }
    js = (dg->isp+1) - ((dg->isp+1)/(dg->Nsd)) * (dg->Nsd); // equiv to fmod
    if (!icheck[threadIdx.x] || js==0 ){dg->inewscanflg=1;} 
  } // if (x<box[0...
}


// constructor
TempField::TempField(const Grid &g)
{
  tInd = 0;
  bmV = g.bmV;
  patternID = g.patternID;
  T0 = g.T0;
  DelT = g.bmDelT;
  bmDX[0]=DelT*bmV; bmDX[1]=g.bhatch; bmDX[2]=g.layerT; // (SD, TD, BD) per layer
  offset[0]=g.offset[0];
  offset[1]=g.offset[1];
  offset[2]=g.offset[2];
  ispvec = (int*)malloc(g.NpT*sizeof(int));
  if (patternID==1 || patternID==3){
    std::iota(ispvec,ispvec+g.NpT,0);
  } // if (patternID==1...
  if (patternID==2 || patternID==4){
    int k;
    for (int j=0;j<g.Ntd;++j){
      k=g.Nsd*j;
      if (fmod(j,2)==0){
	std::iota(ispvec+k,ispvec+k+g.Nsd,k);
      } // if (fmod(j,2)==0...
      else {
	for (int j1=0;j1<g.Nsd;++j1){
	  ispvec[k+j1]=k+g.Nsd-1-j1;
	} // for (int j1=0...
      } // else (fmod(j,2)==0...
    } // for (int j=0...
  } // if (patternID==2...
  int Ntot = g.nX[0]*g.nX[1]*g.nX[2];
  TempCurr = (double*)malloc(Ntot*sizeof(double));
} // end TempField

void TempField::InitializeAnalytic(const Grid &g)
{
  /*
    This this 2 double ellipsoids (one encompassing another) to represent the temperature
    field. This approach is used in Rogers_Madison_Tikare, Comp Mat Sci, 2017.
    SCAN INFORMATION
    patternID=0: scan in +X direction only
    patternID=1: scan in alternating +/- X direction (same direction in
		 layer above and below
    patternID=2: scan in alternating +/- X direction (different direction 
		 in layer above and below
    patternID=3: scan in +X direction for layer i and +Y direction for 
		 layer i+1 ...
    patternID=4: scan in alternating +/- X direction for layer i and 
		 alternating +/- Y direction for layer i+1
   */
  a1[0] = g.meltparam[0];
  a1[1] = g.meltparam[2];
  a1[2] = g.meltparam[3];
  a1[3] = g.meltparam[1];
  a1[4] = a1[1];
  a1[5] = a1[2];
} // end InitializeAnalytic 

void TempField::AnalyticalTempCalcMacro(Grid *dg,TempField *dtempF, double *dtempOut, 
					double *dlcoor,double *dlcoor2,int *ispv,
					int &nThreads, int &nBlocks, int &Ntot)
{

  analyticTempCalcPart1<<<nBlocks,nThreads>>>(dg,dtempF,ispv,dlcoor,dtempOut,Ntot);
//  analyticTempCalcPart2<<<nBlocks,nThreads>>>(dg,dtempF,ispv,dlcoor2,dtempOut,Ntot);
					
}
