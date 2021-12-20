// define member functions of VoxelCA

#include "Grid.cuh"
#include "VoxelsCA.cuh"
#include "iostream"
#include "fstream"
#include "sstream"
#include <math.h>
#include <algorithm>
#include <random>
#include <ctime>
#include "SampleOrientation.cuh"
#include <chrono>
#include <adios2.h>

// constructor
VoxelsCA::VoxelsCA(Grid &g)
{ 
  int Ntot = g.nX[0]*g.nX[1]*g.nX[2];
  gID = (int*)malloc(Ntot*sizeof(int));
  vState = (int*)malloc(Ntot*sizeof(int));
  extents = (double*)malloc(Ntot*sizeof(double));
  centroidOct = (double*)malloc(3*Ntot*sizeof(double));
  memset(gID,0,Ntot*sizeof(int));
  memset(vState,0,Ntot*sizeof(int));
  memset(extents,0,Ntot*sizeof(double));
  memset(centroidOct,0,3*Ntot*sizeof(double));
  seed0= 2132512;
  seed1=2912351;
  // establishes ineighID and ineighptr for convertSolid1 
  /*
  int cc=0;
  std::vector<int> neigh;

  
  ineighptr.assign(Ntot1+1,0);
  for (int j=0; j < Ntot1;++j){
    g.ComputeNeighborhood(_part->icellidLoc[j],g.neighOrder,neigh);
    cc+=neigh.size();
    ineighptr[j+1]=cc;
  }
  ineighID.assign(cc,0);
  cc=0;
  for (int j=0;j<Ntot1;++j){
    g.ComputeNeighborhood(_part->icellidLoc[j],g.neighOrder,neigh);
    for (int j1=0;j1<neigh.size();++j1){
      ineighID[cc] = std::distance(_part->icellidLoc.begin(),
			 std::find(_part->icellidLoc.begin(),_part->icellidLoc.end(),neigh[j1]));
      cc+=1;
    }
  }
  */
} // end constructor

void VoxelsCA::WriteToHDF1(const std::string &filename, const Grid &g, const double *tempcurr)
{
  // writes gID, vState, cTheta per voxel
  int Ntot = g.nX[0]*g.nX[1]*g.nX[2],icase;
  std::string hdf5Filename = filename + ".h5";
  std::vector< float> TempOut(Ntot,0),IPFmapBD(3*Ntot,0), IPFmapx(3*Ntot,0), IPFmapy(3*Ntot,0),cth(4*nGrain,0);
  double vBD[3]={0.0,0.0,1.0},omega,ax[3],vCD[3],rRot[3][3],
    vX[3]={1.0,0.0,0.0},vY[3]={0.0,1.0,0.0},xp,yp,x0,y0,m,a,b,c,H,S,V,sMax,ff,p,q,t;
  std::vector<std::vector<double>> triPts(2,std::vector<double>(3,0));
  triPts[0][0]=0.0;
  triPts[0][1]=2./pow(2,.5)/(1.+1./pow(2,.5));
  triPts[0][2]=2./pow(3,.5)/(1.+1./pow(3,.5));
  triPts[1][0]=0.0;
  triPts[1][1]=0.0;
  triPts[1][2]=2./pow(3,.5)/(1.+1./pow(3,.5));
  m=tan(1./2.*atan2(triPts[1][2],triPts[0][2]));
  a=pow(pow(triPts[1][2]-triPts[1][1],2.)+pow(triPts[0][2]-triPts[0][1],2.),.5);
  b=pow(pow(triPts[0][1],2.)+pow(triPts[1][1],2.) ,0.5);
  c=pow(pow(triPts[0][2],2.)+pow(triPts[1][2],2.),0.5);
  y0=1./2.*pow((b+c-a)*(c+a-b)*(a+b-c)/(a+b+c),.5);
  x0=y0/m;
  sMax=pow(pow(x0,2.)+pow(y0,2.),.5);
  for (int j=0;j<Ntot;++j){
    TempOut[j] = tempcurr[j];
    if (gID[j]<1){
      IPFmapBD[3*j] = 0.0;
      IPFmapBD[3*j+1] = 0.0;
      IPFmapBD[3*j+2] = 0.0;
      IPFmapx[3*j] = 0.0;
      IPFmapx[3*j+1] = 0.0;
      IPFmapx[3*j+2] = 0.0;
      IPFmapy[3*j] = 0.0;
      IPFmapy[3*j+1] = 0.0;
      IPFmapy[3*j+2] = 0.0;
    } else {
      omega = cTheta[4*(gID[j]-1)];
      ax[0]= cTheta[4*(gID[j]-1)+1];
      ax[1]= cTheta[4*(gID[j]-1)+2];
      ax[2]= cTheta[4*(gID[j]-1)+3];
      // matrix is local->global; need to multiply by transpose for global->local            
      rRot[0][0] = cos(omega) + pow(ax[0],2.0)*(1-cos(omega));
      rRot[0][1] = ax[0]*ax[1]*(1-cos(omega)) - ax[2]*sin(omega);
      rRot[0][2] = ax[0]*ax[2]*(1-cos(omega)) + ax[1]*sin(omega);
      rRot[1][0] = ax[0]*ax[1]*(1-cos(omega)) + ax[2]*sin(omega);
      rRot[1][1] = cos(omega) + pow(ax[1],2.0)*(1-cos(omega));
      rRot[1][2] = ax[1]*ax[2]*(1-cos(omega)) - ax[0]*sin(omega);
      rRot[2][0] = ax[2]*ax[0]*(1-cos(omega)) - ax[1]*sin(omega);
      rRot[2][1] = ax[2]*ax[1]*(1-cos(omega)) + ax[0]*sin(omega);
      rRot[2][2] = cos(omega) + pow(ax[2],2.0)*(1-cos(omega));
      vCD[0] = std::fabs(rRot[0][0]*vBD[0]+rRot[1][0]*vBD[1]+rRot[2][0]*vBD[2]);
      vCD[1] = std::fabs(rRot[0][1]*vBD[0]+rRot[1][1]*vBD[1]+rRot[2][1]*vBD[2]);
      vCD[2] = std::fabs(rRot[0][2]*vBD[0]+rRot[1][2]*vBD[1]+rRot[2][2]*vBD[2]);
      std::sort(vCD,vCD+3);
      std::swap(vCD[0],vCD[1]);
      vCD[2] = std::min(vCD[2],1.0);
      xp = 2.*vCD[0]/(1.+vCD[2]);
      yp = 2.*vCD[1]/(1.+vCD[2]);
      H=atan( (yp-y0)/(xp-x0))*180./M_PI;
      xp < x0 ? H+=180: H;
      H = H+240-atan((triPts[1][2]-y0)/(triPts[0][2]-x0))*180/M_PI;
      V=1.;
      S=pow(pow(xp-x0,2.)+pow(yp-y0,2.),.5);
      S=S/sMax*0.8 + 0.2;
      H>=360.0 ? H=0.0 : H;
      icase = floor(H/60.0);
      ff= H/6.0 - icase;
      p=V*(1.-S);
      q=V*(1.-S*ff);
      t=V*(1.-(S*(1.-ff)));
      if (S<=0.0){
	IPFmapBD[3*j]=0.0;
	IPFmapBD[3*j+1]=0.0;
	IPFmapBD[3*j+2]=0.0;
      } else {
	if (icase==0){
	  IPFmapBD[3*j]=V;
	  IPFmapBD[3*j+1]=t;
	  IPFmapBD[3*j+2]=p;
	}
	if (icase==1){
	  IPFmapBD[3*j]=q;
	  IPFmapBD[3*j+1]=V;
	  IPFmapBD[3*j+2]=p;
	}
	if (icase==2){
	  IPFmapBD[3*j]=p;
	  IPFmapBD[3*j+1]=V;
	  IPFmapBD[3*j+2]=t;
	}
	if (icase==3){
	  IPFmapBD[3*j]=p;
	  IPFmapBD[3*j+1]=q;
	  IPFmapBD[3*j+2]=V;
	}
	if (icase==4){
	  IPFmapBD[3*j]=t;
	  IPFmapBD[3*j+1]=p;
	  IPFmapBD[3*j+2]=V;
	}
	if (icase==5){
	  IPFmapBD[3*j]=V;
	  IPFmapBD[3*j+1]=p;
	  IPFmapBD[3*j+2]=q;
	}
      }
      // x dir
      vCD[0] = std::fabs(rRot[0][0]*vX[0]+rRot[1][0]*vX[1]+rRot[2][0]*vX[2]);
      vCD[1] = std::fabs(rRot[0][1]*vX[0]+rRot[1][1]*vX[1]+rRot[2][1]*vX[2]);
      vCD[2] = std::fabs(rRot[0][2]*vX[0]+rRot[1][2]*vX[1]+rRot[2][2]*vX[2]);
      std::sort(vCD,vCD+3);
      std::swap(vCD[0],vCD[1]);
      vCD[2]=std::min(vCD[2],1.0);
      xp = 2.*vCD[0]/(1.+vCD[2]);
      yp = 2.*vCD[1]/(1.+vCD[2]);
      H=atan( (yp-y0)/(xp-x0))*180./M_PI;
      xp < x0 ? H+=180: H;
      H = H+240-atan((triPts[1][2]-y0)/(triPts[0][2]-x0))*180/M_PI;
      V=1.;
      S=pow(pow(xp-x0,2.)+pow(yp-y0,2.),.5);
      S=S/sMax*0.8 + 0.2;
      H>=360.0 ? H=0.0 : H;
      icase = floor(H/60.0);
      ff= H/6.0 - icase;
      p=V*(1.-S);
      q=V*(1.-S*ff);
      t=V*(1.-(S*(1.-ff)));
      if (S<=0.0){
	IPFmapx[3*j]=0.0;
	IPFmapx[3*j+1]=0.0;
	IPFmapx[3*j+2]=0.0;
      } else {
	if (icase==0){
	  IPFmapx[3*j]=V;
	  IPFmapx[3*j+1]=t;
	  IPFmapx[3*j+2]=p;
	}
	if (icase==1){
	  IPFmapx[3*j]=q;
	  IPFmapx[3*j+1]=V;
	  IPFmapx[3*j+2]=p;
	}
	if (icase==2){
	  IPFmapx[3*j]=p;
	  IPFmapx[3*j+1]=V;
	  IPFmapx[3*j+2]=t;
	}
	if (icase==3){
	  IPFmapx[3*j]=p;
	  IPFmapx[3*j+1]=q;
	  IPFmapx[3*j+2]=V;
	}
	if (icase==4){
	  IPFmapx[3*j]=t;
	  IPFmapx[3*j+1]=p;
	  IPFmapx[3*j+2]=V;
	}
	if (icase==5){
	  IPFmapx[3*j]=V;
	  IPFmapx[3*j+1]=p;
	  IPFmapx[3*j+2]=q;
	}
      }
      // y dir 
      vCD[0] = std::fabs(rRot[0][0]*vY[0]+rRot[1][0]*vY[1]+rRot[2][0]*vY[2]);
      vCD[1] = std::fabs(rRot[0][1]*vY[0]+rRot[1][1]*vY[1]+rRot[2][1]*vY[2]);
      vCD[2] = std::fabs(rRot[0][2]*vY[0]+rRot[1][2]*vY[1]+rRot[2][2]*vY[2]);
      std::sort(vCD,vCD+3);
      std::swap(vCD[0],vCD[1]);
      vCD[2]=std::min(vCD[2],1.0);
      xp = 2.*vCD[0]/(1.+vCD[2]);
      yp = 2.*vCD[1]/(1.+vCD[2]);
      H = H+240-atan((triPts[1][2]-y0)/(triPts[0][2]-x0))*180/M_PI;
      V=1.;
      S=pow(pow(xp-x0,2.)+pow(yp-y0,2.),.5);
      S=S/sMax*0.8 + 0.2;
      H>=360.0 ? H=0.0 : H;
      icase = floor(H/60.0);
      ff= H/6.0 - icase;
      p=V*(1.-S);
      q=V*(1.-S*ff);
      t=V*(1.-(S*(1.-ff)));
      if (S<=0.0){
	IPFmapy[3*j]=0.0;
	IPFmapy[3*j+1]=0.0;
	IPFmapy[3*j+2]=0.0;
      } else {
	if (icase==0){
	  IPFmapy[3*j]=V;
	  IPFmapy[3*j+1]=t;
	  IPFmapy[3*j+2]=p;
	}
	if (icase==1){
	  IPFmapy[3*j]=q;
	  IPFmapy[3*j+1]=V;
	  IPFmapy[3*j+2]=p;
	}
	if (icase==2){
	  IPFmapy[3*j]=p;
	  IPFmapy[3*j+1]=V;
	  IPFmapy[3*j+2]=t;
	}
	if (icase==3){
	  IPFmapy[3*j]=p;
	  IPFmapy[3*j+1]=q;
	  IPFmapy[3*j+2]=V;
	}
	if (icase==4){
	  IPFmapy[3*j]=t;
	  IPFmapy[3*j+1]=p;
	  IPFmapy[3*j+2]=V;
	}
	if (icase==5){
	  IPFmapy[3*j]=V;
	  IPFmapy[3*j+1]=p;
	  IPFmapy[3*j+2]=q;
	}
      } // if (S<0.0...
    } //     if (gID[j]<1){
  } // for (int j..
  for (int j=0;j<4*nGrain;++j){cth[j]=cTheta[j];}
  unsigned int nVoxT, js,jc,ncth;
  nVoxT = g.nX[0]*g.nX[1]*g.nX[2];
  ncth=4*nGrain;
  js = 0;
  jc = nVoxT;
  std::vector<float> dX(g.dX,g.dX+3);
  adios2::ADIOS adios;
  adios2::IO hdf5IO = adios.DeclareIO("HDFFileIO");
  hdf5IO.SetEngine("HDF5");
  // global array : name, { shape (total) }, { start (local) }, { count (local) }  all are constant dimensions
  adios2::Variable<int> dimsa = hdf5IO.DefineVariable<int>(
	      "dims", {3}, {0}, {3});
  adios2::Variable<float> dxa = hdf5IO.DefineVariable<float>(
	      "VoxelDX", {3}, {0}, {3});
  adios2::Variable<int> gida = hdf5IO.DefineVariable<int>(
	      "gID", {nVoxT}, {js}, {jc});
  adios2::Variable<int> vStatea = hdf5IO.DefineVariable<int>(
	      "vState", {nVoxT}, {js}, {jc});
  adios2::Variable<float> TempOuta = hdf5IO.DefineVariable<float>(
	      "Temperature", {nVoxT}, {js}, {jc});
  adios2::Variable<float> IPFmapBDa = hdf5IO.DefineVariable<float>(
	      "IPFz", {3*nVoxT}, {3*js}, {3*jc});
  adios2::Variable<float> IPFmapxa = hdf5IO.DefineVariable<float>(
	      "IPFx", {3*nVoxT}, {3*js}, {3*jc});
  adios2::Variable<float> IPFmapya = hdf5IO.DefineVariable<float>(
	      "IPFy", {3*nVoxT}, {3*js}, {3*jc});
  adios2::Variable<float> angAx = hdf5IO.DefineVariable<float>(
	      "angleAxis", {ncth}, {0}, {ncth});
  adios2::Engine hdf5Writer =
      hdf5IO.Open(hdf5Filename, adios2::Mode::Write);
  hdf5Writer.Put<int>(dimsa, g.nX);
  hdf5Writer.Put<float>(dxa, dX.data());
  hdf5Writer.Put<int>(gida, gID);
  hdf5Writer.Put<int>(vStatea, vState);
  hdf5Writer.Put<float>(TempOuta, TempOut.data());
  hdf5Writer.Put<float>(IPFmapBDa, IPFmapBD.data());
  hdf5Writer.Put<float>(IPFmapxa, IPFmapx.data());
  hdf5Writer.Put<float>(IPFmapya, IPFmapy.data());
  hdf5Writer.Put<float>(angAx, cth.data());
  hdf5Writer.Close();
} // end WriteToHDF1
