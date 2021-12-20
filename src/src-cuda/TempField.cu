// member functions for TempField

#include "Grid.cuh"
#include "TempField.cuh"
#include "fstream"
// #include "iostream"
#include "math.h"
#include "numeric"
#include <algorithm>

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
