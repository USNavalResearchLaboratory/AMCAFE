// member functions for TempField

#include "Grid.h"
#include "TempField.h"
#include "fstream"
// #include "iostream"
#include "math.h"
#include "Partition.h"
#include "BasePlate.h"
#include "numeric"
#include <algorithm>
#include "mpi.h"

// constructor
TempField::TempField(Grid &g, Partition & part, BasePlate &bp)
{
  _xyz = &g;
  _part = &part;
  _bp = &bp;
  TempCurr.resize(_part->ncellLoc+_part->nGhost,0.0);
  DDtTemp.assign(_part->ncellLoc+_part->nGhost,0.0);
  ilaserLoc = _bp->Nzh + 1; // make it one higher than base plate
  tInd = 0;
  //ilaserLoc = _xyz->nX.back(); // make z (3rd) direction the build direction   

} // end TempField

void TempField::InitializeMoose(std::string &filnambase,
		     const int &NtIn, const double &dtMIn, const std::vector<int> &nXMIn,
		     const std::vector<double> &dXMIn)
{
  filnambaseptr = &filnambase;
  NtM = NtIn;
  dtM = dtMIn;
  nXM = nXMIn;
  dXM = dXMIn;
  Temp.resize(2,std::vector<double>(_part->ncellLoc + _part->nGhost,0.0));
  indexM = 0;
} // end InitializeMoose()
void TempField::Test2ComputeTemp(double T20, double T10, double a,double tcurr)
{
  int Ntot = _part->nGhost+_part->ncellLoc, j1,j2,j3,j;
  double T1,T2,x,y,z,Lx[3];
  for (int j=0;j<3;++j){Lx[j] = _xyz->dX[j]*_xyz->nX[j];}
  T2 = T20 - a*tcurr;
  T1 = T10 - a*tcurr;
  TempCurr.assign(Ntot,0.0);
  for (int j=0;j<Ntot;++j){
    j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
    j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
    j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
    x = (j1+.5)*_xyz->dX[0];
    y = (j2+.5)*_xyz->dX[1];
    z = (j3+.5)*_xyz->dX[2];
    //TempCurr[j]  = (T2-T1)/3*(x/Lx[0]+y/Lx[1]+z/Lx[2]) + T1;
    TempCurr[j]  = (T2-T1)*(z/Lx[2]) + T1;
    DelT = .000006;
    DelT = .0006;
    //DelT = .001;
    //DelT = .0125;
  } // for (int j...     
}
void TempField::InitializeAnalytic()
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
  bmV = _xyz->bmV;
  patternID = _xyz->patternID;
  a1.resize(6);
  Qp.resize(2);
  a1[0] = _xyz->meltparam[0];
  a1[1] = _xyz->meltparam[2];
  a1[2] = _xyz->meltparam[3];
  a1[3] = _xyz->meltparam[1];
  a1[4] = a1[1];
  a1[5] = a1[2];
  T0 = _xyz->T0;
  zlaserOff=1.0; // 1.0 (this specifies where laser z value is - see SchwalbachTempCurr)
  DelT = _xyz->bmDelT;
  bmDX = {DelT*bmV,_xyz->bhatch,_xyz->layerT}; // (SD, TD, BD) per layer
  offset=_xyz->offset;
  shiftL={a1[3],0.0,0.0}; // (SD, TD, BD) per layer
  bmLx={_xyz->LX[0]+shiftL[0],_xyz->LX[1],_xyz->LX[2]}; // (SD, TD, BD) per layer
  nTTemp = {_xyz->nTsd==2 ? _xyz->nTsd : int(ceil(bmLx[0]/bmDX[0] ))+1,
	    bmDX[1]>_xyz->LX[1] ? 1 : int(ceil(bmLx[1]/bmDX[1]))+1,
            bmDX[2]<std::numeric_limits<double>::epsilon() ? 1: int(ceil(bmLx[2]/bmDX[2]))};
  bmLx={(nTTemp[0]-1)*bmDX[0],(nTTemp[1]-1)*bmDX[1],(nTTemp[2]-1)*bmDX[2]};
  
} // end InitializeAnalytic 

void TempField::AnalyticTempCurr(double tcurr,std::vector<double> & TempOut, std::vector<int> &icellid, int Ntot)
{
  //computes temp field based on Schwalbach et al
  int j1,j2,j3,iplay;
  double x0,y0,z0,x,y,z,dsq,dsq2,bx,by,xi;
  std::vector<double> rij1(3),xs1(3),xs2(3);
  ilaserLoc = _bp->Nzh + (floor( round(tcurr/DelT)/(nTTemp[0]*nTTemp[1]))+1)*_xyz->nZlayer;
  iplay=_xyz->nX[0]*_xyz->nX[1]*ilaserLoc;
  TempOut.assign(Ntot,T0);
  xi = _xyz->tL*(1+std::numeric_limits<double>::epsilon() );
  if (patternID==1){
    int js1, js2;
    std::vector<double> a1m(6);
    for (int j=0;j<6;++j){a1m[j]=a1[j];}
    // x,y,z spatial location of source
    y = floor(fmod(round(tcurr/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
    z = (floor((tcurr/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height -
      ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
    js1 = fmod( round(tcurr/DelT),nTTemp[0]*nTTemp[1]);
    js2 = fmod(floor(js1/nTTemp[0])+1,2);
    x = (.5*pow(-1,js2)+.5)*bmLx[0] -
      pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] -
      (pow(-1,js2)+1)/2.0*shiftL[0] + pow(-1,js2)*offset[0];
    for (int j=0;j<Ntot;++j){
      j3 = floor(icellid[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (icellid[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = icellid[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      if (z0>z){continue;}
      rij1[0] = x0-x;
      rij1[1] = y0-y;
      rij1[2] = z0-z;
      if (rij1[0]*pow(-1,(js2+1))  >=0){
	dsq = pow(rij1[1]/a1m[1],2.0)+pow(rij1[2]/a1m[2],2.0);
	if (dsq<1.0 && (fabs(rij1[0])<bmDX[0]) ){
	TempOut[j]=xi;
	} else {
	  TempOut[j] = _xyz->tS;
	}	
      } else {
	dsq = pow(rij1[0]/a1m[3],2.0)+pow(rij1[1]/a1m[4],2.0)+pow(rij1[2]/a1m[5],2.0);
	if (dsq<1.0){
	  TempOut[j]=xi;
	} else {
	  TempOut[j] = _xyz->tS;
	}
      } // if (rij1[0]>0
    } // for (int j...
  } // if (patternID==1 ...
  if (patternID==4){
    int js1, js2,js3;
    std::vector<double> xy(2);
    std::vector<int> jxy(2),ja(6);
    // x,y,z spatial location of source
    z = (floor((tcurr/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height -
      ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
    js1 = fmod( round(tcurr/DelT),nTTemp[0]*nTTemp[1]);
    js2 = fmod(floor(js1/nTTemp[0])+1,2);
    js3 = fmod( floor( round(tcurr/DelT)/(nTTemp[0]*nTTemp[1])),2);
    if (js3==0){
      jxy={0,1};
      ja={0,1,2,3,4,5};
    } else {
      jxy={1,0};
      ja={1,0,2,4,3,5};
    }
    xy[0] = (.5*pow(-1,js2)+.5)*bmLx[0] -
	pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] -
	(pow(-1,js2)+1)/2.0*shiftL[0] + pow(-1,js2)*offset[0];
    xy[1]=floor(fmod(round(tcurr/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
    for (int j=0;j<Ntot;++j){
      j3 = floor(icellid[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (icellid[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = icellid[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      if (z0>z){continue;}
      rij1[jxy[0]] = x0-xy[jxy[0]];
      rij1[jxy[1]] = y0-xy[jxy[1]];
      rij1[2] = z0-z;
      if (rij1[0]*pow(-1,(js2+1))  >=0){
	dsq = pow(rij1[1]/a1[1],2.0)+pow(rij1[2]/a1[2],2.0);
	if (dsq<1.0 && (fabs(rij1[0])<bmDX[0]) ){
	TempOut[j]=xi;
	} else {
	  TempOut[j] = _xyz->tS;
	}	
      } else {
	dsq = pow(rij1[0]/a1[3],2.0)+pow(rij1[1]/a1[4],2.0)+pow(rij1[2]/a1[5],2.0);
	if (dsq<1.0){
	  TempOut[j]=xi;
	} else {
	  TempOut[j] = _xyz->tS;
	}
      } // if (rij1[0]>0
    } // for (int j...
  } // if (patternID==4 ...
} // end AnalyticTempCurr()            

void TempField::ComputeStartTime()
{
  /*
    determines the last time a voxel's temperature goes below the solidification temperature. This
    pertains to voxels that will get melted in time period [tcurr, tcurr+Delta t]. The start time,
    array tBeg, is used in the UpdateVoxels to determine time until. This should allow for much
    larger time steps. This time is determined by when the voxel center point is on the isosurface
    of the melt pool
   */
  int j1,j2,j3,iplay,js0,js1,js2,Ntot= _part->ncellLoc+_part->nGhost;
  double x0,x1,ds1,tcurr,xs1,tcheck1;
  std::vector<int> iloc(3);
  std::vector<double> xs(3),rij1(3),xs2(3),sigm(3),xlt(6);
  tBeg0 = 1e6;
  tcheck1=1.5*tBeg0;
  tcurr=_xyz->time;
  ilaserLoc = _bp->Nzh + (floor( round(tcurr/DelT)/(nTTemp[0]*nTTemp[1]))+1)*_xyz->nZlayer;
  iplay=_xyz->nX[0]*_xyz->nX[1]*ilaserLoc;
  sigm[0] = a1[3];
  sigm[1] = a1[4];
  sigm[2] = a1[5];
  if (patternID==1){
    // xlt spatial location of source at t and t+DelT
    xlt[1] = floor(fmod(round(tcurr/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
    xlt[2] = (floor((tcurr/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height -
      ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
    js1 = fmod( round(tcurr/DelT),nTTemp[0]*nTTemp[1]);
    js2 = fmod(floor(js1/nTTemp[0])+1,2); // js2=1: laser going +x; js2=0; laser going -x;
    xlt[0] = (.5*pow(-1,js2)+.5)*bmLx[0] -
      pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] -
      (pow(-1,js2)+1)/2.0*shiftL[0] + pow(-1,js2)*offset[0];
    xlt[4] = xlt[1];
    xlt[5] = xlt[2];
    xlt[3] = xlt[0] -  pow(-1,js2)*bmDX[0];
    iloc={0,1,2};
  }
  if (patternID==2){
    // xlt spatial location of source at t and t+DelT
    xlt[2] = (floor((tcurr/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height -
      ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
    xlt[5] = xlt[2];
    js1 = fmod( round(tcurr/DelT),nTTemp[0]*nTTemp[1]);
    js0 = fmod( floor( round(tcurr/DelT)/(nTTemp[0]*nTTemp[1])),2);
    if (js0==0){
      xlt[1] = floor(fmod(round(tcurr/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];;
      js2 = fmod(floor(js1/nTTemp[0])+1,2);
      xlt[0] = (.5*pow(-1,js2)+.5)*bmLx[0] -
        pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] -
        (pow(-1,js2)+1)/2.0*shiftL[0] + pow(-1,js2)*offset[0];
      xlt[4] = xlt[1];
      xlt[3] = xlt[0] - pow(-1,js2)*bmDX[0];
      iloc={0,1,2};
    } else {
      xlt[0] = floor(fmod(round(tcurr/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];;
      js2 = fmod(floor(js1/nTTemp[0])+1,2);
      xlt[1] = (.5*pow(-1,js2)+.5)*bmLx[0] -
        pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] -
        (pow(-1,js2)+1)/2.0*shiftL[0] + pow(-1,js2)*offset[0];
      xlt[3] = xlt[0];
      xlt[4] = xlt[1] - pow(-1,js2)*bmDX[0];
      iloc={1,0,2};
      sigm[0] = a1[4];
      sigm[1] = a1[3];
    }
  }
  tBeg.assign(Ntot,tBeg0);
  x0=std::min(xlt[iloc[0]],xlt[iloc[0]+3]);
  x1=std::max(xlt[iloc[0]],xlt[iloc[0]+3]);
  ds1 = (xlt[iloc[0]+3] - xlt[iloc[0]])/std::fabs(xlt[iloc[0]+3] - xlt[iloc[0]]);
  for (int j=0;j<Ntot;++j){
    if (_part->icellidLoc[j] >= iplay){continue;}
    j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
    j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
    j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
    xs[0]=(double(j1)+.5)*_xyz->dX[0];
    xs[1]=(double(j2)+.5)*_xyz->dX[1];
    xs[2]=(double(j3)+.5)*_xyz->dX[2];
    rij1[1] = std::fabs(xs[iloc[1]] - xlt[iloc[1]+3]);
    rij1[2] = std::fabs(xs[iloc[2]] - xlt[iloc[2]+3]);
    if (rij1[1] < sigm[1] && rij1[2] < sigm[2]){
      rij1[0] = sigm[0]*pow(1 - pow(rij1[1]/sigm[1],2.0)- pow(rij1[2]/sigm[2],2.0),.5);
      xs1 = xs[iloc[0]] + ds1*rij1[0];
      if (xs1 >=x0 &&  xs1<=x1){
        tBeg[j] = ds1*(xs1 - xlt[iloc[0]])/bmV;
      }
      if ( (js2==1 && xs1>x1)){tBeg[j] = tcheck1;}
      if ( (js2==0 && xs1<x0)){tBeg[j] = tcheck1;}
    }
  } // for (int j=0...
} // end StartTime                             

void TempField::InitializeSchwalbach()
{
  /*
    This follows the model: Schwalbach, Edwin J., et al. "A discrete source model of powder 
    bed fusion additive manufacturing thermal history." Additive Manufacturing 25 (2019): 485-498.
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
  bmSTD = _xyz->beamSTD;
  bmV = _xyz->bmV;
  T0targ = _xyz->T0targ;
  bmEta = _xyz->beamEta;
  patternID = _xyz->patternID;
  T0 = _xyz->T0;
  zlaserOff=1.0; // 1.0 (this specifies where laser z value is - see SchwalbachTempCurr)
  double tb,minTemp;
  if (patternID==0 || patternID==1 || patternID==2 || patternID==3 || patternID==4){
    DelT = 4.0/3.0*bmSTD[0]/bmV;
    //bmDX = {DelT*bmV,4.0/3.0*bmSTD[1],_xyz->layerT};
    bmDX = {DelT*bmV,2.7*bmSTD[1],_xyz->layerT};
    double x1 = pow( pow(bmSTD[0],2.0)*pow(bmSTD[1],2.0)*pow(bmSTD[2],2.0),.5);
    bmP = T0targ*_xyz->cP*_xyz->rho*pow(2.0,.5)*pow(M_PI,1.5)*x1/DelT;
    rcut = pow( -2* pow(*std::min_element(bmSTD.begin(),bmSTD.end()),2.0)*log(.001),.5);
    Ci = bmEta*bmP*DelT/(_xyz->rho*_xyz->cP*pow(2.0,.5)*pow(M_PI,1.5));
    alpha = _xyz->kappa/_xyz->cP/_xyz->rho;
    minTemp = 0.250; // cut off change in temperature for tcut
    tcut = pow(Ci/(2*alpha*minTemp),2.0/3.0);		       
    nSource = ceil(tcut/DelT);
    offset=_xyz->offset;
    bmLx={_xyz->LX[0]+1*bmDX[0],_xyz->LX[1]+1*bmDX[1],_xyz->LX[2]};
    if (patternID==1){
      offset={0.0,-_xyz->nX[1]*_xyz->dX[1]/2.0,0.0}; // positive value means starting outside domain
      shiftL={3*bmDX[0],0.0,0.0};
      bmLx={_xyz->LX[0]+shiftL[0],_xyz->LX[1]+1*bmDX[1],_xyz->LX[2]};
    } // if (patternID==1...
    nTTemp = {int(floor(bmLx[0]/bmDX[0] ))+1,int(floor(bmLx[1]/bmDX[1]))+1,
            int(floor(bmLx[2]/bmDX[2]))};
    bmLx={(nTTemp[0]-1)*bmDX[0],(nTTemp[1]-1)*bmDX[1],(nTTemp[2]-1)*bmDX[2]};
  }

} // end InitializeSchwalbach

void TempField::SchwalbachTempCurr()
{
  //computes temp field based on Schwalbach et al
  int Ntot = _part->nGhost+_part->ncellLoc, j1,j2,j3,iplay;
  double x0,y0,z0,rij,x,y,z,tc,xst;
  std::vector<double> lam(3);
  ilaserLoc = _bp->Nzh + floor( (floor( round(_xyz->time/DelT)/(nTTemp[0]*nTTemp[1]))+1)*_xyz->layerT/_xyz->dX[2]);
  iplay=_xyz->nX[0]*_xyz->nX[1]*ilaserLoc;
  TempCurr.assign(Ntot,T0);
  if (patternID==0){
    for (int j=0;j<Ntot;++j){
      //if (_part->icellidLoc[j] >= iplay){continue;}
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = _xyz->time - jt*DelT;
	x = fmod(round(tc/DelT),(nTTemp[0]))*bmDX[0] - offset[0];
	y = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	z = (floor((tc/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
	  ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
	lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempCurr[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==0 ...
  if (patternID==1){
    int js1, js2;
    y = floor(fmod(round(_xyz->time/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
    z = (floor((_xyz->time/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
      ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
    for (int j=0;j<Ntot;++j){
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = jt*DelT;
	js1 = fmod( round(_xyz->time/DelT),nTTemp[0]*nTTemp[1]);
	js2 = fmod(floor(js1/nTTemp[0])+1,2);
	x = (.5*pow(-1,js2)+.5)*bmLx[0] - 
	  pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] + pow(-1,js2)*tc*bmV - 
	  (pow(-1,js2)+1)/2.0*shiftL[0] + pow(-1,js2)*offset[0] ;
	//if (_part->myid==0){if (j==0 & jt==0){std::cout<< tInd<<","<<x<<","<<y<<","<<z<<std::endl;}}
	lam = {pow(bmSTD[0],2.0)+2*alpha*(tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(tc),
	       pow(bmSTD[2],2.0)+2*alpha*(tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempCurr[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==1 ...
  if (patternID==2){
    int js1, js2,js3;
    for (int j=0;j<Ntot;++j){
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = _xyz->time - jt*DelT;
	js1 = fmod( round(tc/DelT),nTTemp[0]*nTTemp[1]);
	js3 = fmod( floor( round(tc/DelT)/(nTTemp[0]*nTTemp[1]) ),2);
	js2 = fmod(floor(js1/nTTemp[0])+1+js3,2);
	x = (.5*pow(-1,js2)+.5)*bmLx[0] - 
	  pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] - offset[0];
	y = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	z = (floor((tc/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
	  ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
	//if (_part->myid==0){if (j==0 & jt==0){std::cout<< tInd<<","<<x<<","<<y<<","<<z<<std::endl;}}
	lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempCurr[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==2 ...
  if (patternID==3){
    int js;
    for (int j=0;j<Ntot;++j){
      //if (_part->icellidLoc[j] >= iplay){continue;}
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = _xyz->time - jt*DelT;
	js = fmod( floor( round(tc/DelT)/(nTTemp[0]*nTTemp[1])),2);
	if (js==0){
	  x = fmod(round(tc/DelT),(nTTemp[0]))*bmDX[0] - offset[0];
	  y = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	} else {
	  y = fmod(round(tc/DelT),(nTTemp[0]))*bmDX[0] - offset[0];
	  x = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	}
	z = (floor((tc/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
	  ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
	lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempCurr[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==3 ...
  if (patternID==4){
    int js1, js2,js3;
    for (int j=0;j<Ntot;++j){
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = _xyz->time - jt*DelT;
	js1 = fmod( round(tc/DelT),nTTemp[0]*nTTemp[1]);
	js2 = fmod(floor(js1/nTTemp[0])+1,2);
	js3 = fmod( floor( round(tc/DelT)/(nTTemp[0]*nTTemp[1])),2);
	if (js3==0){
	  x = (.5*pow(-1,js2)+.5)*bmLx[0] - 
	    pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] - offset[0];
	  y = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	} else {
	  y = (.5*pow(-1,js2)+.5)*bmLx[0] - 
	    pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] - offset[0];
	  x = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	}
	z = (floor((tc/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
	  ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
	//if (_part->myid==0){if (j==0 & jt==0){std::cout<< tInd<<","<<x<<","<<y<<","<<z<<std::endl;}}
	lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempCurr[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==4 ...
} // end SchwalbachTempCurr()

void TempField::SchwalbachTempCurr(double tcurr,std::vector<double> & TempOut )
{
  //computes temp field based on Schwalbach et al
  int Ntot = _part->nGhost+_part->ncellLoc, j1,j2,j3,iplay;
  double x0,y0,z0,rij,x,y,z,tc;
  std::vector<double> lam(3);
  if (patternID==0){
    for (int j=0;j<Ntot;++j){
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = tcurr  - jt*DelT;
	x = fmod((tc/DelT),(nTTemp[0]))*bmDX[0] - offset[0];
	y = floor(fmod((tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	z = (floor((tc/DelT)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
	  ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
	lam = {pow(bmSTD[0],2.0)+2*alpha*(tcurr - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(tcurr - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(tcurr - tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempOut[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==0 ...
  if (patternID==1){
    int js1, js2;
    y = floor(fmod(round(_xyz->time/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
    z = (floor((_xyz->time/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
      ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
    for (int j=0;j<Ntot;++j){
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = jt*DelT;
	js1 = fmod( round(_xyz->time/DelT),nTTemp[0]*nTTemp[1]);
	js2 = fmod(floor(js1/nTTemp[0])+1,2);
	x = (.5*pow(-1,js2)+.5)*bmLx[0] - 
	  pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] + pow(-1,js2)*tc*bmV - 
	  (pow(-1,js2)+1)/2.0*shiftL[0] + pow(-1,js2)*offset[0] ;
	lam = {pow(bmSTD[0],2.0)+2*alpha*(tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(tc),
	       pow(bmSTD[2],2.0)+2*alpha*(tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempOut[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==1 ...
  if (patternID==2){
    int js1, js2,js3;
    for (int j=0;j<Ntot;++j){
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = _xyz->time - jt*DelT;
	js1 = fmod( round(tc/DelT),nTTemp[0]*nTTemp[1]);
	js3 = fmod( floor( round(tc/DelT)/(nTTemp[0]*nTTemp[1]) ),2);
	js2 = fmod(floor(js1/nTTemp[0])+1+js3,2);
	x = (.5*pow(-1,js2)+.5)*bmLx[0] - 
	  pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] - offset[0];
	y = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	z = (floor((tc/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
	  ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
	//if (_part->myid==0){if (j==0 & jt==0){std::cout<< tInd<<","<<x<<","<<y<<","<<z<<std::endl;}}
	lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempOut[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==2 ...
  if (patternID==3){
    int js;
    for (int j=0;j<Ntot;++j){
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = _xyz->time - jt*DelT;
	js = fmod( floor( round(tc/DelT)/(nTTemp[0]*nTTemp[1])),2);
	if (js==0){
	  x = fmod(round(tc/DelT),(nTTemp[0]))*bmDX[0] - offset[0];
	  y = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	} else {
	  y = fmod(round(tc/DelT),(nTTemp[0]))*bmDX[0] - offset[0];
	  x = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	}
	z = (floor((tc/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
	  ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
	lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempOut[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==3 ...
  if (patternID==4){
    int js1, js2,js3;
    for (int j=0;j<Ntot;++j){
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      for (int jt=0;jt<nSource;++jt){
	// tc,x,y,z space-time location of source
	tc = _xyz->time - jt*DelT;
	js1 = fmod( round(tc/DelT),nTTemp[0]*nTTemp[1]);
	js2 = fmod(floor(js1/nTTemp[0])+1,2);
	js3 = fmod( floor( round(tc/DelT)/(nTTemp[0]*nTTemp[1])),2);
	if (js3==0){
	  x = (.5*pow(-1,js2)+.5)*bmLx[0] - 
	    pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] - offset[0];
	  y = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	} else {
	  y = (.5*pow(-1,js2)+.5)*bmLx[0] - 
	    pow(-1,js2)*fmod(js1,nTTemp[0])*bmDX[0] - offset[0];
	  x = floor(fmod(round(tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
	}
	z = (floor((tc/DelT+DelT*1e-6)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
	  ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];
	lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempOut[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==4 ...
} // end SchwalbachTempCurr()

void TempField::SchwalbachDDtTemp()
{
  //computes temp field time derivative based on Schwalbach et al
  int Ntot = _part->ncellLoc+_part->nGhost,j1,j2,j3,iplay;
  double x0,y0,z0,rij,x,y,z,tc,tmin,rdij;
  std::vector<double> lam(3);
  DDtTemp.assign(Ntot,0.0);
  iplay=_xyz->nX[0]*_xyz->nX[1]*ilaserLoc;
  for (int j=0;j<Ntot;++j){
    if (_part->icellidLoc[j] > iplay){continue;}
    j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
    j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
    j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
    x0 = (double(j1)+.5)*(_xyz->dX[0]);
    y0 = (double(j2)+.5)*(_xyz->dX[1]);
    z0 = (double(j3)+.5)*(_xyz->dX[2]);
    for (int jt=0;jt<nSource;++jt){
      // tc,x,y,z space-time location of source
      tc = tmin + jt*DelT;
      x = fmod((tc/DelT),(nTTemp[0]))*bmDX[0] - offset[0];
      y = floor(fmod((tc/DelT),(nTTemp[0]*nTTemp[1]))/nTTemp[0])*bmDX[1]-offset[1];
      z = (floor((tc/DelT)/(nTTemp[0]*nTTemp[1]))+zlaserOff)*bmDX[2] + _bp->height - 
	  ceil(offset[2]/_xyz->dX[2])*_xyz->dX[2];

      lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	     pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	     pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
      rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	pow(z-z0,2.0)/2/lam[2];
      rdij = pow((x-x0)/lam[0],2.0) +pow((y-y0)/lam[1],2.0) +
	pow((z-z0)/lam[2],2.0) - 1/lam[0] - 1/lam[1] - 1/lam[2];
      DDtTemp[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij) * 
	alpha*rdij;
    } // for (int j1...
  } // for (int j...
} // end SchwalbachDDtTemp()

void TempField::SchwalbachGradTemp()
{
  //computes temp field gradient based on Schwalbach et al

} // end SchwalbachGradTemp()

void TempField::ReadCSVMoose2()
{
  // reads in time index and finds appropriate moose file to open
  std::vector<double> tempLoc;
  indexM = int(floor(_xyz->time/dtM));
  for (int jtime=0;jtime<2;++jtime){
    if (_part->myid ==0){

      std::string filename,cc;
      std::ifstream file;
      int indM,jx,jy,jn,jz,nM,cc1;
      double tmp,xM,yM,zM;
      nM = std::accumulate(nXM.begin(),nXM.end(),1,std::multiplies<int>());
      std::vector<double> tempM(nM);
      if (jtime==0){
	indM = int(floor(_xyz->time/dtM));
	filename = *filnambaseptr+std::to_string(indM)+".csv";
	file.open(filename);
	while (!file.good()){
	  indM -=1;
	  filename = *filnambaseptr+std::to_string(indM)+".csv";
	  file.open(filename);
	  if (indM < 0){throw std::runtime_error("no start time for temperature found");}
	}
      } else {
	indM+=1;
	filename = *filnambaseptr+std::to_string(indM)+".csv";
	file.open(filename);
	cc1=0;
	while (!file.good()){
	  indM +=1;
	  cc1+=1;
	  filename = *filnambaseptr+std::to_string(indM)+".csv";
	  file.open(filename);
	  if (cc1 > 50){throw std::runtime_error("cant find temperature field for tInd+1");}
	}
      }
      getline(file,cc,'\n');    
      for (int jM=0;jM<nM;++jM){
	getline(file,cc,',');
	tmp = std::stod(cc);
	getline(file,cc,',');
	getline(file,cc,',');
	getline(file,cc,',');
	xM = std::stod(cc);
	getline(file,cc,',');
	yM = std::stod(cc);
	getline(file,cc,'\n');
	zM = std::stod(cc);
	jx = round(xM/dXM[0]);
	jy = round(yM/dXM[1]);
	jz = round(yM/dXM[2]);
	jn = nXM[0]*nXM[1]*jz+nXM[0]*jy+jx;
	tempM[jn] = tmp;
      } // end for jM
      file.close();
      // interpolate to CA grid
      for (int j=0;j<_part->nprocs;++j){
	tempLoc.assign(_part->icellidLocAll[j].size(),0.0);
	InterpolateToGrid2(tempM,_part->icellidLocAll[j],tempLoc);
	if (j==0){
	  for (int j1=0;j1<(_part->nGhost+_part->ncellLoc);++j1){
	    Temp[jtime][j1] = tempLoc[j1];
	  }
	} else {
	  MPI_Send(&tempLoc[0],_part->icellidLocAll[j].size(),MPI_DOUBLE,j,0,MPI_COMM_WORLD);
	}
      } // for (int = j ...
    } // if (_part->myid==0)
    if (_part->myid > 0){
      tempLoc.assign(_part->nGhost+_part->ncellLoc,0);
      MPI_Recv(&tempLoc[0],(_part->nGhost+_part->ncellLoc),MPI_DOUBLE,0,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
      for (int j1=0;j1<(_part->nGhost+_part->ncellLoc);++j1){
	Temp[jtime][j1] = tempLoc[j1];
      }
    } // if (_part->myid >0)
  } // for (int jtime=0...
  double t1,t2;
  t1 = floor(_xyz->time/dtM)*dtM;
  t2 = (floor(_xyz->time/dtM) + 1 )*dtM;
  for (int j=0;j< (_part->ncellLoc + _part->nGhost);++j){
    TempCurr[j] = (Temp[1][j]-Temp[0][j])/(t2-t1)*(_xyz->time-t1) + Temp[0][j];
  }
  MPI_Barrier(MPI_COMM_WORLD);
} // end ReadVTKMoose2

void TempField::Test2(){
  // sets the temperature field to 1584
  int NxyM = std::accumulate(nXM.begin(),nXM.end(),1,std::multiplies<int>());
  int Ntot = _part->nGhost+_part->ncellLoc;
  std::vector<double> tempM(NxyM,0),tempLoc(Ntot,0);
  for (int jM=0;jM<NxyM;++jM){
    tempM[jM] = 1584;
  } // end for jM
  InterpolateToGrid2(tempM,_part->icellidLoc,tempLoc);
  for (int j1=0;j1<(Ntot);++j1){
    Temp[0][j1] = tempLoc[j1];
    Temp[1][j1] = tempLoc[j1];
    TempCurr[j1] = tempLoc[j1];
  }
} // end test2

void TempField::InterpolateToGrid2(const std::vector<double> &tempM,
				   const std::vector<int> &icellid, std::vector<double> &tloc)
{
  // order temperature since paraview order is strange
  double x,y,z,xi1,xi2,xi3;
  int j1,j2,j3,jxM0,jxM1,jyM0,jyM1,jzM0,jzM1,ns1,ns2,nsM1,nsM2;
  std::vector<double> tnode(_xyz->nnodePerCell),Nshape(_xyz->nnodePerCell);
  ns1 = _xyz->nX[0]*_xyz->nX[1];
  ns2 = _xyz->nX[0];
  nsM1 = nXM[0]*nXM[1];
  nsM2 = nXM[0];
  for (int j=0; j < icellid.size(); ++j){
    j3 = floor(icellid[j] /ns1 );
    j2 = floor( (icellid[j] - ns1*j3)/ns2);
    j1 = icellid[j] - ns1*j3 - ns2*j2;
    x = j1*(_xyz->dX[0]);
    y = j2*(_xyz->dX[1]);
    z = j3*(_xyz->dX[2]);
    jxM0 = floor(x/dXM[0]);
    jxM1 = jxM0+1;
    jyM0 = floor(y/dXM[1]);
    jyM1 = jyM0+1; 
    jzM0 = floor(z/dXM[2]);
    jzM1 = jzM0+1; 
    tnode[0] = tempM[nsM1*jzM0 + nsM2*jyM0+jxM0];
    tnode[1] = tempM[nsM1*jzM0 + nsM2*jyM0+jxM1];
    tnode[2] = tempM[nsM1*jzM0 + nsM2*jyM1+jxM1];
    tnode[3] = tempM[nsM1*jzM0 + nsM2*jyM1+jxM0];
    tnode[4] = tempM[nsM1*jzM1 + nsM2*jyM0+jxM0];
    tnode[5] = tempM[nsM1*jzM1 + nsM2*jyM0+jxM1];
    tnode[6] = tempM[nsM1*jzM1 + nsM2*jyM1+jxM1];
    tnode[7] = tempM[nsM1*jzM1 + nsM2*jyM1+jxM0];
    xi1 = 1/dXM[0]*(2*x - (jxM0+jxM1)*dXM[0]);
    xi2 = 1/dXM[1]*(2*y - (jyM0+jyM1)*dXM[1]);
    xi3 = 1/dXM[2]*(2*z - (jzM0+jzM1)*dXM[2]);
    Nshape[0] = .125*(1-xi1)*(1-xi2)*(1-xi3);
    Nshape[1] = .125*(1+xi1)*(1-xi2)*(1-xi3);
    Nshape[2] = .125*(1+xi1)*(1+xi2)*(1-xi3);
    Nshape[3] = .125*(1-xi1)*(1+xi2)*(1-xi3);
    Nshape[4] = .125*(1-xi1)*(1-xi2)*(1+xi3);
    Nshape[5] = .125*(1+xi1)*(1-xi2)*(1+xi3);
    Nshape[6] = .125*(1+xi1)*(1+xi2)*(1+xi3);
    Nshape[7] = .125*(1-xi1)*(1+xi2)*(1+xi3);
    for (int i1=0;i1<_xyz->nnodePerCell;++i1){
	tloc[j] += tnode[i1]*Nshape[i1];
    } // end i1
  } // end for j
} // end InterpolateToGrid2

void TempField::ComputeDDtTempDiscrete()
{
  for (int j=0;j<(_part->ncellLoc+_part->nGhost);++j){
    DDtTemp[j] = (Temp[1][j]-Temp[0][j])/dtM;
  } // for (int j...)
} // end ComputeDDtTemp
