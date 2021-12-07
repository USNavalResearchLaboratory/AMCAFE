// member functions for TempField

#include "Grid.h"
#include "TempField.h"
#include "fstream"
// #include "iostream"
#include "math.h"
#include "Partition.h"
#include "BasePlate.h"
#include "Utilities.h"
#include "numeric"
#include <algorithm>
#include "mpi.h"

// constructor
TempField::TempField(Grid &g, Partition & part, BasePlate &bp, Utilities &ut)
{
  _xyz = &g;
  _part = &part;
  _bp = &bp;
  _ut = &ut;
  TempCurr.resize(_part->ncellLoc+_part->nGhost,0.0);
  DDtTemp.assign(_part->ncellLoc+_part->nGhost,0.0);
  tInd = 0;
  bmV = _xyz->bmV;
  patternID = _xyz->patternID;
  T0 = _xyz->T0;
  zlaserOff=1.0; // 1.0 (this specifies where laser z value is - see SchwalbachTempCurr)
  DelT = _xyz->bmDelT;
  bmDX = {DelT*bmV,_xyz->bhatch,_xyz->layerT}; // (SD, TD, BD) per layer
  offset=_xyz->offset;
  ispvec.assign(_xyz->NpT,0);
  if (patternID==1 || patternID==3){
    std::iota(ispvec.begin(),ispvec.end(),0);
  } // if (patternID==1...
  if (patternID==2 || patternID==4){
    int k;
    for (int j=0;j<_xyz->Ntd;++j){
      k=_xyz->Nsd*j;
      if (fmod(j,2)==0){
	std::iota(ispvec.begin()+k,ispvec.begin()+k+_xyz->Nsd,k);
      } // if (fmod(j,2)==0...
      else {
	for (int j1=0;j1<_xyz->Nsd;++j1){
	  ispvec[k+j1]=k+_xyz->Nsd-1-j1;
	} // for (int j1=0...
      } // else (fmod(j,2)==0...
    } // for (int j=0...
  } // if (patternID==2...
} // end TempField

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
  a1.resize(6);
  Qp.resize(2);
  a1[0] = _xyz->meltparam[0];
  a1[1] = _xyz->meltparam[2];
  a1[2] = _xyz->meltparam[3];
  a1[3] = _xyz->meltparam[1];
  a1[4] = a1[1];
  a1[5] = a1[2];
} // end InitializeAnalytic 

void TempField::EASM(std::vector<double> &TempOut, std::vector<int> &icellid, int Ntot){
  int j1, j2, j3, iplay;
  double x0,y0,x,y,z,dsq,dsq2,bx,by,xi,xp,yp,dirp,zp;
  //  ilaserLoc = _bp->Nzh + (_xyz->indlayer)*_xyz->nZlayer;
  iplay = _xyz->nX[0]*_xyz->nX[1]*_xyz->ilaserLoc;
  TempOut.assign(Ntot,0);
  xi = _xyz->tL*(1+std::numeric_limits<double>::epsilon() );
  int js1,js2;
  x=_xyz->lcoor[2*ispvec[_xyz->isp]]-offset[0];  // Laser x  
  y=_xyz->lcoor[2*ispvec[_xyz->isp]+1]-offset[1];  // Laser y                                                                                                                           
  z=_xyz->ilaserLoc*_xyz->dX[2]-offset[2]; // Laser z                                                                                                                                            
  for (int j=0;j<Ntot;++j){
    j3 = floor(icellid[j]/(_xyz->nX[0]*_xyz->nX[1]));
    j2 = floor( (icellid[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
    j1 = icellid[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
    x0 = (double(j1)+.5)*(_xyz->dX[0]);
    y0 = (double(j2)+.5)*(_xyz->dX[1]);
    zp = (double(j3)+.5)*(_xyz->dX[2]);
    if (zp>z){continue;}
    //    if (zp < 90.0e-6){continue;}
    xp=cos(_xyz->gth)*(x0-_xyz->LX[0]/2.) + sin(_xyz->gth)*(y0-_xyz->LX[1]/2.)+_xyz->LX[0]/2.;
    yp=-sin(_xyz->gth)*(x0-_xyz->LX[0]/2.) + cos(_xyz->gth)*(y0-_xyz->LX[1]/2.)+_xyz->LX[1]/2.; 

    std::vector<double> R = {xp, yp, zp};
    std::vector<double> L = {x, y, z};
    if (x0 > x+0.00020){continue;}

    TempOut[j] = _ut->EASM_Temp_LG(L, R, TempOut[j]);
    }//
  // check if last simulation of scan
  
  double tmelt=_xyz->tL;
  int n1=_part->ncellLoc,icheck,ichecktmp;
  x=_xyz->lcoor2[2*ispvec[_xyz->isp]]-offset[0];
  y=_xyz->lcoor2[2*ispvec[_xyz->isp]+1]-offset[1];
  if (x<_xyz->gbox[0] || x>_xyz->gbox[1] || y<_xyz->gbox[2] || y>_xyz->gbox[3]){
    ichecktmp=std::any_of(TempOut.begin(),TempOut.begin()+n1,[&tmelt]
                          (double tchk){return tchk >= tmelt;});
    MPI_Allreduce(&ichecktmp,&icheck,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (icheck==0 || fmod(_xyz->isp+1,_xyz->Nsd)==0){
      _xyz->inewscanflg=1;
    } // if (icheck==0...                                                                                                                                                                 
    if (_xyz->isp==(_xyz->NpT-1)){
      _xyz->inewscanflg=1;
      _xyz->inewlayerflg=1;
      _xyz->isp=0;
      _xyz->indlayer+=1;
    }
  } // if (x<box[0...
}
  


void TempField::AnalyticTempCurr(double tcurr,std::vector<double> & TempOut, std::vector<int> &icellid, int Ntot)
{
  //computes temp field based on Schwalbach et al
  int j1,j2,j3,iplay;
  double x0,y0,x,y,z,dsq,dsq2,bx,by,xi,xp,yp,dirp,zp;
  std::vector<double> rij1(3),xs1(3),xs2(3);
  iplay=_xyz->nX[0]*_xyz->nX[1]*_xyz->ilaserLoc;
  TempOut.assign(Ntot,T0);
  xi = _xyz->tL*(1+std::numeric_limits<double>::epsilon() );
  int js1, js2;
  std::vector<double> a1m(6);
  for (int j=0;j<6;++j){a1m[j]=a1[j];}
  // x,y,z spatial location of source (in grid reference frame)
    x=_xyz->lcoor[2*ispvec[_xyz->isp]]-offset[0];
  y=_xyz->lcoor[2*ispvec[_xyz->isp]+1]-offset[1];
  z=_xyz->ilaserLoc*_xyz->dX[2]-offset[2];
  for (int j=0;j<Ntot;++j){
    j3 = floor(icellid[j]/(_xyz->nX[0]*_xyz->nX[1]));
    j2 = floor( (icellid[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
    j1 = icellid[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
    x0 = (double(j1)+.5)*(_xyz->dX[0]);
    y0 = (double(j2)+.5)*(_xyz->dX[1]);
    zp = (double(j3)+.5)*(_xyz->dX[2]);
    if (zp>z){continue;}
    xp=cos(_xyz->gth)*(x0-_xyz->LX[0]/2.)+
      sin(_xyz->gth)*(y0-_xyz->LX[1]/2.)+_xyz->LX[0]/2.;
    yp=-sin(_xyz->gth)*(x0-_xyz->LX[0]/2.) + 
      cos(_xyz->gth)*(y0-_xyz->LX[1]/2.)+_xyz->LX[1]/2.;
    if (fmod(_xyz->isp,_xyz->Nsd)==0){
      dirp=(_xyz->lcoor[2*ispvec[_xyz->isp+1]]-_xyz->lcoor[2*ispvec[_xyz->isp]]);
    } else {
      dirp=(_xyz->lcoor[2*ispvec[_xyz->isp]]-_xyz->lcoor[2*ispvec[_xyz->isp-1]]);
    }
    rij1[0] = xp-x;
    rij1[1] = yp-y;
    rij1[2] = zp-z;
    if (dirp*rij1[0]>0){ 
      //xp,yp is in front of laser
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
    } // if (dirp*rij1[0]>0...
  } // for (int j=0...
  // check if last simulation of scan
  double tmelt=_xyz->tL;
  int n1=_part->ncellLoc,icheck,ichecktmp;
  x=_xyz->lcoor2[2*ispvec[_xyz->isp]]-offset[0];
  y=_xyz->lcoor2[2*ispvec[_xyz->isp]+1]-offset[1];
  if (x<_xyz->gbox[0] || x>_xyz->gbox[1] || y<_xyz->gbox[2] || y>_xyz->gbox[3]){
    ichecktmp=std::any_of(TempOut.begin(),TempOut.begin()+n1,[&tmelt]
			  (double tchk){return tchk >= tmelt;});
    MPI_Allreduce(&ichecktmp,&icheck,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (icheck==0 || fmod(_xyz->isp+1,_xyz->Nsd)==0){ 
      _xyz->inewscanflg=1;
    } // if (icheck==0...
    if (_xyz->isp==(_xyz->NpT-1)){
      _xyz->inewscanflg=1;
      _xyz->inewlayerflg=1;
      _xyz->isp=0;
      _xyz->indlayer+=1;
    }
  } // if (x<box[0...
} // end AnalyticTempCurr()            
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
  iplay=_xyz->nX[0]*_xyz->nX[1]*_xyz->ilaserLoc;
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

