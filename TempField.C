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

void TempField::InitializeSchwalbach(int & patternIDIn, std::vector<double> & beamSTDIn, 
				     double & beamSpacingIn, double & beamVelocityIn,double & beamPowerIn,
				     double & beamEtaIn, std::vector<double> & LxAllIn, double & T0In)
{
  /*
    This follows the model: Schwalbach, Edwin J., et al. "A discrete source model of powder 
    bed fusion additive manufacturing thermal history." Additive Manufacturing 25 (2019): 485-498.
   */
  bmSTD = beamSTDIn;
  bmS = beamSpacingIn;
  bmV = beamVelocityIn;
  bmP = beamPowerIn;
  bmEta = beamEtaIn;
  patternID = patternIDIn;
  bmLx = LxAllIn;
  T0 = T0In;
  bmX0 = {(bmLx[0] - (_xyz->nX[0]*_xyz->dX[0]))/2.0,(bmLx[1] - (_xyz->nX[1]*_xyz->dX[1]))/2.0,
	  (bmLx[2] - (_xyz->nX[2]*_xyz->dX[2]))/2.0};
  double tb,minTemp;
  if (patternID==0){
    DelT = 4.0/3.0*bmSTD[0]/bmV;
    bmDX = {DelT*bmV,bmS,_xyz->layerT};
    bmPeriod = {bmLx[0]/bmV,(floor(bmLx[1]/bmS)+1)};
  }
  rcut = pow( -2* pow(*std::min_element(bmSTD.begin(),bmSTD.end()),2.0)*log(.001),.5);
  Ci = bmEta*bmP*DelT/(_xyz->rho*_xyz->cP*pow(2.0,.5)*pow(M_PI,1.5));
  alpha = _xyz->kappa/_xyz->cP/_xyz->rho;
  minTemp = 0.50; // cut of change in temperature for tcut
  tcut = pow(Ci/(2*alpha*minTemp),2.0/3.0);		       
} // end InitializeSchwalbach

void TempField::SchwalbachTempCurr()
{
  //computes temp field based on Schwalbach et al
  int Ntot = _part->nGhost+_part->ncellLoc, j1,j2,j3,iplay;
  double x0,y0,z0,rij,x,y,z,tc,tmin;
  std::vector<double> lam(3);
  //tmin = std::max(0.0,_xyz->time-tcut);
  tmin = _xyz->time-tcut;
  nSource = ceil( (_xyz->time - tmin) / DelT);
  ilaserLoc = _bp->Nzh + floor(floor(_xyz->time/(bmPeriod[1]*bmPeriod[0]))*_xyz->layerT/_xyz->dX[2]) + 1;
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
	x = bmV*fmod(tc,bmPeriod[0]);
	y = fmod(floor(tc/bmPeriod[0]),bmPeriod[1])*bmS;
	z = floor(tc/(bmPeriod[1]*bmPeriod[0]))*_xyz->layerT + (_bp->Nzh+1)*_xyz->dX[2];
	lam = {pow(bmSTD[0],2.0)+2*alpha*(_xyz->time - tc), 
	       pow(bmSTD[1],2.0)+2*alpha*(_xyz->time - tc),
	       pow(bmSTD[2],2.0)+2*alpha*(_xyz->time - tc)};
	rij = pow(x-x0,2.0)/2/lam[0] +pow(y-y0,2.0)/2/lam[1] +
	  pow(z-z0,2.0)/2/lam[2];
	TempCurr[j] += Ci/pow(lam[0]*lam[1]*lam[2],.5)*exp(-rij);
      } // for (int j1...
    } // for (int j...
  } // if (patternID==0 ...
} // end SchwalbachTempCurr()

void TempField::SchwalbachDDtTemp()
{
  //computes temp field time derivative based on Schwalbach et al
  int Ntot = _part->ncellLoc+_part->nGhost,j1,j2,j3,iplay;
  double x0,y0,z0,rij,x,y,z,tc,tmin,rdij;
  std::vector<double> lam(3);
  DDtTemp.assign(Ntot,0.0);
  tmin = std::max(0.0,_xyz->time-tcut);
  nSource = ceil( (_xyz->time - tmin) / DelT);
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
      x = bmV*fmod(tc,bmPeriod[0]);
      y = fmod(floor(tc/bmPeriod[0]),bmPeriod[1])*bmS;
      z = floor(tc/(bmPeriod[1]*bmPeriod[0]))*_xyz->layerT + (_bp->Nzh+1)*_xyz->dX[2];
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
  int Ntot = _part->nGhost+_part->ncellLoc;
  ilaserLoc = _xyz->nX.back(); // make z (3rd) direction the build direction   
  for (int j1=0;j1<(Ntot);++j1){
    TempCurr[j1] = 1584;
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
