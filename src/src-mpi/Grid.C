// member function definitions for Grid.C

#include "Grid.h"
#include "math.h"
#include "iostream"
#include "string"
#include <algorithm>
#include <vector>
#include "mpi.h"
#include "fstream"
#include "sstream"
#include "numeric"

//constructor
Grid::Grid(std::string &filIn, int & myidIn, int & nprocsIn)
{
  // initialize default settings which will change if specified input file
  time =0.0;
  tInd =0;
  dt=0;
  filInput = filIn;
  myid = myidIn;
  nprocs = nprocsIn;
  nX = {128,128,64};
  dX = {1e-6,1e-6,1e-6};
  nDim = 3;
  outint = 1;
  patternID = 1;
  tL = 1620;
  tS = 1531.5;
  dTempM = (tL-tS)*.75; //7.5; // 2.5 // K (mean undercooling for nucleation)
  dTempS = (tL-tS)/3.0; //5.0; // 1.0 // K (standard dev undercooling for nucleation)
  mL = -10.9; // (K / wt%)
  dL = 3e-9; // (m^2/s)
  Gamma = 1e-7;  // (K m)
  muN = 9e-2; // 9 // 9e-2;
  dP = .48;
  c0 = 4.85; // (wt %)
  neighOrder = "first"; // can only equal "first"
  ntype = "Moore"; // either "Moore" or "VonNeumann"
  nnodePerCell = pow(2,nDim);
  rho = 8000.0; // kg /m^3
  cP = 502; // 502.0; // J/kg-K)
  kappa = 18; //18; //18.0; // W/(m-K)
  beamSTD = {7.5e-5,7.5e-5,7.5e-5};
  layerT = 25e-6;
  bpH = 0.0;
  mu = .01; // note that this is rate/ (\mu m)^3 for base plate tessellation
  beamEta =1.0;
  T0targ = 1500;
  T0 = 300;
  ictrl = 3;
  gsize={0,0};
  EASM_Flag = 0; // 0 for Analytical, 1 for EASM
  meltparam = {75e-6,162.75e-6,75e-6,75e-6};
  bhatch = 1.53*meltparam[2];
  bmDelT = 4.0/3.0*meltparam[0]/bmV;
  bmV = 500e-3;
  rNmax = .002;
  offset={0.0,0.0,0.0}; // positive value means starting outside of domain
  outNL = 0;
  gth0=0.;
  gth=0.;
  Avel=0.;
  nvel=0.;
  //read data from input file
  readInputFile();
  //if (gsize[0]==0){gsize={nX[0]*dX[0]*2,nX[1]*dX[1]*2};}
  LX = {nX[0]*dX[0],nX[1]*dX[1],nX[2]*dX[2]};
  bpH< std::numeric_limits<double>::epsilon() ? bpH=layerT : bpH=bpH;
  nZlayer = round(layerT/dX[2]);
  deltaXmin = *std::min_element(dX.begin(),dX.end());
  deltaTcheck = pow(deltaXmin,2.0)/dL;
  ethresh = .01;
  deltaThresh=.95;
  isp=0;
  inewscanflg=1;
  inewlayerflg=1;
  Nzhg = std::min(nX[2],int(ceil(bpH*.9999/dX[2])));
  indlayer=1;
  ilaserLoc= Nzhg + indlayer*nZlayer;
  Ntd=int(ceil(gsize[1]/bhatch))+1;
  Nsd=int(ceil(gsize[0]/(bmV*bmDelT)))+1;
  NpT=Nsd*Ntd;
  lcoor.assign(2*NpT,0.);
  lcoor2.assign(2*NpT,0.);
  double xlmin=(LX[0]-gsize[0])/2.,ylmin=(LX[1]-gsize[1])/2.;
  int k;
  double gmid[2];
  gmid[0]=LX[0]/2.;
  gmid[1]=LX[1]/2.;
  gth+=gth0;
  for (int j2=0;j2<Ntd;++j2){
    for (int j1=0;j1<Nsd;++j1){
      k=Nsd*j2+j1;
      lcoor[2*k]=j1*bmDelT*bmV+xlmin;
      lcoor[2*k+1]=j2*bhatch+ylmin;
      lcoor2[2*k] = cos(gth)*(lcoor[2*k]-gmid[0])-
	sin(gth)*(lcoor[2*k+1]-gmid[1])+gmid[0];
      lcoor2[2*k+1] = sin(gth)*(lcoor[2*k]-gmid[0])+
	cos(gth)*(lcoor[2*k+1]-gmid[1])+gmid[1];
    } // j1
  } // j2
  //grid box
  gbox[0]=-bhatch/2.;
  gbox[1]=LX[0]+bhatch/2.;
  gbox[2]=-bhatch/2.;
  gbox[3]=LX[1]+bhatch/2.;
} // end constructor
void Grid::UpdateLaser(){
  int itmp,iflg=0,irep=0;
  double x,y;

  // update laser to new scan if at last point in scan
  if (fmod(isp+1,Nsd)==0){inewscanflg=1;}
  // if at last step of layer, new layer gets updated here
  if (isp==(NpT-1)){
      inewscanflg=1;
      inewlayerflg=1;
      isp=0;
      indlayer+=1;
  }
  // otherwise, update laser location in general
  while(irep==0 || isp==0){
    irep+=1;
    itmp=isp;
    if (inewscanflg==0 && itmp<(NpT-1)){
      isp+=1;
    } else {
      while (itmp<(NpT-1) && iflg==0){
	itmp+=1;
	x=lcoor2[2*itmp];
	y=lcoor2[2*itmp+1];
	if (x>gbox[0] && x<gbox[1] && y>gbox[2] && y<gbox[3]){iflg=1;}
      } // while (itmp< ...
      if (itmp<(NpT-1)){
	isp=itmp;
	inewscanflg=0;
      } else {
	inewlayerflg=1;
	inewscanflg=1;
	indlayer+=1;
	ilaserLoc= Nzhg + indlayer*nZlayer;
	isp=0;
	// update grid
	double gmid[2];
	gmid[0]=LX[0]/2.;
	gmid[1]=LX[1]/2.;
	int k;
	gth+=gth0;
	for (int j2=0;j2<Ntd;++j2){
	  for (int j1=0;j1<Nsd;++j1){
	    k=Nsd*j2+j1;
	    lcoor2[2*k] = cos(gth)*(lcoor[2*k]-gmid[0])-
	      sin(gth)*(lcoor[2*k+1]-gmid[1])+gmid[0];
	    lcoor2[2*k+1] = sin(gth)*(lcoor[2*k]-gmid[0])+
	      cos(gth)*(lcoor[2*k+1]-gmid[1])+gmid[1];
	  } // j1
	} // j2
      } // if/else (itmp<NpT-1...
    } // if (inewscanflg==0...
  } // while(irep==0...
  if (irep==1){inewlayerflg=0;}
} // end UpdateLaser
void Grid::readInputFile()
{
  std::ifstream filIn;
  std::string inputData,keyword;
  int k=0,n=0;
  if (myid==0){
    filIn.open(filInput.c_str());
    while (!filIn.eof()){
      char c;
      filIn.get(c);
      if (c == '#') {
        // skip comments indicated by "#" in input file
	while (filIn && c != '\n') {
          filIn.get(c);
        }
      }
      if (k >= n) {
        n = 2*n;
        const int m = 10000;
        if (n < m) n = m;
        inputData.resize(n);
      }
      inputData[k++] = c;
    }
    filIn.close();
    inputData.resize(k);
  }
  MPI_Bcast(&k,1,MPI_INT,0,MPI_COMM_WORLD);
  inputData.resize(k);
  MPI_Bcast(&inputData[0],k,MPI_CHAR,0,MPI_COMM_WORLD);
  std::istringstream simInput(inputData);
  simInput >> keyword;
  while(simInput){
    std::transform(keyword.begin(),keyword.end(),keyword.begin(),
                   [](unsigned char s){return std::tolower(s);});
    if (keyword=="dx") {
      dX.resize(3);
      simInput >> keyword;
      dX[0]=std::stod(keyword);
      dX[1]=dX[0];
      dX[2]=dX[0];
    }
    if (keyword=="bmv") {
      simInput >> keyword;
      bmV=std::stod(keyword);
    }
    if (keyword=="bmp") {
      simInput >> keyword;
      bmP=std::stod(keyword);
    }
    if (keyword=="bmdelt") {
      simInput >> keyword;
      bmDelT=std::stod(keyword);
    }
    if (keyword=="nx") {
      nX.resize(3);
      simInput >> keyword;
      nX[0]=std::stoi(keyword);
      simInput >> keyword;
      nX[1]=std::stoi(keyword);
      simInput >> keyword;
      nX[2]=std::stoi(keyword);
    }
    if (keyword=="r") {
      simInput >> keyword;
      rNmax=std::stod(keyword);
    }
    if (keyword=="lt") {
      simInput >> keyword;
      layerT=std::stod(keyword);
    }
    if (keyword=="tl") {
      simInput >> keyword;
      tL=std::stod(keyword);
    }
    if (keyword=="ts") {
      simInput >> keyword;
      tS=std::stod(keyword);
    }
    if (keyword=="bhatch") {
      simInput >> keyword;
      bhatch=std::stod(keyword);
    }
    if (keyword=="mu") {
      simInput >> keyword;
      mu=std::stod(keyword);
    }
    if (keyword=="bph") {
      simInput >> keyword;
      bpH=std::stod(keyword);
    }
    if (keyword=="lrate") {
      simInput >> keyword;
      lrate=std::stod(keyword);
    }
    if (keyword=="outnl") {
      simInput >> keyword;
      outNL=std::stoi(keyword);
    }
    if (keyword=="ntsd") {
      simInput >> keyword;
      nTsd=std::stoi(keyword);
    }
    if (keyword=="offset") {
      simInput >> keyword;
      offset[0]=std::stod(keyword);
      simInput >> keyword;
      offset[1]=std::stod(keyword);
      simInput >> keyword;
      offset[2]=std::stod(keyword);
    }
    if (keyword=="EASM") {
      simInput >> keyword;
      EASM_Flag = std::stod(keyword);
    }
    if (keyword=="meltparam") {
      meltparam.resize(4);
      simInput >> keyword;
      meltparam[0]=std::stod(keyword);
      simInput >> keyword;
      meltparam[1]=std::stod(keyword);
      simInput >> keyword;
      meltparam[2]=std::stod(keyword);
      simInput >> keyword;
      meltparam[3]=std::stod(keyword);
    }
    if (keyword=="patternid") {
      simInput >> keyword;
      patternID=std::stoi(keyword);
    }
    if (keyword=="outint") {
      simInput >> keyword;
      outint=std::stoi(keyword);
    }
    if (keyword=="gridsize") {
      simInput >> keyword;
      gsize[0]=std::stod(keyword);
      simInput >> keyword;
      gsize[1]=std::stod(keyword);
    }
    if (keyword=="gridtheta") {
      simInput >> keyword;
      gth0=std::stod(keyword);
    }
    if (keyword=="dendritegrowth2") {
      simInput >> keyword;
      Avel=std::stod(keyword);
      simInput >> keyword;
      nvel=std::stod(keyword);
    }

    simInput >> keyword;
  } // while(simInput)
} // readInputFile
void Grid::UpdateTime(const double &vmax)
{
  dt = .25*std::min(deltaTcheck,deltaXmin/vmax);
  time +=dt;
  tInd +=1;
} // end UpdateTime

void Grid::UpdateTimeTest1(const double &vmax)
{
  //dt = .25*std::min(deltaTcheck,deltaXmin/vmax);
  dt = 1.50*std::min(deltaTcheck,deltaXmin/vmax);
  time +=dt;
  tInd +=1;
} // end UpdateTime

void Grid::UpdateTime2(const double &dtIn)
{
  dt = dtIn;
  time +=dt;
  tInd +=1;
} // end UpdateTime2

void Grid::SkipTime(const double &DelT)
{
  // DelT is the time increment of the temperature field
  //dt does not change
  time = ceil(time/DelT)*DelT ;
  //tInd +=1;
} // end UpdateTime2

void Grid::ComputeNeighborhood(int &j, std::string &nO,std::vector<int> &nn){
  // determines neighborhood of voxel j
  int nSize,nSize1;
  if (nO.compare("first")==0){
    ComputeNeighborhoodFirst(j,ntype,nn);
  } else if (nO.compare("second")==0){
    std::vector<int> nntmp;
    ComputeNeighborhoodFirst(j,ntype,nn);
    nSize = nn.size();
    for (int j1 =0;j1<nSize;++j1){
      ComputeNeighborhoodFirst(nn[j1],ntype,nntmp);
      nn.insert(nn.end(),nntmp.begin(),nntmp.end());
    } // for int j1
  } else if (nO.compare("third")==0){
    std::vector<int> nntmp1,nntmp2;
    ComputeNeighborhoodFirst(j,ntype,nn);
    nSize=nn.size();
    for (int j1 =0;j1<nSize;++j1){
      ComputeNeighborhoodFirst(nn[j1],ntype,nntmp1);
      nn.insert(nn.end(),nntmp1.begin(),nntmp1.end());
      for (int j2=0;j2<nntmp1.size();++j2){
	ComputeNeighborhoodFirst(nntmp1[j2],ntype,nntmp2);
	nn.insert(nn.end(),nntmp2.begin(),nntmp2.end());
      }
    } // for int j1
  }// end if nO.compare
  std::sort(nn.begin(),nn.end());
  std::vector<int>::iterator it = std::unique(nn.begin(),nn.end());
  nn.resize(std::distance(nn.begin(),it));
  nn.erase(std::remove(nn.begin(),nn.end(),j),nn.end());
}; // end ComputeNeighborhood

void Grid::ComputeNeighborhoodFirst(int &j,std::string &ntype, std::vector<int> &nn){
  // determines neighborhood of voxel j
  nn.assign(0,0);
  if (ntype.compare("Moore")==0){
    if (nDim ==2){
      int j2,j1,jst;
      j2 = floor(j/( nX[0] ));
      j1 = j - (nX[0])*j2;
      std::vector<int> itmp = {-1,0,1};
      for (int i2 =0;i2<3;++i2){
	if ( (j2+itmp[i2]<0) || (j2+itmp[i2]>=nX[1])){continue;}
	for (int i1 =0;i1<3;++i1){
	  if ( (j1+itmp[i1]<0) || (j1+itmp[i1]>=nX[0])){continue;}
	  jst = nX[0]*(j2+itmp[i2])+j1+itmp[i1];
	  if (jst !=j){nn.push_back(jst);}
	}
      }
    } else {
      int j3,j2,j1,jst;
      j3 = floor(j /( nX[0]*nX[1]) );
      j2 = floor( (j - nX[0]*nX[1]*j3)/nX[0] );
      j1 = j - nX[0]*nX[1]*j3 - nX[0]*j2;
      std::vector<int> itmp = {-1,0,1};
      for (int i3 =0;i3<3;++i3){
	if ( (j3+itmp[i3]<0) || (j3+itmp[i3]>=nX[2])){continue;}
	for (int i2 =0;i2<3;++i2){
	  if ( (j2+itmp[i2]<0) || (j2+itmp[i2]>=nX[1])){continue;}
	  for (int i1 =0;i1<3;++i1){
	    if ( (j1+itmp[i1]<0) || (j1+itmp[i1]>=nX[0])){continue;}
	    jst = nX[0]*nX[1]*(j3+itmp[i3])+nX[0]*(j2+itmp[i2])+j1+itmp[i1];
	    if (jst !=j){nn.push_back(jst);}
	  }
	}
      }
    } // if (nDim==2 ...
  } else { // this is for VonNeumann
    if (nDim ==2){
      int j2,j1,jst;
      j2 = floor(j/( nX[0] ));
      j1 = j - (nX[0])*j2;
      std::vector<int> itmp = {-1,1};
      for (int i1 =0;i1<2;++i1){
	if ( (j2+itmp[i1]<0) || (j2+itmp[i1]>=nX[1])){continue;}
	  jst = nX[0]*(j2+itmp[i1])+j1;
	  nn.push_back(jst);
      } // for int i1 ...
      for (int i1 =0;i1<2;++i1){
	if ( (j1+itmp[i1]<0) || (j1+itmp[i1]>=nX[0])){continue;}
	  jst = nX[0]*j2 +j1+itmp[i1];
	  nn.push_back(jst);
      } // for int i1 ...
    } else {
      int j3,j2,j1,jst;
      j3 = floor(j /( nX[0]*nX[1]) );
      j2 = floor( (j - nX[0]*nX[1]*j3)/nX[0] );
      j1 = j - nX[0]*nX[1]*j3 - nX[0]*j2;
      std::vector<int> itmp = {-1,1};
      for (int i1=0;i1<2;++i1){
	if ( (j3+itmp[i1]<0) || (j3+itmp[i1]>=nX[2])){continue;}
	jst = nX[0]*nX[1]*(j3+itmp[i1])+nX[0]*j2+j1;
	nn.push_back(jst);
      } // for (int i1...
      for (int i1=0;i1<2;++i1){
	if ( (j2+itmp[i1]<0) || (j2+itmp[i1]>=nX[1])){continue;}
	jst = nX[0]*nX[1]*j3+nX[0]*(j2+itmp[i1])+j1;
	nn.push_back(jst);
      } // for (int i1...
      for (int i1=0;i1<2;++i1){
	if ( (j1+itmp[i1]<0) || (j1+itmp[i1]>=nX[0])){continue;}
	jst = nX[0]*nX[1]*j3+nX[0]*j2+j1+itmp[i1];
	nn.push_back(jst);
      } // for (int i1...
    } // if (nDim==2 ...
  } // if (ntype.compare...

}; // end ComputeNeighborhoodFirst
