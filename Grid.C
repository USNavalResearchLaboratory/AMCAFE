// member function definitions for Grid.C

#include "Grid.h"
#include "math.h"
#include "iostream"
#include "string"
#include <algorithm>
#include <vector>

// constructor
Grid::Grid(const std::vector<double> & dxIn, const std::vector<int> & NxIn, 
	   const double &tLIn,
	   const double &mLIn,const double &c0In,const double &GammaIn,
	   const double &kPIn,const double &dLIn,const double &muNIn, 
	   const double & rhoIn, const double & cPIn, const double & kappaIn,
	   const double &layerThicknessIn, const std::string neighOrderIn,
	   const double &dTempMIn, const double &dTempSIn, 
	   const double &rNmaxIn, const int &nDimIn)
{
  // read in domain and  material parameters
  nDim = nDimIn;
  dX.resize(nDim);
  nX.resize(nDim);
  for (int j=0;j<nDim;++j){
    dX[j] = dxIn[j];
    nX[j] = NxIn[j];
  }
  tL = tLIn;
  mL = mLIn;
  c0 = c0In;
  Gamma = GammaIn;
  kP = kPIn;
  dL = dLIn;
  muN = muNIn;
  rho = rhoIn;
  cP = cPIn;
  kappa = kappaIn;
  dTempM = dTempMIn;
  dTempS = dTempSIn;
  rNmax = rNmaxIn;
  layerT = floor(layerThicknessIn/dX[2])*dX[2];
  deltaXmin = *std::min_element(dX.begin(),dX.end());
  deltaTcheck = pow(deltaXmin,2.0)/dL;
  time =0.0;
  tInd =0;
  dt=0;
  neighOrder = neighOrderIn;
  nnodePerCell = pow(2,nDim);

  // this should change into a user parameter
  ethresh = .01; 
  deltaThresh=.95;  
} // end constructor

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
    ComputeNeighborhoodFirst(j,nn);
  } else if (nO.compare("second")==0){
    std::vector<int> nntmp;
    ComputeNeighborhoodFirst(j,nn);
    nSize = nn.size();
    for (int j1 =0;j1<nSize;++j1){
      ComputeNeighborhoodFirst(nn[j1],nntmp);
      nn.insert(nn.end(),nntmp.begin(),nntmp.end());
    } // for int j1
  } else if (nO.compare("third")==0){
    std::vector<int> nntmp1,nntmp2;
    ComputeNeighborhoodFirst(j,nn);
    nSize=nn.size();
    for (int j1 =0;j1<nSize;++j1){
      ComputeNeighborhoodFirst(nn[j1],nntmp1);
      nn.insert(nn.end(),nntmp1.begin(),nntmp1.end());
      for (int j2=0;j2<nntmp1.size();++j2){
	ComputeNeighborhoodFirst(nntmp1[j2],nntmp2);
	nn.insert(nn.end(),nntmp2.begin(),nntmp2.end());
      }
    } // for int j1
  }// end if nO.compare
  std::sort(nn.begin(),nn.end());
  std::vector<int>::iterator it = std::unique(nn.begin(),nn.end());
  nn.resize(std::distance(nn.begin(),it));
  nn.erase(std::remove(nn.begin(),nn.end(),j),nn.end());
}; // end ComputeNeighborhood

void Grid::ComputeNeighborhoodFirst(int &j, std::vector<int> &nn){
  // This is the Moore neighborhood
  // determines neighborhood of voxel j
  nn.assign(0,0);
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
}; // end ComputeNeighborhoodFirst

