#ifndef VOXELSCA_CUH
#define VOXELSCA_CUH

#include "Grid.cuh"
#include "TempField.cuh"
#include "BasePlate.cuh"
#include "vector"
#include "string"
#include <random>

// #include <math.h>


class VoxelsCA
{
 public:
  VoxelsCA(Grid &);
  void WriteToHDF1(const std::string &filname, const Grid &, const double *);
  inline void loadS(std::vector<std::vector<double>>&S,std::vector<std::vector<int>> &sInd)
  {
    // this is for decentered octahedron method: 
    //      S is 6 corners of octahedron in local coor and sInd gives the 3 corner 
    //      indices for a given octant 
    S.assign(6,std::vector<double>(3));
    S[0] = {1,0,0};
    S[1] = {0,1,0};
    S[2] = {0,0,1};
    S[3] = {-1,0,0};
    S[4] = {0,-1,0};
    S[5] = {0,0,-1};
    sInd.assign(8,std::vector<int>(3));
    sInd[0] = {0,1,2};
    sInd[1] = {1,2,3};
    sInd[2] = {0,2,4};
    sInd[3] = {2,3,4};
    sInd[4] = {0,1,5};
    sInd[5] = {1,3,5};
    sInd[6] = {0,4,5};
    sInd[7] = {3,4,5};
  }// end inline void loadS
  inline void projectPointLine(double *A, double *x0, double *x1, double *xproj)
  {
    // computes the coordinate of projecting a point A on line connecting x0 and x1 = xproj
    // note that A,x0,x1,xproj are all length 3 arrays
    double n[3],snorm,d,t;
    snorm = pow(pow(x1[0]-x0[0],2)+pow(x1[1]-x0[1],2)+pow(x1[2]-x0[2],2),.5);
    n[0] = (x1[0]-x0[0])/snorm;
    n[1] = (x1[1]-x0[1])/snorm;
    n[2] = (x1[2]-x0[2])/snorm;
    d = -(n[0]*A[0]+n[1]*A[1]+n[2]*A[2]);
    t = (-d - n[0]*x0[0] - n[1]*x0[1] - n[2]*x0[2])/
      (n[0]*(x1[0]-x0[0])+n[1]*(x1[1]-x0[1])+n[2]*(x1[2]-x0[2]));
    xproj[0] = x0[0] + (x1[0]-x0[0])*t;
    xproj[1] = x0[1] + (x1[1]-x0[1])*t;
    xproj[2] = x0[2] + (x1[2]-x0[2])*t;
  } //end inline void projectPointLine...
  inline void loadRotMat(double omega, double *ax, std::vector<std::vector<double>> &rRot)
  {
    // loads the rotation matrix from (omega,ax), note that
    // ax is a 3x1  and rRot is a 3x3 static arrays
    rRot[0][0] = cos(omega) + pow(ax[0],2.0)*(1-cos(omega));
    rRot[0][1] = ax[0]*ax[1]*(1-cos(omega)) - ax[2]*sin(omega);
    rRot[0][2] = ax[0]*ax[2]*(1-cos(omega)) + ax[1]*sin(omega);
    rRot[1][0] = ax[0]*ax[1]*(1-cos(omega)) + ax[2]*sin(omega);
    rRot[1][1] = cos(omega) + pow(ax[1],2.0)*(1-cos(omega));
    rRot[1][2] = ax[1]*ax[2]*(1-cos(omega)) - ax[0]*sin(omega);
    rRot[2][0] = ax[2]*ax[0]*(1-cos(omega)) - ax[1]*sin(omega);
    rRot[2][1] = ax[2]*ax[1]*(1-cos(omega)) + ax[0]*sin(omega);
    rRot[2][2] = cos(omega) + pow(ax[2],2.0)*(1-cos(omega));
  } // end inline void loadRotMat
  inline void getNumPowderGrains(const Grid &g,int &numPG)
  {
    double rate = g.lrate* (g.dX[0]*1e6)*(g.dX[1]*1e6)*
      g.nX[0]*g.nX[1]* (g.layerT*1e6);
    std::poisson_distribution<int> ng(rate);
    numPG =0;
    while (numPG == 0){numPG = ng(genlayer);} // end while
  } // end getNumPowderGrains
  int *gID,*ineighID,*ineighptr,*vState;
  double *cTheta,*extents,*centroidOct;
  double vmax;
  int nGrain,seed0,seed1,NzhBP;
  std::default_random_engine genlayer;
}; // end class VoxelCA


__global__ void addlayer1part1(Grid *g,VoxelsCA *vx,double *xs, double *troids,
                          int *gD, int *vs, int *itmp, int numPG);

__global__ void addlayer1part2(Grid *g,VoxelsCA *vx, double *cTheta,
                          int *gD, int *itmp, int numPG);

__global__ void addlayer1part3(const Grid *g,int *gD, int *vs);

__global__ void copyGlobal(double *x1,double *x0, int n);

__global__ void getSites(Grid *g,VoxelsCA *vx,double *xs, int numPG);

void resizeGlobalArray(double **y, int &n0, int &n1);
void resizeArray(double **y, int &n);

#endif
