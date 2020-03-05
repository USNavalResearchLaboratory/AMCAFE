#ifndef VOXELSCA_H
#define VOXELSCA_H

#include "Grid.h"
#include "TempField.h"
#include "BasePlate.h"
#include "Partition.h"
#include "vector"
#include "string"
// #include <math.h>


class VoxelsCA
{
 public:
  VoxelsCA(Grid &, TempField &, Partition &);
  void ConvertSolid(const int &iswitch);
  void ConvertSolid1(const int &iswitch);
  void InitializeVoxels(BasePlate &bp);
  void InitializeTest1();
  void InitializeTest2();
  void SetLiquid();
  void SetLiquid2();
  void SetLiquid3();
  void SetLiquid4();
  void ZeroVoxels();
  void ZeroVoxels1();
  void ComputeVoxelCapture();
  void ComputeVoxelCapture2();
  void ComputeExtents();
  void UpdateVoxels();
  void UpdateVoxels2();
  void UpdateVoxels3();
  void UpdateVoxels4();
  void UpdateVoxels5();
  void UpdateVoxels6();
  void UpdateVoxels7();
  void UpdateLayer(std::string &filCSV);
  void NucleateGrains(std::vector<int> &nucInd, std::vector<double> &tnuc);
  void ComputeNucleation1();
  void ExtentsInitialize();
  void CheckTimeSkip();
  void WriteToVTU0(const std::string &filname);
  void WriteToVTU1(const std::string &filname);
  void WriteCSVData(const std::string &filname);
  void WriteCSVData1(const std::string &filname);
  void WriteToPVD(const std::string &filname, const std::vector<int> &, const std::vector<double> &);
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
  inline double getVelocity(double &tL,double &mL,double &kP,double &Gamma,double &c0,
			  double &T)
  {
    /*
    double v=(5.51*pow(M_PI,2.0)*pow((-mL)*(1-kP),1.5)*
       (Gamma))*( pow((tL - T),2.5)/pow(c0,1.5));
    */
    double v=(5.51*pow(M_PI,2.0)*pow((-mL)*(1-kP),1.5)*
       (Gamma))*( pow((tL - T),2.5)/pow(c0,1.5));
    return v;
  } // end inline void getVelocity

  std::vector<int> gID,gNucleus,ineighID,ineighptr;
  std::vector<int> vState; // 0=uninitialized; 1=liquid; 2=mushy; 3=solid
  std::vector<double> cTheta,extents,centroidOct;
  double vmax;
  int nGrain,seed0,seed1,NzhBP;
 private:
  double extentsInitialValue,extentsInitialValue2;
  Grid *_xyz;
  BasePlate *bp;
  TempField *_temp;
  Partition *_part;
  std::vector<int> nn;
}; // end class VoxelCA
#endif
