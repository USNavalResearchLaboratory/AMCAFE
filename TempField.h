// temperature class (obtained from Moose simulation)

#ifndef TEMPFIELD_H
#define TEMPFIELD_H

#include "Grid.h"
#include <string>
#include <vector>
#include "Partition.h"
#include "BasePlate.h"

class TempField
{
 public:
  // define constructor
  TempField(Grid &g, Partition &, BasePlate &);
  void InitializeMoose(std::string & filnambase, const int &NtIn, 
		       const double &dtM, const std::vector<int> &nXM, const std::vector<double> &dXM);
  void InitializeSchwalbach(int & patternId, std::vector<double> & beamSTD, 
			    double & beamVelocity, double & beamPower, double & beamEta, std::vector<double> & Lx,
			    double & T0);
  void ReadCSVMoose1();
  void ReadCSVMoose2();
  void InterpolateToGrid1(const std::vector<double> &tempM,const std::vector<double> &xM,
			 const std::vector<double> &yM, const int &jt);
  void InterpolateToGrid2(const std::vector<double> &tempM, const std::vector<int> &,
			  std::vector<double> &);
  void ComputeDDtTempDiscrete();
  void ComputeDDtTempAnalytical();
  void SchwalbachTempCurr();
  void SchwalbachTempCurr(double tcurr, std::vector<double> &TempOut);
  void SchwalbachDDtTemp();
  void SchwalbachGradTemp();
  void Test2();
  std::vector<std::vector<double>> Temp;
  std::vector<double> DDtTemp,dXM,TempCurr,lamXYZ,bmSTD,bmLx,bmX0,bmDX,bmPeriod,offset;
  std::vector<int> nXM,nTTemp;
  int NtM,indexM,patternID,nSource,ilaserLoc,tInd;
  double dtM,bmV,bmP,bmEta,rcut,tcut,T0,Ci,DelT,alpha;

 private:
  Grid *_xyz;
  Partition *_part;
  BasePlate *_bp;
  std::string *filnambaseptr;
}; // end class TempField

#endif
