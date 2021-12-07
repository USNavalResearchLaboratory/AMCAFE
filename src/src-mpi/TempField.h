// temperature class (obtained from Moose simulation)

#ifndef TEMPFIELD_H
#define TEMPFIELD_H

#include "Grid.h"
#include <string>
#include <vector>
#include "Partition.h"
#include "BasePlate.h"
#include "Utilities.h"

class TempField
{
 public:
  // define constructor
  TempField(Grid &g, Partition &, BasePlate &, Utilities &);
  void InitializeSchwalbach();
  void InitializeAnalytic();
  void EASM(std::vector<double> &TempOut, std::vector<int> &icellid, int Ntot);
  void SchwalbachTempCurr();
  void AnalyticTempCurr(double tcurr,std::vector<double> &TempOut,std::vector<int> &icellid,int Ntot);    
  void SchwalbachTempCurr(double tcurr, std::vector<double> &TempOut);
  std::vector<std::vector<double>> Temp;
  std::vector<double> DDtTemp,dXM,TempCurr,lamXYZ,bmSTD,bmLx,bmX0,bmDX,bmPeriod,offset,shiftL;
  std::vector<int> nXM,nTTemp,ispvec;
  int NtM,indexM,patternID,nSource,tInd;
  double dtM,bmV,bmP,bmEta,rcut,tcut,T0,Ci,DelT,alpha,T0targ,zlaserOff,tBeg0;
  std::vector<double> a1,a2,Qp,tBeg; // 
 private:
  Grid *_xyz;
  Partition *_part;
  BasePlate *_bp;
  Utilities *_ut;
  std::string *filnambaseptr;
}; // end class TempField

#endif
