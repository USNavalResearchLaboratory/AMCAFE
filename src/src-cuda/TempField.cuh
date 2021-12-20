// temperature class (obtained from Moose simulation)

#ifndef TEMPFIELD_CUH
#define TEMPFIELD_CUH

#include "Grid.cuh"
#include <string>
#include <vector>

class TempField
{
 public:
  // define constructor
  TempField(const Grid &);
  void InitializeAnalytic(const Grid &);
  double bmSTD[3],bmLx[3],bmDX[3],offset[3],a1[6],bmV,bmP,T0,DelT,*TempCurr;
  int patternID,tInd,*ispvec;
}; // end class TempField

#endif
