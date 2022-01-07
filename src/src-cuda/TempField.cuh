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
  void AnalyticalTempCalcMacro(Grid *dg,TempField *dtempF, double *dtempOut,
                                        double *dlcoor,double *dlcoor2,int *ispv,
                                        int &nThreads, int &nBlocks, int &Ntot);
  double bmSTD[3],bmLx[3],bmDX[3],offset[3],a1[6],bmV,bmP,T0,DelT,*TempCurr;
  int patternID,tInd,*ispvec;
}; // end class TempField

__global__ void calcTemptInd(const Grid *g,TempField *temp);
__global__ void analyticTempCalcPart1(Grid *dg,TempField *dtempF, int *ispv,
                                 double *dlcoor, double *dtempOut,
                                 int Ntot);
__global__ void analyticTempCalcPart2(Grid *dg, TempField *dtempT,int *ispv,
				      double *dlcoor2,double *dtempOut, int Ntot);
#endif
