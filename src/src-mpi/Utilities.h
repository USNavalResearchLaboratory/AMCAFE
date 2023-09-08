// temperature class (obtained from Moose simulation)

#ifndef UTILITIES_H
#define UTILITIES_H

#include "Grid.h"
#include <string>
#include <vector>
#include "Partition.h"
#include "BasePlate.h"
#include <gsl/gsl_integration.h>
#include <math.h>

class Utilities
{
 public:
  // define constructor
  Utilities(Grid &g, Partition &, BasePlate &);
  void Build();
  double rho_Ti64(double T);
    double cP_Ti64(double T);
    double k_Ti64(double T);
    double alpha_Ti64(double T);
    std::vector<double> Materials_Ti64(double T);
    double linter(double x, std::vector<double> ti, std::vector<double> yi, int imax);
    struct my_f_params;
    double rho_316(double T);
    double cP_316(double T);
    double k_316(double T);
    double alpha_316(double T);
    void Init_Util();
    void InitializeUT();
    double Integral(std::vector<double> L, double T, std::vector<double> R);
    double EASM_Temp_LG(std::vector<double> L, std::vector<double> R, double T);
    static double Temp_Fun_G(double x, void *p);
    static double Temp_fun(double x, void *p);
    std::vector<double> Materials_316L(double T);
    std::vector<double> bmSTD,bmLx,bmDX,offset,shiftL, dX;
    std::vector<int> nTTemp;
    int ilaserLoc, tInd, testvar;
    double bmV, DelT, Q, zlaserOff, bp_height, lat_shift;
    void LayerUpdate(double t);
    double tLayer;
    int tLC;    

 private:
    Grid *_xyz;
    Partition *_part;
    BasePlate *_bp;
    //std::string *filnambaseptr;
}; // end class TempField

#endif
