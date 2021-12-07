// member functions for TempField

#include "Grid.h"
#include "Utilities.h"
#include "fstream"
#include "math.h"
#include "Partition.h"
#include "BasePlate.h"
#include "numeric"
#include <algorithm>
#include "mpi.h"
#include <gsl/gsl_integration.h>


// constructor
Utilities::Utilities(Grid &g, Partition &part, BasePlate &bp)
//Utilities::Utilities()
{
    _xyz = &g;
    _part = &part;
    _bp = &bp;
    double Q;
    double bmV;
    double bmP;
}

double Utilities::linter(double x, std::vector<double> ti, std::vector<double> yi, int imax){
    double y;
    int j;
    // if x is ouside the xi[] interval
    if (x <= ti[0])
        return y = yi[0];
    if (x >= ti[imax-1])
        return y = yi[imax-1];
    // loop to find j so that x[j-1] < x < x[j]
    j = 0;
    while (j <= imax-1){
        if (ti[j] >= x) break;
        j = j + 1;
    }
    y = yi[j-1]+(yi[j]-yi[j-1])*(x-ti[j-1])/(ti[j]-ti[j-1]);
    return y;
}

void Utilities::Build()
{
}

struct Utilities::my_f_params{
  double t;
  double bmV;
  std::vector<double> R;
  double T;
  double Q;
  double lat_shift;
  std::vector<double> M;
  std::vector<double> L;
};

double Utilities::rho_316(double T){
  if (T < 1660)
    return 8084.2 - 0.42086*T - 3.8942*pow(10,-5)*T*T;
  else if (T >= 1660 && T < 3090)
    return 7432.7 + 3.9338*pow(10,-2)*T - 1.8007*pow(10,-4)*T*T;
  else
    return 6124.72;
}

double Utilities::cP_316(double T){
  if (T < 1660)
    return 458.984 + 0.1328*T;
  else
    return 769.855;
}

double Utilities::k_316(double T){
  if (T < 1660)
    return 9.248 + 1.571*pow(10,-2)*T;
  else if (T >= 1660 && T < 3090)
    return 12.41 + 3.279*pow(10,-3)*T;
  else
    return 21.6465;
}

double Utilities::alpha_316(double T){
  return k_316(T)/(rho_316(T) * cP_316(T));
}

std::vector<double> Utilities::Materials_316L(double T){
    std::vector<double> M(4);
    M[0] = cP_316(T);
    M[2] = rho_316(T);
    M[1] = k_316(T);
    M[3] = alpha_316(T);
    return M;
}

void Utilities::LayerUpdate(double t){
  if (tLC != 0){
    if (tLayer == 0){tLayer = t;}
  }
  tLC = tLC + 1;
}

double Utilities::Integral(std::vector<double> L, double T, std::vector<double> R){
  double Temp;
  double t = 150.0e-5;
  double to = 0.0;
  //  std::vector<double> M = Materials_316L(T);
  std::vector<double> M = Materials_Ti64(T);
  struct my_f_params params = {t, bmV, R, T, Q, lat_shift, M, L};
  gsl_integration_fixed_workspace * w;
  const gsl_integration_fixed_type * Y = gsl_integration_fixed_legendre;

  gsl_function G;
  
  w = gsl_integration_fixed_alloc (Y, 15, to, t, 0.0, 0.0);
  G.function = &Utilities::Temp_Fun_G;
  G.params = &params;

  gsl_integration_fixed (&G, &Temp,  w);
  gsl_integration_fixed_free (w);
  return 475 + Temp;  // Make this T0 + Temp
}



double Utilities::EASM_Temp_LG(std::vector<double> L, std::vector<double> R, double T){
  double Told = T;
  int n=0, nlim=100;
  double tol = 10;
  double Tnew = Integral(L, T, R);
  while (abs(Told-Tnew)>tol && n < nlim){
    Told = Tnew;
    double T_temp = Integral(L, Told, R);
    Tnew = 0.5*(Told+T_temp);
    n++;
  }
  return Tnew;
}

void Utilities::InitializeUT()
{
  tLayer = 0;
  tLC = 0;
  bmV = _xyz->bmV;
  Q = _xyz->bmP;
  zlaserOff = 1.0; // 1.0 (this specifies where laser z value is - see SchwalbachTempCurr)
  bp_height = _bp->height; 
  DelT = _xyz->bmDelT;
  dX = _xyz->dX;
  bmDX = {DelT*bmV, _xyz->bhatch,_xyz->layerT}; // (SD, TD, BD) per layer                                                                                                                                                                  
  offset = _xyz->offset;

  // For Angled Structures:    ---------------------------------------------------------
  //  shiftL = {_xyz->meltparam[1], 0.0, 0.0}; // (SD, TD, BD) per layer
  //  lat_shift =_xyz->layerT*tan(_xyz->angle);
  bmLx = {_xyz->LX[0] + 2*offset[0], _xyz->LX[1], _xyz->LX[2]};
  //  if (_xyz->patternID == 9){bmLx[1]=_xyz->struct_width/sin(_xyz->angle);}  
  //  else if (_xyz->patternID == 24){bmLx[1] = _xyz->struct_width/sin(_xyz->angle);}
  //  ------------------------------------------------------------------------------------
  
  nTTemp = {_xyz->nTsd==2 ? _xyz->nTsd : int(ceil(bmLx[0]/bmDX[0] ))+1,
            bmDX[1]>_xyz->LX[1] ? 1 : int(ceil(bmLx[1]/bmDX[1]))+1,
            bmDX[2]<std::numeric_limits<double>::epsilon() ? 1: int(ceil(bmLx[2]/bmDX[2]))};
  bmLx={(nTTemp[0]-1)*bmDX[0],(nTTemp[1]-1)*bmDX[1],(nTTemp[2]-1)*bmDX[2]};
  //  printf("nTTemp[0] TESTING : \n   _xyz->nTsd = Hold   \n bmLx[0] = %f  \n bmDX[0] = %f \n  END nTTemp[0] TEST", bmLx[0], bmDX[0]);
} // end InitializeEASM   



double Utilities::rho_Ti64(double T){
    if (T < 24.0)
        return 4453.153;
    else if (T >= 24.0 && T < 300.0)
        return 4453 + 0.02869485*T - 6.448869*pow(10,-4)*pow(T,2) + 9.646377*pow(10,-7)*pow(T,3);
    else if (T >= 300.0 && T < 1033.0)
        return 4467.094 - 0.119171*T - 1.275079*pow(10,-5)*pow(T,2);
    else
        return 3684.36;

}


double Utilities::cP_Ti64(double T){
    if (T < 21.0)
        return 9.43019;
    else if (T >= 21.0 && T < 95.0)
    return 15.63498 - 2.331676*T + 0.1090673*pow(T,2) - 5.764251*pow(10,-4)*pow(T,3);
    else if (T >= 95.0 && T < 303.0)
        return -167.1732 + 6.754298*T - 0.02352377*pow(T,2) + 2.956257*pow(10,-5)*pow(T,3);
    else if (T >= 303.0 && T < 1850.0)
        return 383.3514 + 0.6708818*T - 5.35234*pow(10,-4)*pow(T,2) + 1.635172*pow(10,-7)*pow(T,3);
    else
        return 827.9;
}



double Utilities::k_Ti64(double T){
    if (T < 4.0)
        return 0.32;
    else if (T >= 4.0 && T < 311.0)
        return -0.08721429 + 0.1026071*T - 1.785714*pow(10,-4)*pow(T,2);
    else if (T >= 311.0 && T < 811.0)
        return 8.114005 - 0.01485211*T + 4.468662*pow(10,-5)*pow(T,2) - 2.273481*pow(10,-8)*pow(T,3);
    else
        return 13.333;
}



double Utilities::alpha_Ti64(double T){
   return k_Ti64(T)/(rho_Ti64(T) * cP_Ti64(T));
}


std::vector<double> Utilities::Materials_Ti64(double T){
    std::vector<double> M(4);
    M[0] = cP_Ti64(T);
    M[2] = rho_Ti64(T);
    M[1] = k_Ti64(T);
    M[3] = alpha_Ti64(T);
    return M;
}


double Utilities::Temp_Fun_G(double x, void *p){
  double sigma = 100.0e-6;
  double eta = 0.56;
  struct my_f_params *params = (struct my_f_params *)p;
  //From params struct
  double t = (params->t); //t is final time (set to 100.0e-4) whereas x is current time in integral
  std::vector<double> R = (params->R);
  double Q = (params->Q);
  double T = (params->T);
  double bmV = (params->bmV);
  double lat_shift = (params->lat_shift);
  std::vector<double> M = (params->M);
  std::vector<double> L = (params->L);
  double Lx, Ly, Lz;
  double C = eta*Q/(M_PI*M[2]*M[0]*sqrt(4*M_PI*M[3]));
  double x0 = L[0] - bmV*t; //Starting x for integral is wherever laser was t time ago
  Lx = x0+bmV*x; //x is current time(wrt integral)  here
  Ly = L[1]; // y position
  Lz = L[2]; // z position
  return C*pow(t-x,-0.5)/(2*M[3]*(t-x) + pow(sigma,2))*exp(-1*(pow(R[0]-Lx,2)+pow(R[1]-Ly,2)) /
                                                             (4*M[3]*(t-x)+2*pow(sigma,2)) - pow(abs(R[2]-Lz),2)/(4*M[3]*(t-x)));
}

