// temperature class (obtained from Moose simulation)

#ifndef PARAM_H
#define PARAM_H

#include <string>
#include <vector>

class Param
{
 public:
  // define constructor
  Param(std::string & filIn, int & myid, int & nprocs);
  void readFile(std::string & filIn, int & myid, int & nprocs);

  int patternID,outint;
  double bmV,bmP, layerT,tL,tS,bhatch,rnuc;
  std::vector<int> nX;
  std::vector<double> meltparam,dX;
}; // end class Param

#endif
