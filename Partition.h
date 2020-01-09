#ifndef PARTITION_H
#define PARTITION_H

#include "Grid.h"
#include <vector>

class Partition
{
 public:
  Partition(Grid &, int &, int &);
  void PartitionGraph();
  void PartitionGraph2();
  void PassInformation(std::vector<int> &);
  void PassInformation(std::vector<double> &);
  std::vector<int> icellidLoc,ineighProcId,ineigh2LocVals,ineigh2Locptr,
    iloc2NeighVals,iloc2Neighptr,ipointidLoc,iconnectivityLoc;
  std::vector<std::vector<int>> icellidLocAll; // this stores all the local icellidLoc in myid=0
  int ncellLoc, nGhost,nneighProc,myid,nprocs,npointLoc;
 private:
  Grid *_xyz;

};

#endif
