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
  void InitializeVoxels(BasePlate &bp);
  void InitializeTest1();
  void SetLiquid();
  void SetLiquid2();
  void ZeroVoxels();
  void ComputeVoxelCapture();
  void ComputeVoxelCapture2();
  void ComputeExtents();
  void UpdateVoxels();
  void UpdateVoxels2();
  void UpdateVoxels3();
  void NucleateGrain();
  void ComputeNucleation1();
  void ExtentsInitialize();
  void CheckTimeSkip();
  void WriteToVTU1(const std::string &filname);
  void WriteCSVData(const std::string &filname);
  void WriteCSVDataTest(const std::string &filname); // MUST DELETE AFTER TESTING
  void WriteToPVD(const std::string &filname, const std::vector<int> &, const std::vector<double> &);
  std::vector<int> gID,gNucleus;
  std::vector<int> vState; // 0=uninitialized; 1=liquid; 2=mushy; 3=solid
  std::vector<double> cTheta,extents;
  std::vector<double> velocityStore; // MUST DELETE AFTER TESTING
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
