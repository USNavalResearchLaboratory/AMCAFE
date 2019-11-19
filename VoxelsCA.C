// define member functions of VoxelCA

#include "Grid.h"
#include "BasePlate.h"
#include "VoxelsCA.h"
#include "Partition.h"
#include "iostream"
#include "fstream"
#include <math.h>
#include <algorithm>
#include <random>
#include <ctime>
#include "mpi.h"
#include "SampleOrientation.h"

// constructor
VoxelsCA::VoxelsCA(Grid &g,TempField &tf, Partition &part)
{ 
  _xyz = &g;
  _temp = &tf;
  _part = &part;
  int Ntot = (_part->ncellLoc)+(_part->nGhost);
  int Ntot1 = _part->ncellLoc;
  gID.resize(Ntot,0);
  vState.resize(Ntot,0);
  extentsInitialValue = 0.0; // (xyz->dx)/2.0
  extents.assign(Ntot1,extentsInitialValue);
  centroidOct.assign(3*Ntot1,0.0);
  seed0= 2132512;
  seed1=2912351;
} // end constructor

void VoxelsCA::InitializeVoxels(BasePlate &bp){
  // generate baseplate and add to voxel data
  bp.GenerateVoxel();
  int Ntot = (_part->ncellLoc)+(_part->nGhost);
  nGrain = bp.Ngrain;
  NzhBP = bp.Nzh;
  gNucleus.resize(bp.Ngrain,0);
  for (int j = 0;j < bp.Ngrain;++j){ gNucleus[j] = bp.gNucleus[j];}

  for (int j=0; j<bp.jVals.size() ;++j){
    gID[j] = bp.X[j];
    if (  _part->icellidLoc[j] < (bp.Nzh - 1)*_xyz->nX[0]*_xyz->nX[1]){
      vState[j] = 3; // base plate to be solid except top layer
    } else {
      vState[j] = 2; // baseplate top layer to be mushy
      //vState[j] = 0; // baseplate top layer to be mushy (will be enforced after it gets melted)
    }
  } // for j
  cTheta.assign(bp.Ngrain*4,0);
  for (int j1=0;j1<bp.Ngrain*4;++j1){cTheta[j1] = bp.cTheta[j1];} // end for j1
  // share vState information among processors and fill ghost cell information
  _part->PassInformation(vState);
  _part->PassInformation(gID);
  int i2 = _xyz->nX[0], i3 = _xyz->nX[0]*_xyz->nX[1],jn[3];
  for (int j=0;j<_part->ncellLoc;++j){
    if (gID[j]==0){continue;}
    jn[2] = floor(_part->icellidLoc[j]/i3);
    jn[1] =floor( (_part->icellidLoc[j]- i3*jn[2])/i2);
    jn[0] = _part->icellidLoc[j] -i3*jn[2] - i2*jn[1];
    centroidOct[3*j]=(double(jn[0])+.5)*_xyz->dX[0];
    centroidOct[3*j+1]=(double(jn[1])+.5)*_xyz->dX[1];
    centroidOct[3*j+2]=(double(jn[2])+.5)*_xyz->dX[2];
  } // for (int j ...
} // end InitializeVoxles

void VoxelsCA::InitializeTest1()
{
  int Ntot = _part->ncellLoc + _part->nGhost,jst,j2;
  nGrain = 1;
  gID.assign(Ntot,0);
  jst = _xyz->nX[0]*_xyz->nX[1]*floor(_xyz->nX[2]/2.0)+_xyz->nX[0]*floor(_xyz->nX[1]/2.0)+floor(_xyz->nX[0]/2.0);
  gNucleus.resize(1,jst);
  cTheta = {M_PI*60.0/180.0,1.0/pow(3.0,.5),1.0/pow(3.0,.5),1.0/pow(3.0,.5)};
  for (int j =0; j <Ntot; ++j){ vState[j] = 1;}
  j2 = std::distance(_part->icellidLoc.begin(),
		     std::find(_part->icellidLoc.begin(),_part->icellidLoc.end(),jst));
  if (j2 != _part->icellidLoc.size()){
    vState[j2] = 2;
    gID[j2] = 1;
  }
}

void VoxelsCA::UpdateVoxels()
{
  // set liquid if temperature > TL
  SetLiquid();
  // zero voxels above laser (no material there)
  ZeroVoxels();
  // SetLiquids, ZeroVoxels: no need to pass information for  b/c only local evaluation)
  // solid (vState=3) to mushy (vState=2) if one neighbor liquid (vState=1)
  ConvertSolid(1);
  _part->PassInformation(vState);
  // update extents of growing grains 
  ComputeExtents();
  _part->PassInformation(extents);
  // (ComputeExtents: no information pass: only local evaluation)
  // determine if liquid voxel is captured by neighbor
  ComputeVoxelCapture();
  _part->PassInformation(vState);
  _part->PassInformation(gID);
  _part->PassInformation(extents);
  // mushy (vState=2) to solid (vState=3) if all neighbors 2 
  ConvertSolid(0);
  _part->PassInformation(vState);
  _part->PassInformation(extents);
};
void VoxelsCA::UpdateVoxels2()
{
  /*
    This is the modified updating scheme that is described in notes
   */
  // set liquid if temperature > TL
  SetLiquid2();
  // zero voxels above laser (no material there)
  ZeroVoxels();
  // SetLiquids, ZeroVoxels: no need to pass information for  b/c only local evaluation)
  // solid (vState=3) to mushy (vState=2) if one neighbor liquid (vState=1)
  ConvertSolid(1);
  _part->PassInformation(vState);
  // capture all undercooled liquid voxels by growing grains
  int i0,iT=10;
  while (iT>0){
    ComputeVoxelCapture2();
    _part->PassInformation(vState);
    _part->PassInformation(gID);
    // mushy (vState=2) to solid (vState=3) if all neighbors 2 
    ConvertSolid(0);
    _part->PassInformation(vState);
    i0=0;
    for (int j=0;j<_part->ncellLoc;++j){
      if (vState[j]==1 && _temp->TempCurr[j]< _xyz->tL){i0+=1;}
    }
    MPI_Allreduce(&i0,&iT,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  } // while (iT>0...

};
void VoxelsCA::UpdateVoxels3()
{
  /*
    This is the modified updating scheme that is described in notes
   */
  // set liquid if temperature > TL 
  SetLiquid2();
  // zero voxels above laser (no material there)
  ZeroVoxels();
  // solid (vState=3) to mushy (vState=2) if one neighbor liquid (vState=1)
  ConvertSolid(1);
  _part->PassInformation(vState);
  // assemble arrays to be used to compute grain growth
  int NlocA=0,Na=0,NvecA[_part->nprocs],Ntot,cc,cc1,jn[3],itmp[3]={-1,0,1};
  Ntot = _part->ncellLoc;
  // start:determine nucleation: location and time
  std::vector<int> nucInd; // index of where nucleation
  std::vector<double> tnuc; // time within [t,t+DT] when nucleation occurs
  NucleateGrains(nucInd,tnuc);
  int Nnuc = nucInd.size(),Nnucvec[_part->nprocs],NnucA=0, jnuc0[_part->nprocs];
  std::vector<int> :: iterator i1p;
  MPI_Allgather(&Nnuc,1,MPI_INT,&Nnucvec[0],1,MPI_INT,MPI_COMM_WORLD);
  for (int j=0;j<_part->nprocs;++j){NnucA+=Nnucvec[j];}
  jnuc0[0]=0;
  for (int j=1;j<_part->nprocs;++j){jnuc0[j] = jnuc0[j-1]+Nnucvec[j-1];}
  int nucA[NnucA];
  double tnucA[NnucA];
  // end:determine nucleation: location and time
  std::vector<std::vector<double>> aa;
  SampleOrientation sa;
  unsigned int sdloc; 
  double nA[3] = {1.0/pow(3.0,.5),1.0/pow(3.0,.5),1.0/pow(3.0,.5)},
    pA[3] = {1.0,0.0,0.0};
  for (int j=0;j<Ntot;++j){
    if (vState[j]==2 || (vState[j]==1 && _temp->TempCurr[j]< _xyz->tL)){NlocA+=1;}
  }
  MPI_Allgather(&NlocA,1,MPI_INT,&NvecA[0],1,MPI_INT,MPI_COMM_WORLD);
  for (int j=0;j<_part->nprocs;++j){Na+=NvecA[j];}
  double T[Na];
  int vS[Na],V[Na][24],G[Na],vI[Na],j0[_part->nprocs],
    Nneigh[Na],vneigh[26],iv[_part->nprocs+1];
  std::fill(Nneigh,Nneigh+Na,0); // initialize all to value 0
  j0[0]=0;
  for (int j=1;j<_part->nprocs;++j){j0[j]= j0[j-1]+NvecA[j-1];}
  cc=0;
  cc1=0;
  for (int j=0;j<Ntot;++j){
    if (vState[j]==2 || (vState[j]==1 && _temp->TempCurr[j]< _xyz->tL)){
      vS[j0[_part->myid]+cc] = vState[j];
      vI[j0[_part->myid]+cc] = _part->icellidLoc[j];
      G[j0[_part->myid]+cc] = gID[j];
      T[j0[_part->myid]+cc] = _temp->TempCurr[j];
      i1p=std::find(nucInd.begin(),nucInd.end(),j);
      if (i1p !=nucInd.end()){
	nucA[jnuc0[_part->myid]+cc1] = cc + j0[_part->myid];
	tnucA[jnuc0[_part->myid]+cc1] = tnuc[std::distance(nucInd.begin(),i1p)];
	cc1+=1;
      } // if (std::any_of ...
      cc+=1;
    } // if (vState[j]...
  } // for (int j...
  for (int j=0;j<_part->nprocs;++j){
    if (NvecA[j]==0){continue;}
    MPI_Bcast(&vS[j0[j]],NvecA[j],MPI_INT,j,MPI_COMM_WORLD);
    MPI_Bcast(&vI[j0[j]],NvecA[j],MPI_INT,j,MPI_COMM_WORLD);
    MPI_Bcast(&G[j0[j]],NvecA[j],MPI_INT,j,MPI_COMM_WORLD);
    MPI_Bcast(&T[j0[j]],NvecA[j],MPI_DOUBLE,j,MPI_COMM_WORLD);
  } // for (int j ..
  for (int j=0;j<_part->nprocs;++j){
    if (Nnucvec[j]==0){continue;}
    MPI_Bcast(&nucA[jnuc0[j]],Nnucvec[j],MPI_INT,j,MPI_COMM_WORLD);
    MPI_Bcast(&tnucA[jnuc0[j]],Nnucvec[j],MPI_DOUBLE,j,MPI_COMM_WORLD);
  } // for (int j ..
  int i1=_xyz->nX[0],i2=_xyz->nX[0]*_xyz->nX[1],jst;
  for (int j=0;j<Na;++j){
    std::fill(vneigh,vneigh+26,-1);
    jn[2] = floor(vI[j]/i2);
    jn[1] = floor((vI[j]- i2*jn[2])/i1);
    jn[0] = vI[j] - i2*jn[2] - i1*jn[1];
    cc=0;
    for (int j3=0;j3<3;++j3){
      if ( (jn[2]+itmp[j3]<0) || (jn[2]+itmp[j3]>=_xyz->nX[2])){continue;}      
      for (int j2=0;j2<3;++j2){
	if ( (jn[1]+itmp[j2]<0) || (jn[1]+itmp[j2]>=_xyz->nX[1])){continue;}      
	for (int j1=0;j1<3;++j1){
	  if ( (jn[0]+itmp[j1]<0) || (jn[0]+itmp[j1]>=_xyz->nX[0])){continue;}      
          jst = i2*(jn[2]+itmp[j3])+i1*(jn[1]+itmp[j2])+jn[0]+itmp[j1];
	  if (jst !=vI[j]){vneigh[cc]=jst; cc+=1;}
	} // for (int j1...
      } // for (int j2...
    } // for (int j3...    
    for (int j1=0;j1<26;++j1){
      cc=std::distance(vI,std::find(vI,vI+Na,vneigh[j1]));
      if (cc <Na){
	V[j][Nneigh[j]] = cc;
	Nneigh[j]+=1;
      } // if (cc < ...
    } // for (int j1...
  } // for (int j=0...
  i1 = ceil( (double)Na / (double)_part->nprocs);
  i2 = floor( (double)Na / (double)i1);
  for (int j=0;j<(_part->nprocs+1);++j){
    if (j < (i2+1)){iv[j] = i1*j;}
    if (j>i2 && j<_part->nprocs){iv[j] = i2*i1 + 
	floor( (double)(Na-i2*i1)/(double)(_part->nprocs-i2));}
    if (j==_part->nprocs){iv[j]=Na;}
  } // for (int j...
  // end assemble arrays to be used to compute grain growth
  // capture all undercooled liquid voxels by growing grains
  int js,j1s,jx[3],jy[3],jTs[_part->nprocs], j1Ts[_part->nprocs],i3,countS=Na;
  double DtT[_part->nprocs],velX,velY,omega,ax[3],
    rRot[3][3],dlocX[3],dlocYnorm,dlocY[3],dr,vhat,tmn1[Na],
    timeUntil,dnx[3],dny[3],tmp1,tinc=0.0;
  std::fill(tmn1,tmn1+Na,0.0);
  i2 = _xyz->nX[0]; i3 = _xyz->nX[0]*_xyz->nX[1];
  cc=0;
  while (std::any_of(vS,vS+Na,[](int n){return n==1;}) && std::count(vS,vS+Na,2)!=countS) {
    double DtMin=1e6;
    countS = std::count(vS,vS+Na,2);
    for (int j=iv[_part->myid];j<iv[_part->myid+1];++j){
      if (vS[j]==1){
	jx[2] = floor(vI[j]/i3);
	jx[1] = floor((vI[j]- i3*jx[2])/i2);
	jx[0] = vI[j] - i3*jx[2] - i2*jx[1];
	velX=(5.51*pow(M_PI,2.0)*pow((- _xyz->mL)*(1-_xyz->kP),1.5)*
		 (_xyz->Gamma))*( pow((_xyz->tL - T[j]),2.5)/pow(_xyz->c0,1.5));
	for (int j1=0;j1<Nneigh[j];++j1){
	  if (vS[V[j][j1]]==1){continue;}
	  jy[2] = floor(vI[V[j][j1]]/i3);
	  jy[1] =floor( (vI[V[j][j1]]- i3*jy[2])/i2);
	  jy[0] = vI[V[j][j1]] -i3*jy[2] - i2*jy[1];
	  jn[2] = floor(gNucleus[G[V[j][j1]]-1]/i3);
	  jn[1] =floor( (gNucleus[G[V[j][j1]]-1]- i3*jn[2])/i2);
	  jn[0] = gNucleus[G[V[j][j1]]-1] -i3*jn[2] - i2*jn[1];
	  dnx[0] = (jx[0]-jn[0])*_xyz->dX[0];
	  dnx[1] = (jx[1]-jn[1])*_xyz->dX[1];
	  dnx[2] = (jx[2]-jn[2])*_xyz->dX[2];
	  dny[0] = (jy[0]-jn[0])*_xyz->dX[0];
	  dny[1] = (jy[1]-jn[1])*_xyz->dX[1];
	  dny[2] = (jy[2]-jn[2])*_xyz->dX[2];
	  tmp1 = (dnx[0]*dny[0]+dnx[1]*dny[1]+dnx[2]*dny[2])/
	    (dnx[0]*dnx[0]+dnx[1]*dnx[1]+dnx[2]*dnx[2]);
	  omega = cTheta[4*(G[V[j][j1]]-1)];
	  ax[0]=cTheta[4*(G[V[j][j1]]-1)+1];
	  ax[1]=cTheta[4*(G[V[j][j1]]-1)+2];
	  ax[2]=cTheta[4*(G[V[j][j1]]-1)+3];
	  // matrix is local->global; need to multiply by transpose for global->local            
	  rRot[0][0] = cos(omega) + pow(ax[0],2.0)*(1-cos(omega));
	  rRot[0][1] = ax[0]*ax[1]*(1-cos(omega)) - ax[2]*sin(omega);
	  rRot[0][2] = ax[0]*ax[2]*(1-cos(omega)) + ax[1]*sin(omega);
	  rRot[1][0] = ax[0]*ax[1]*(1-cos(omega)) + ax[2]*sin(omega);
	  rRot[1][1] = cos(omega) + pow(ax[1],2.0)*(1-cos(omega));
	  rRot[1][2] = ax[1]*ax[2]*(1-cos(omega)) - ax[0]*sin(omega);
	  rRot[2][0] = ax[2]*ax[0]*(1-cos(omega)) - ax[1]*sin(omega);
	  rRot[2][1] = ax[2]*ax[1]*(1-cos(omega)) + ax[0]*sin(omega);
	  rRot[2][2] = cos(omega) + pow(ax[2],2.0)*(1-cos(omega));
	  // put into 1st quadrant b/c of symmetry     
	  dlocX[0] = std::fabs(rRot[0][0]*dnx[0]+rRot[1][0]*dnx[1]+rRot[2][0]*dnx[2]);
	  dlocX[1] = std::fabs(rRot[0][1]*dnx[0]+rRot[1][1]*dnx[1]+rRot[2][1]*dnx[2]);
	  dlocX[2] = std::fabs(rRot[0][2]*dnx[0]+rRot[1][2]*dnx[1]+rRot[2][2]*dnx[2]);
	  dlocY[0] = std::fabs(rRot[0][0]*dny[0]+rRot[1][0]*dny[1]+rRot[2][0]*dny[2]);
	  dlocY[1] = std::fabs(rRot[0][1]*dny[0]+rRot[1][1]*dny[1]+rRot[2][1]*dny[2]);
	  dlocY[2] = std::fabs(rRot[0][2]*dny[0]+rRot[1][2]*dny[1]+rRot[2][2]*dny[2]);
	  dlocYnorm = pow(pow(dlocY[0],2.0)+pow(dlocY[1],2.0)+pow(dlocY[2],2.0),.5);
          dr = - (nA[0]*pA[0]+nA[1]*pA[1]+nA[2]*pA[2])*dlocYnorm;
	  dr = std::fabs(nA[0]*dlocX[0]+nA[1]*dlocX[1]+nA[2]*dlocX[2] + dr);	  
	  velY=(5.51*pow(M_PI,2.0)*pow((- _xyz->mL)*(1-_xyz->kP),1.5)*
		 (_xyz->Gamma))*( pow((_xyz->tL - T[V[j][j1]]),2.5)/pow(_xyz->c0,1.5));
	  vhat = .5*(velY+velX)/pow(3.0,.5);
	  timeUntil = (dr-vhat*tmn1[V[j][j1]])/vhat;
	  if (timeUntil < DtMin){
	    DtMin = timeUntil;
	    js = j;
	    j1s = j1;
	  } // if (timeUntil < ...
	} // for (int j1
      } // if (vS[j...
    } // for (int j...
    MPI_Allgather(&DtMin,1,MPI_DOUBLE,&DtT[0],1,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Allgather(&js,1,MPI_INT,&jTs[0],1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&j1s,1,MPI_INT,&j1Ts[0],1,MPI_INT,MPI_COMM_WORLD);
    DtMin = *std::min_element(DtT,DtT+_part->nprocs);
    double DtMin2 = 1.01*DtMin;
    for (int j=0;j<NnucA;++j){
      if ( (tnucA[j] > tinc) && (tnucA[j] <= tinc+DtMin)){
	if (DtMin2> tnucA[j]-tinc){
	  DtMin2 = tnucA[j]-tinc;
	  js=nucA[j];
	} // if (DtMin2...
      } // if ( (tnucA ...
    } // for int j ...
    if (DtMin2 < DtMin){
      vS[js] = 2;
      nGrain+=1;
      G[js] = nGrain;
      sdloc= seed0 + 32*_xyz->tInd +64*nGrain;
      sa.GenerateSamples(1,sdloc,aa);
      cTheta.insert(cTheta.end(), &aa[0][0],&aa[0][0]+4);
      gNucleus.push_back(vI[js]);
      for (int j=0;j<Na;++j){tmn1[j]+=DtMin2;}
      tmn1[js] = 0.0;
    } else {
      i1 = std::distance(DtT,std::min_element(DtT,DtT+_part->nprocs));
      js = jTs[i1];
      j1s = j1Ts[i1];
      vS[js] = 2;
      G[js] = G[V[js][j1s]];
      for (int j=0;j<Na;++j){tmn1[j]+=DtMin;}
      tmn1[js] = 0.0;
    } // if (DtMin2 < ...
    cc+=1;
    tinc += std::min(DtMin2,DtMin);
  } // while (std::any
  // end capture all undercooled liquid voxels by growing grains
  // bring global variables back to local variables
  cc1=0;
  for (int j=0;j<Ntot;++j){
    if (vState[j]==2 || (vState[j]==1 && _temp->TempCurr[j]< _xyz->tL)){
      vState[j] = vS[j0[_part->myid]+cc1];
      gID[j] = G[j0[_part->myid]+cc1];
      cc1+=1;
    } // if (vState[j]...
  } // for (int j...
  // pass information 
  for (int j=Ntot;j<_part->nGhost+Ntot;++j){
    cc1=std::distance(vI,std::find(vI,vI+Na,_part->icellidLoc[j]));
    if (cc1<Na){
      vState[j] = vS[cc1];
      gID[j] = G[cc1];
    } // if (cc1<
  } // for (int j...
  // mushy (vState=2) to solid (vState=3) if all neighbors 2 
  ConvertSolid(0);
  _part->PassInformation(vState);
}; // end UpdateVoxels3
void VoxelsCA::UpdateVoxels4()
{
  /*
    This is the version that includes the decentered octahedral method
   */
  // set liquid if temperature > TL 
  SetLiquid3();
  // zero voxels above laser (no material there)
  ZeroVoxels1();
  // solid (vState=3) to mushy (vState=2) if one neighbor liquid (vState=1)
  ConvertSolid1(1);
  _part->PassInformation(vState);
  // assemble arrays to be used to compute grain growth
  int NlocA=0,Na=0,NvecA[_part->nprocs],Ntot,cc,cc1,jn[3];
  Ntot = _part->ncellLoc;
  // start:determine nucleation: location and time
  std::vector<int> nucInd; // index of where nucleation
  std::vector<double> tnuc; // time within [t,t+DT] when nucleation occurs
  NucleateGrains(nucInd,tnuc);
  int Nnuc = nucInd.size(),Nnucvec[_part->nprocs],NnucA=0, jnuc0[_part->nprocs];
  std::vector<int> :: iterator i1p;
  MPI_Allgather(&Nnuc,1,MPI_INT,&Nnucvec[0],1,MPI_INT,MPI_COMM_WORLD);
  for (int j=0;j<_part->nprocs;++j){NnucA+=Nnucvec[j];}
  jnuc0[0]=0;
  for (int j=1;j<_part->nprocs;++j){jnuc0[j] = jnuc0[j-1]+Nnucvec[j-1];}
  int nucA[NnucA];
  double tnucA[NnucA];
  // end:determine nucleation: location and time
  std::vector<std::vector<double>> aa,sdiag0,sdiag(6,std::vector<double>(3));
  std::vector<std::vector<int>> sInd;
  loadS(sdiag0,sInd); // this is for the decentered octohedron algorithm
  SampleOrientation sa;
  unsigned int sdloc; 
  double nA[3] = {1.0/pow(3.0,.5),1.0/pow(3.0,.5),1.0/pow(3.0,.5)},
    pA[3] = {1.0,0.0,0.0};
  for (int j=0;j<Ntot;++j){
    if (vState[j]==2 || (vState[j]==1 && _temp->TempCurr[j]< _xyz->tL)){NlocA+=1;}
  }
  MPI_Allgather(&NlocA,1,MPI_INT,&NvecA[0],1,MPI_INT,MPI_COMM_WORLD);
  for (int j=0;j<_part->nprocs;++j){Na+=NvecA[j];}
  double T[Na],ExtA[Na],CentroidA[3*Na];
  int vS[Na],V[Na][24],G[Na],vI[Na],j0[_part->nprocs],
    Nneigh[Na],iv[_part->nprocs+1];
  std::fill(Nneigh,Nneigh+Na,0); // initialize all to  0
  std::fill(ExtA,ExtA+Na,0.0); // initialize all to 0 for each time step
  j0[0]=0;
  for (int j=1;j<_part->nprocs;++j){j0[j]= j0[j-1]+NvecA[j-1];}
  cc=0;
  cc1=0;
  for (int j=0;j<Ntot;++j){
    if (vState[j]==2 || (vState[j]==1 && _temp->TempCurr[j]< _xyz->tL)){
      vS[j0[_part->myid]+cc] = vState[j];
      vI[j0[_part->myid]+cc] = _part->icellidLoc[j];
      G[j0[_part->myid]+cc] = gID[j];
      T[j0[_part->myid]+cc] = _temp->TempCurr[j];
      //ExtA[j0[_part->myid]+cc] = extents[j];
      CentroidA[3*j0[_part->myid]+3*cc]=centroidOct[3*j];
      CentroidA[3*j0[_part->myid]+3*cc+1]=centroidOct[3*j+1];
      CentroidA[3*j0[_part->myid]+3*cc+2]=centroidOct[3*j+2];
      i1p=std::find(nucInd.begin(),nucInd.end(),j);
      if (i1p !=nucInd.end()){
	nucA[jnuc0[_part->myid]+cc1] = cc + j0[_part->myid];
	tnucA[jnuc0[_part->myid]+cc1] = tnuc[std::distance(nucInd.begin(),i1p)];
	cc1+=1;
      } // if (std::any_of ...
      cc+=1;
    } // if (vState[j]...
  } // for (int j...
  for (int j=0;j<_part->nprocs;++j){
    if (NvecA[j]==0){continue;}
    MPI_Bcast(&vS[j0[j]],NvecA[j],MPI_INT,j,MPI_COMM_WORLD);
    MPI_Bcast(&vI[j0[j]],NvecA[j],MPI_INT,j,MPI_COMM_WORLD);
    MPI_Bcast(&G[j0[j]],NvecA[j],MPI_INT,j,MPI_COMM_WORLD);
    MPI_Bcast(&T[j0[j]],NvecA[j],MPI_DOUBLE,j,MPI_COMM_WORLD);
    //MPI_Bcast(&ExtA[j0[j]],NvecA[j],MPI_DOUBLE,j,MPI_COMM_WORLD);
    MPI_Bcast(&CentroidA[3*j0[j]],3*NvecA[j],MPI_DOUBLE,j,MPI_COMM_WORLD);
  } // for (int j ..
  for (int j=0;j<_part->nprocs;++j){
    if (Nnucvec[j]==0){continue;}
    MPI_Bcast(&nucA[jnuc0[j]],Nnucvec[j],MPI_INT,j,MPI_COMM_WORLD);
    MPI_Bcast(&tnucA[jnuc0[j]],Nnucvec[j],MPI_DOUBLE,j,MPI_COMM_WORLD);
  } // for (int j ..
  int i1=_xyz->nX[0],i2=_xyz->nX[0]*_xyz->nX[1],jst;
  if (_xyz->ntype.compare("Moore")){
    int itmp[3]={-1,0,1},vneigh[26];
    for (int j=0;j<Na;++j){
      std::fill(vneigh,vneigh+26,-1);
      jn[2] = floor(vI[j]/i2);
      jn[1] = floor((vI[j]- i2*jn[2])/i1);
      jn[0] = vI[j] - i2*jn[2] - i1*jn[1];
      cc=0;
      for (int j3=0;j3<3;++j3){
	if ( (jn[2]+itmp[j3]<0) || (jn[2]+itmp[j3]>=_xyz->nX[2])){continue;}      
	for (int j2=0;j2<3;++j2){
	  if ( (jn[1]+itmp[j2]<0) || (jn[1]+itmp[j2]>=_xyz->nX[1])){continue;}      
	  for (int j1=0;j1<3;++j1){
	    if ( (jn[0]+itmp[j1]<0) || (jn[0]+itmp[j1]>=_xyz->nX[0])){continue;}      
	    jst = i2*(jn[2]+itmp[j3])+i1*(jn[1]+itmp[j2])+jn[0]+itmp[j1];
	    if (jst !=vI[j]){vneigh[cc]=jst; cc+=1;}
	  } // for (int j1...
	} // for (int j2...
      } // for (int j3...    
      for (int j1=0;j1<26;++j1){
	cc=std::distance(vI,std::find(vI,vI+Na,vneigh[j1]));
	if (cc <Na){
	  V[j][Nneigh[j]] = cc;
	  Nneigh[j]+=1;
	} // if (cc < ...
      } // for (int j1...
    } // for (int j=0...
  } else {
    int itmp[2]={-1,1},vneigh[6];
    for (int j=0;j<Na;++j){
      std::fill(vneigh,vneigh+6,-1);
      jn[2] = floor(vI[j]/i2);
      jn[1] = floor((vI[j]- i2*jn[2])/i1);
      jn[0] = vI[j] - i2*jn[2] - i1*jn[1];
      cc=0;
      for (int j1=0;j1<2;++j1){
	if ( (jn[2]+itmp[j1]<0) || (jn[2]+itmp[j1]>=_xyz->nX[2])){continue;}
	jst = i2*(jn[2]+itmp[j1])+i1*jn[1]+jn[0];
	vneigh[cc]=jst; 
	cc+=1;
      } // for (int j1...
      for (int j1=0;j1<2;++j1){
	if ( (jn[1]+itmp[j1]<0) || (jn[1]+itmp[j1]>=_xyz->nX[1])){continue;}
	jst = i2*jn[2]+i1*(jn[1]+itmp[j1])+jn[0];
	vneigh[cc]=jst; 
	cc+=1;
      } // for (int j1...
      for (int j1=0;j1<2;++j1){
	if ( (jn[0]+itmp[j1]<0) || (jn[0]+itmp[j1]>=_xyz->nX[0])){continue;}
	jst = i2*jn[2]+i1*jn[1]+jn[0]+itmp[j1];
	vneigh[cc]=jst; 
	cc+=1;
      } // for (int j1...
      for (int j1=0;j1<6;++j1){
	cc=std::distance(vI,std::find(vI,vI+Na,vneigh[j1]));
	if (cc <Na){
	  V[j][Nneigh[j]] = cc;
	  Nneigh[j]+=1;
	} // if (cc < ...
      } // for (int j1...
    } // for (int j=0...
  } // if (_xyz->ntype.compare...
  i1 = ceil( (double)Na / (double)_part->nprocs);
  i2 = floor( (double)Na / (double)i1);
  for (int j=0;j<(_part->nprocs+1);++j){
    if (j < (i2+1)){iv[j] = i1*j;}
    if (j>i2 && j<_part->nprocs){iv[j] = i2*i1 + 
	floor( (double)(Na-i2*i1)/(double)(_part->nprocs-i2));}
    if (j==_part->nprocs){iv[j]=Na;}
  } // for (int j...
  // end assemble arrays to be used to compute grain growth
  // capture all undercooled liquid voxels by growing grains
  int js,j1s,jx[3],jy[3],jTs[_part->nprocs], j1Ts[_part->nprocs],i3,countS=Na;
  double DtT[_part->nprocs],velX,velY,omega,ax[3],
    rRot[3][3],dlocX[3],locX[3],dr,vhat,tmn1[Na],
    timeUntil,dnx[3],tinc=0.0;
  std::fill(tmn1,tmn1+Na,0.0);
  i2 = _xyz->nX[0]; i3 = _xyz->nX[0]*_xyz->nX[1];
  cc=0;
  while (std::any_of(vS,vS+Na,[](int n){return n==1;}) && 
	 std::count(vS,vS+Na,2)!=countS && tinc<_temp->DelT) {
    double DtMin=1e6;
    countS = std::count(vS,vS+Na,2);
    for (int j=iv[_part->myid];j<iv[_part->myid+1];++j){
      if (vS[j]==1){
	jx[2] = floor(vI[j]/i3);
	jx[1] = floor((vI[j]- i3*jx[2])/i2);
	jx[0] = vI[j] - i3*jx[2] - i2*jx[1];
	velX=(5.51*pow(M_PI,2.0)*pow((- _xyz->mL)*(1-_xyz->kP),1.5)*
		 (_xyz->Gamma))*( pow((_xyz->tL - T[j]),2.5)/pow(_xyz->c0,1.5));
	for (int j1=0;j1<Nneigh[j];++j1){
	  if (vS[V[j][j1]]==1){continue;}
	  dnx[0] = (double(jx[0])+.5)*_xyz->dX[0] - CentroidA[3*V[j][j1]];
	  dnx[1] = (double(jx[1])+.5)*_xyz->dX[1] - CentroidA[3*V[j][j1]+1];
	  dnx[2] = (double(jx[2])+.5)*_xyz->dX[2] - CentroidA[3*V[j][j1]+2];
	  omega = cTheta[4*(G[V[j][j1]]-1)];
	  ax[0]=cTheta[4*(G[V[j][j1]]-1)+1];
	  ax[1]=cTheta[4*(G[V[j][j1]]-1)+2];
	  ax[2]=cTheta[4*(G[V[j][j1]]-1)+3];
	  // matrix is local->global; need to multiply by transpose for global->local            
	  rRot[0][0] = cos(omega) + pow(ax[0],2.0)*(1-cos(omega));
	  rRot[0][1] = ax[0]*ax[1]*(1-cos(omega)) - ax[2]*sin(omega);
	  rRot[0][2] = ax[0]*ax[2]*(1-cos(omega)) + ax[1]*sin(omega);
	  rRot[1][0] = ax[0]*ax[1]*(1-cos(omega)) + ax[2]*sin(omega);
	  rRot[1][1] = cos(omega) + pow(ax[1],2.0)*(1-cos(omega));
	  rRot[1][2] = ax[1]*ax[2]*(1-cos(omega)) - ax[0]*sin(omega);
	  rRot[2][0] = ax[2]*ax[0]*(1-cos(omega)) - ax[1]*sin(omega);
	  rRot[2][1] = ax[2]*ax[1]*(1-cos(omega)) + ax[0]*sin(omega);
	  rRot[2][2] = cos(omega) + pow(ax[2],2.0)*(1-cos(omega));
	  // put into 1st quadrant b/c of symmetry     
	  dlocX[0] = std::fabs(rRot[0][0]*dnx[0]+rRot[1][0]*dnx[1]+rRot[2][0]*dnx[2]);
	  dlocX[1] = std::fabs(rRot[0][1]*dnx[0]+rRot[1][1]*dnx[1]+rRot[2][1]*dnx[2]);
	  dlocX[2] = std::fabs(rRot[0][2]*dnx[0]+rRot[1][2]*dnx[1]+rRot[2][2]*dnx[2]);
	  dr = dlocX[0]+dlocX[1]+dlocX[2] - ExtA[V[j][j1]];
	  T[V[j][j1]]>=_xyz->tL ? velY=0.0: velY=(5.51*pow(M_PI,2.0)*pow((- _xyz->mL)*(1-_xyz->kP),1.5)*
		 (_xyz->Gamma))*( pow((_xyz->tL - T[V[j][j1]]),2.5)/pow(_xyz->c0,1.5));
	  //vhat = .5*(velX+velY);
	  vhat = velY;
	  //timeUntil = (dr-vhat*tmn1[V[j][j1]])/vhat;
	  timeUntil = dr/vhat;
	  if (timeUntil < DtMin){
	    DtMin = timeUntil;
	    js = j;
	    j1s = j1;
	  } // if (timeUntil < ...
	} // for (int j1
      } // if (vS[j...
    } // for (int j...
    MPI_Allgather(&DtMin,1,MPI_DOUBLE,&DtT[0],1,MPI_DOUBLE,MPI_COMM_WORLD);
    MPI_Allgather(&js,1,MPI_INT,&jTs[0],1,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather(&j1s,1,MPI_INT,&j1Ts[0],1,MPI_INT,MPI_COMM_WORLD);
    DtMin = *std::min_element(DtT,DtT+_part->nprocs);
    double DtMin2 = 1.01*DtMin;
    for (int j=0;j<NnucA;++j){
      if ( (tnucA[j] > tinc) && (tnucA[j] <= tinc+DtMin)){
	if (DtMin2> tnucA[j]-tinc){
	  DtMin2 = tnucA[j]-tinc;
	  js=nucA[j];
	} // if (DtMin2...
      } // if ( (tnucA ...
    } // for int j ...
    if (std::abs(DtMin2) < -1e6*std::abs(DtMin)){
      vS[js] = 2;
      nGrain+=1;
      G[js] = nGrain;
      sdloc= seed0 + 32*_xyz->tInd +64*nGrain;
      sa.GenerateSamples(1,sdloc,aa);
      cTheta.insert(cTheta.end(), &aa[0][0],&aa[0][0]+4);
      gNucleus.push_back(vI[js]);
      for (int j=0;j<Na;++j){tmn1[j]+=DtMin2;}
      tmn1[js] = 0.0;
      jx[2] = floor(vI[js]/i3);
      jx[1] = floor((vI[js]- i3*jx[2])/i2);
      jx[0] = vI[js] - i3*jx[2] - i2*jx[1];
      CentroidA[3*js] = (double(jx[0])+.5)*_xyz->dX[0];
      CentroidA[3*js+1] = (double(jx[1])+.5)*_xyz->dX[1];
      CentroidA[3*js+2] = (double(jx[2])+.5)*_xyz->dX[2];
      for (int j1=0;j1<Na;++j1){
	if (vS[j1]==2){
	  int vStmp[Nneigh[j1]];
	  for (int j2=0;j2<Nneigh[j1];++j2){vStmp[j2] = vS[V[j1][j2]];}
	  if (std::any_of(vStmp,vStmp+Nneigh[j1],[](int n){return n==1;})){
	    T[j1]>=_xyz->tL? velY=0.0: velY=(5.51*pow(M_PI,2.0)*pow((- _xyz->mL)*(1-_xyz->kP),1.5)*
				     (_xyz->Gamma))*( pow((_xyz->tL - T[j1]),2.5)/pow(_xyz->c0,1.5));
	    vhat = velY;
	    ExtA[j1]+=vhat*DtMin2;
	  } // if (std:: ...
	} // if (vS ...
      } // for (int j1...
      ExtA[js] = extentsInitialValue;
      std::cout << "nucleation"<<","<<cc<<","<<_part->myid<<std::endl;
    } else {
      double xI[3],xJ[3],d1I,dI2,d1J,dJ3,L12,L13,l,Lmud;
      int j;
      // compute new decentered octohedron return centroid and extents
      i1 = std::distance(DtT,std::min_element(DtT,DtT+_part->nprocs));
      js = jTs[i1];
      j1s = j1Ts[i1];
      if (DtMin<-10000.0){
	T[V[js][j1s]]>=_xyz->tL? velY=0.0: velY=(5.51*pow(M_PI,2.0)*pow((- _xyz->mL)*(1-_xyz->kP),1.5)*
					 (_xyz->Gamma))*( pow((_xyz->tL - T[V[js][j1s]]),2.5)/pow(_xyz->c0,1.5));
	vhat = velY;
	ExtA[V[js][j1s]] += vhat*DtMin;
      } else {
      for (int j1=0;j1<Na;++j1){
	if (vS[j1]==2){
	  int vStmp[Nneigh[j1]];
	  for (int j2=0;j2<Nneigh[j1];++j2){vStmp[j2] = vS[V[j1][j2]];}
	  if (std::any_of(vStmp,vStmp+Nneigh[j1],[](int n){return n==1;})){
	    T[j1]>=_xyz->tL? velY=0.0: velY=(5.51*pow(M_PI,2.0)*pow((- _xyz->mL)*(1-_xyz->kP),1.5)*
				   (_xyz->Gamma))*( pow((_xyz->tL - T[j1]),2.5)/pow(_xyz->c0,1.5));
	    vhat = velY;
	    ExtA[j1]+=vhat*DtMin;
	  } // if (std::...
	} // if (vS ...
      } // for (int j1...
      } // if (DtMin...
      vS[js] = 2;
      G[js] = G[V[js][j1s]];
      for (int j=0;j<Na;++j){tmn1[j]+=DtMin;}
      tmn1[js] = 0.0;
      jx[2] = floor(vI[js]/i3);
      jx[1] = floor((vI[js]- i3*jx[2])/i2);
      jx[0] = vI[js] - i3*jx[2] - i2*jx[1];
      jy[2] = floor(vI[V[js][j1s]]/i3);
      jy[1] =floor( (vI[V[js][j1s]]- i3*jy[2])/i2);
      jy[0] = vI[V[js][j1s]] -i3*jy[2] - i2*jy[1];
      l = pow( pow((jx[0]-jy[0])*_xyz->dX[0],2)+ pow((jx[1]-jy[1])*_xyz->dX[1],2)+
	       pow((jx[2]-jy[2])*_xyz->dX[2],2),.5)/3;
      dnx[0] = (double(jx[0])+.5)*_xyz->dX[0] - CentroidA[3*V[js][j1s]];
      dnx[1] = (double(jx[1])+.5)*_xyz->dX[1] - CentroidA[3*V[js][j1s]+1];
      dnx[2] = (double(jx[2])+.5)*_xyz->dX[2] - CentroidA[3*V[js][j1s]+2];
      omega = cTheta[4*(G[V[js][j1s]]-1)];
      ax[0]=cTheta[4*(G[V[js][j1s]]-1)+1];
      ax[1]=cTheta[4*(G[V[js][j1s]]-1)+2];
      ax[2]=cTheta[4*(G[V[js][j1s]]-1)+3];
      // matrix is local->global; need to multiply by transpose for global->local            
      rRot[0][0] = cos(omega) + pow(ax[0],2.0)*(1-cos(omega));
      rRot[0][1] = ax[0]*ax[1]*(1-cos(omega)) - ax[2]*sin(omega);
      rRot[0][2] = ax[0]*ax[2]*(1-cos(omega)) + ax[1]*sin(omega);
      rRot[1][0] = ax[0]*ax[1]*(1-cos(omega)) + ax[2]*sin(omega);
      rRot[1][1] = cos(omega) + pow(ax[1],2.0)*(1-cos(omega));
      rRot[1][2] = ax[1]*ax[2]*(1-cos(omega)) - ax[0]*sin(omega);
      rRot[2][0] = ax[2]*ax[0]*(1-cos(omega)) - ax[1]*sin(omega);
      rRot[2][1] = ax[2]*ax[1]*(1-cos(omega)) + ax[0]*sin(omega);
      rRot[2][2] = cos(omega) + pow(ax[2],2.0)*(1-cos(omega));
      locX[0] = rRot[0][0]*dnx[0]+rRot[1][0]*dnx[1]+rRot[2][0]*dnx[2];
      locX[1] = rRot[0][1]*dnx[0]+rRot[1][1]*dnx[1]+rRot[2][1]*dnx[2];
      locX[2] = rRot[0][2]*dnx[0]+rRot[1][2]*dnx[1]+rRot[2][2]*dnx[2];
      // signbit returns 0 if positive and 1 if negative
      dr = ExtA[V[js][j1s]];
      for (int j1=0;j1<6;++j1){sdiag[j1]={sdiag0[j1][0]*dr,sdiag0[j1][1]*dr,sdiag0[j1][2]*dr};}
      j = 4*std::signbit(locX[2])+ 2*std::signbit(locX[1])+ std::signbit(locX[0]);
      for (int j1=0;j1<3;++j1){ 
      ax[j1] = pow(locX[0]-sdiag[sInd[j][j1]][0],2.0)+
	pow(locX[1]-sdiag[sInd[j][j1]][1],2.0)+pow(locX[2]-sdiag[sInd[j][j1]][2],2.0);
      } // for (int j1...
      std::iota(jy,jy+3,0);
      std::sort(jy,jy+3,[&ax](int j1,int j2){return ax[j1]<ax[j2];});
      projectPointLine(locX,&sdiag[sInd[j][jy[0]]][0],&sdiag[sInd[j][jy[1]]][0],xI);
      projectPointLine(locX,&sdiag[sInd[j][jy[0]]][0],&sdiag[sInd[j][jy[2]]][0],xJ);
      d1I = pow(pow(sdiag[sInd[j][jy[0]]][0]-xI[0],2.0) + 
	pow(sdiag[sInd[j][jy[0]]][1]-xI[1],2.0) + 
	pow(sdiag[sInd[j][jy[0]]][2]-xI[2],2.0),.5);
      dI2 = pow(pow(sdiag[sInd[j][jy[1]]][0]-xI[0],2.0) + 
	pow(sdiag[sInd[j][jy[1]]][1]-xI[1],2.0) + 
	pow(sdiag[sInd[j][jy[1]]][2]-xI[2],2.0),.5);
      d1J = pow(pow(sdiag[sInd[j][jy[0]]][0]-xJ[0],2.0) + 
	pow(sdiag[sInd[j][jy[0]]][1]-xJ[1],2.0) + 
	pow(sdiag[sInd[j][jy[0]]][2]-xJ[2],2.0),.5);
      dJ3 = pow(pow(sdiag[sInd[j][jy[2]]][0]-xJ[0],2.0) + 
	pow(sdiag[sInd[j][jy[2]]][1]-xJ[1],2.0) + 
	pow(sdiag[sInd[j][jy[2]]][2]-xJ[2],2.0),.5);
      L12 = .5*(std::min(d1I,pow(3.0,.5)*l) + std::min(dI2,pow(3.0,.5)*l) );  
      L13 = .5*(std::min(d1J,pow(3.0,.5)*l) + std::min(dJ3,pow(3.0,.5)*l) );  
      Lmud =  pow(3.0,.5)* (pow(2.0/3.0,.5)*std::max(L12,L13));
      ExtA[js] = Lmud;
      dnx[0] = sdiag[sInd[j][jy[0]]][0] - Lmud*sdiag0[sInd[j][jy[0]]][0];
      dnx[1] = sdiag[sInd[j][jy[0]]][1] - Lmud*sdiag0[sInd[j][jy[0]]][1];
      dnx[2] = sdiag[sInd[j][jy[0]]][2] - Lmud*sdiag0[sInd[j][jy[0]]][2];
      locX[0] = rRot[0][0]*dnx[0]+rRot[0][1]*dnx[1]+rRot[0][2]*dnx[2];
      locX[1] = rRot[1][0]*dnx[0]+rRot[1][1]*dnx[1]+rRot[1][2]*dnx[2];
      locX[2] = rRot[2][0]*dnx[0]+rRot[2][1]*dnx[1]+rRot[2][2]*dnx[2];
      CentroidA[3*js] = CentroidA[3*V[js][j1s]] + locX[0];
      CentroidA[3*js+1] = CentroidA[3*V[js][j1s]+1] + locX[1];
      CentroidA[3*js+2] = CentroidA[3*V[js][j1s]+2] + locX[2];
      // end compute new decentered octohedron return centroid and extents
    } // if (DtMin2 < ...
    cc+=1;
    tinc += std::min(DtMin2,DtMin);
  } // while (std::any
  // end capture all undercooled liquid voxels by growing grains
  // bring global variables back to local variables
  cc1=0;
  for (int j=0;j<Ntot;++j){
    if (vState[j]==2 || (vState[j]==1 && _temp->TempCurr[j]< _xyz->tL)){
      vState[j] = vS[j0[_part->myid]+cc1];
      gID[j] = G[j0[_part->myid]+cc1];
      //extents[j] = ExtA[j0[_part->myid]+cc1];
      centroidOct[3*j] = CentroidA[3*j0[_part->myid]+3*cc1];
      centroidOct[3*j+1] = CentroidA[3*j0[_part->myid]+3*cc1+1];
      centroidOct[3*j+2] = CentroidA[3*j0[_part->myid]+3*cc1+2];
      cc1+=1;
    } // if (vState[j]...
  } // for (int j...
  // pass information 
  
  for (int j=Ntot;j<_part->nGhost+Ntot;++j){
    cc1=std::distance(vI,std::find(vI,vI+Na,_part->icellidLoc[j]));
    if (cc1<Na){
      vState[j] = vS[cc1];
      gID[j] = G[cc1];
    } // if (cc1<
  } // for (int j...
  /*
  _part->PassInformation(vState);
  _part->PassInformation(gID);
  */
  // mushy (vState=2) to solid (vState=3) if all neighbors 2 
  ConvertSolid1(0);
  _part->PassInformation(vState);
  //_part->PassInformation(extents);
}; // end UpdateVoxels4

void VoxelsCA::ConvertSolid(const int &iswitch)
{
  // if iswitch =0 then converts mushy to solid
  // else converts solid to mushy
  int Ntot = (_part->ncellLoc),cc,j2,iplay=_xyz->nX[0]*_xyz->nX[1]*_temp->ilaserLoc,
    iplay2=_xyz->nX[0]*_xyz->nX[1]*(NzhBP-1);
  std::vector<int> vneigh,neigh,i1(_part->ncellLoc,0);
  if (iswitch==0){
    // this checks if all neighbors are >=2 and converts to 3 if currently 2
    cc=0;
    for (int j=0; j < Ntot;++j){
      //if (_part->icellidLoc[j] >= iplay){continue;}
      if (vState[j]==2){
	_xyz->ComputeNeighborhood(_part->icellidLoc[j],_xyz->neighOrder,neigh);
	vneigh.assign(neigh.size(),0);
	for (int j1=0;j1<neigh.size();++j1){
	  j2 = std::distance(_part->icellidLoc.begin(),
			     std::find(_part->icellidLoc.begin(),_part->icellidLoc.end(),neigh[j1]));
	  vneigh[j1] = vState[j2];
	}
	if (std::all_of(vneigh.begin(),vneigh.end(),[](int n){return n>=2;})) {
	  i1[cc] = j;
	  cc +=1;
	} // if
      } // if (vState[j]==2)
    } // for j
    for (int j=0;j<cc;++j){
      vState[i1[j]] = 3;
      extents[i1[j]] = extentsInitialValue;
    } // for j
  } else {
    // this checks if any neighbors are 1 and converts to 2 if currently 3
    cc=0;
    for (int j=0; j < Ntot;++j){
      //if (_part->icellidLoc[j] >= iplay){continue;}
      if (_part->icellidLoc[j] < iplay2){continue;}
      if (vState[j]==3){
	_xyz->ComputeNeighborhood(_part->icellidLoc[j],_xyz->neighOrder,neigh);
	vneigh.assign(neigh.size(),0);
	for (int j1=0;j1<neigh.size();++j1){
	  j2 = std::distance(_part->icellidLoc.begin(),
			     std::find(_part->icellidLoc.begin(),_part->icellidLoc.end(),neigh[j1]));
	  vneigh[j1] = vState[j2];
	}
	if (std::any_of(vneigh.begin(),vneigh.end(),[](int n){return n==1;})) {
	  i1[cc] = j;
	  cc += 1;
	} // end if
      } // if (vState[j]==3)
    } // for j
    for (int j=0;j<cc;++j){
      vState[i1[j]] = 2;
    } // for j
  } // if (iswitch==0)
}; // ConvertSolid
void VoxelsCA::ConvertSolid1(const int &iswitch)
{
  // if iswitch =0 then converts mushy to solid
  // else converts solid to mushy
  int Ntot = (_part->ncellLoc),cc,j2,iplay=_xyz->nX[0]*_xyz->nX[1]*_temp->ilaserLoc,
    iplay2=_xyz->nX[0]*_xyz->nX[1]*(NzhBP-1);
  std::vector<int> vneigh,neigh,i1(_part->ncellLoc,0);
  if (iswitch==0){
    // this checks if all neighbors are >=2 and converts to 3 if currently 2
    cc=0;
    for (int j=0; j < Ntot;++j){
      //if (_part->icellidLoc[j] >= iplay){continue;}
      if (vState[j]==2){
	_xyz->ComputeNeighborhood(_part->icellidLoc[j],_xyz->neighOrder,neigh);
	vneigh.assign(neigh.size(),0);
	for (int j1=0;j1<neigh.size();++j1){
	  j2 = std::distance(_part->icellidLoc.begin(),
			     std::find(_part->icellidLoc.begin(),_part->icellidLoc.end(),neigh[j1]));
	  vneigh[j1] = vState[j2];
	}
	if (std::all_of(vneigh.begin(),vneigh.end(),[](int n){return n>=2;})) {
	  i1[cc] = j;
	  cc +=1;
	} // if
      } // if (vState[j]==2)
    } // for j
    for (int j=0;j<cc;++j){
      vState[i1[j]] = 3;
      extents[i1[j]] = extentsInitialValue;
    } // for j
  } else {
    // this checks if any neighbors are 1 and converts to 2 if currently 3
    cc=0;
    for (int j=0; j < Ntot;++j){
      if (vState[j]==3){
	_xyz->ComputeNeighborhood(_part->icellidLoc[j],_xyz->neighOrder,neigh);
	vneigh.assign(neigh.size(),0);
	for (int j1=0;j1<neigh.size();++j1){
	  j2 = std::distance(_part->icellidLoc.begin(),
			     std::find(_part->icellidLoc.begin(),_part->icellidLoc.end(),neigh[j1]));
	  vneigh[j1] = vState[j2];
	}
	if (std::any_of(vneigh.begin(),vneigh.end(),[](int n){return n==1;})) {
	  i1[cc] = j;
	  cc += 1;
	} // end if
      } // if (vState[j]==3)
    } // for j
    for (int j=0;j<cc;++j){
      vState[i1[j]] = 2;
    } // for j
  } // if (iswitch==0)
}; // ConvertSolid1

void VoxelsCA::ComputeVoxelCapture(){
  //determines if liquid is captured by growing mushy neighbor
  int j,j1,j2,j3,j1n,j2n,j3n,jck,Ntot,i1cc,jneigh,iplay,Ntot2;
  double x0,y0,z0,xj,yj,zj,omega,d;
  std::vector<int> neigh,jcklist,ichk;
  std::vector<double> ax(3),avec(3),planeEval(8);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::vector<std::vector<double>> nA(8,std::vector<double>(3)),
    pA(8,std::vector<double>(3)),rRot(3,std::vector<double>(3));

  nA[0] = {1.0/pow(3.0,.5),1.0/pow(3.0,.5),1.0/pow(3.0,.5)};
  nA[1] = {-1.0/pow(3.0,.5),1.0/pow(3.0,.5),1.0/pow(3.0,.5)};
  nA[2] = {1.0/pow(3.0,.5),-1.0/pow(3.0,.5),1.0/pow(3.0,.5)};
  nA[3] = {-1.0/pow(3.0,.5),-1.0/pow(3.0,.5),1.0/pow(3.0,.5)};
  nA[4] = {1.0/pow(3.0,.5),1.0/pow(3.0,.5),-1.0/pow(3.0,.5)};
  nA[5] = {-1.0/pow(3.0,.5),1.0/pow(3.0,.5),-1.0/pow(3.0,.5)};
  nA[6] = {1.0/pow(3.0,.5),-1.0/pow(3.0,.5),-1.0/pow(3.0,.5)};
  nA[7] = {-1.0/pow(3.0,.5),-1.0/pow(3.0,.5),-1.0/pow(3.0,.5)};
  pA[0] = {1.0,0.0,0.0};
  pA[1] = {-1.0,0.0,0.0};
  pA[2] = {0.0,-1.0,0.0};
  pA[3] = {0.0,-1.0,0.0};
  pA[4] = {0.0,0.0,-1.0};
  pA[5] = {0.0,0.0,-1.0};
  pA[6] = {1.0,0.0,0.0};
  pA[7] = {0.0,-1.0,0.0};

  Ntot = _part->ncellLoc;
  Ntot2 = _part->ncellLoc + _part->nGhost;
  ichk.assign(Ntot2,0);
  iplay =_xyz->nX[0]*_xyz->nX[1]*_temp->ilaserLoc;
  for (int j=0;j<Ntot;++j){
    if (_part->icellidLoc[j] >= iplay){continue;}
    if (vState[j]==1 && _temp->TempCurr[j]<_xyz->tL){
      _xyz->ComputeNeighborhood(_part->icellidLoc[j],_xyz->neighOrder,neigh);
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x0 = (double(j1)+.5)*(_xyz->dX[0]);
      y0 = (double(j2)+.5)*(_xyz->dX[1]);
      z0 = (double(j3)+.5)*(_xyz->dX[2]);
      jcklist.assign(neigh.size(),0);
      i1cc=0;
      for (int jj=0;jj<neigh.size();++jj){
	jneigh = std::distance(_part->icellidLoc.begin(),
			   std::find(_part->icellidLoc.begin(),_part->icellidLoc.end(),neigh[jj]));
	if (vState[jneigh]==2){
	  j3n = floor(gNucleus[gID[jneigh]-1]/(_xyz->nX[0]*_xyz->nX[1]));
	  j2n = floor( (gNucleus[gID[jneigh]-1]-_xyz->nX[0]*_xyz->nX[1]*j3n)  /_xyz->nX[0]);
	  j1n = gNucleus[gID[jneigh]-1] - _xyz->nX[0]*_xyz->nX[1]*j3n - _xyz->nX[0]*j2n;
	  xj = (double(j1n)+.5)*(_xyz->dX[0]);
	  yj = (double(j2n)+.5)*(_xyz->dX[1]);
	  zj = (double(j3n)+.5)*(_xyz->dX[2]);
	  omega = cTheta[4*(gID[jneigh]-1)];
	  ax = {cTheta[4*(gID[jneigh]-1)+1],cTheta[4*(gID[jneigh]-1)+2], cTheta[4*(gID[jneigh]-1)+3]};
	  rRot[0][0] = cos(omega) + pow(ax[0],2.0)*(1-cos(omega));
	  rRot[0][1] = ax[0]*ax[1]*(1-cos(omega)) - ax[2]*sin(omega);
	  rRot[0][2] = ax[0]*ax[2]*(1-cos(omega)) + ax[1]*sin(omega);
	  rRot[1][0] = ax[0]*ax[1]*(1-cos(omega)) + ax[2]*sin(omega);
	  rRot[1][1] = cos(omega) + pow(ax[1],2.0)*(1-cos(omega));
	  rRot[1][2] = ax[1]*ax[2]*(1-cos(omega)) - ax[0]*sin(omega);
	  rRot[2][0] = ax[2]*ax[0]*(1-cos(omega)) - ax[1]*sin(omega);
	  rRot[2][1] = ax[2]*ax[1]*(1-cos(omega)) + ax[0]*sin(omega);
	  rRot[2][2] = cos(omega) + pow(ax[2],2.0)*(1-cos(omega));
	  for (int j2=0;j2<8;++j2){
	    avec[0] = rRot[0][0]*nA[j2][0]+rRot[0][1]*nA[j2][1]+rRot[0][2]*nA[j2][2];
	    avec[1] = rRot[1][0]*nA[j2][0]+rRot[1][1]*nA[j2][1]+rRot[1][2]*nA[j2][2];
	    avec[2] = rRot[2][0]*nA[j2][0]+rRot[2][1]*nA[j2][1]+rRot[2][2]*nA[j2][2];
	    d = - (xj*avec[0] + yj*avec[1] + zj*avec[2] + 
		   (nA[j2][0]*pA[j2][0]+nA[j2][1]*pA[j2][1]+nA[j2][2]*pA[j2][2])*extents[jneigh]);
	    planeEval[j2] = avec[0]*x0+avec[1]*y0+avec[2]*z0 + d;
	  } // for (int j2....
	  if (std::all_of(planeEval.begin(),planeEval.end(),[](double x){return x<0;})){
	    if (ichk[jneigh]==0){
	      jcklist[i1cc] = jj;
	      i1cc +=1;
	    }
	  }       
	} // if (vState[jneigh]==2)
      } // end for jj
      if (i1cc > 0){
	std::uniform_int_distribution<> ir(0,i1cc-1);
	std::default_random_engine gen1(seed0+10*_part->icellidLoc[j]+15*_xyz->tInd);
     	vState[j]=2;
	jck = ir(gen);
	ichk[j] = 1;
	jneigh = std::distance(_part->icellidLoc.begin(),
	       std::find(_part->icellidLoc.begin(),_part->icellidLoc.end(),neigh[jcklist[jck]]));
	gID[j] = gID[jneigh];
	extents[j] = extents[jneigh];
      } // if (i1cc>0
    } // if vState[j]==1
  } // for j
}; // end ComputeVoxelCapture

void VoxelsCA::ComputeVoxelCapture2(){
  int Ntot = _part->ncellLoc, j,j1,jx1,jx2,jx3,jy1,jy2,jy3,ncomp,cc,ncompT;
  double dmax,dmin,omega,vhat,dlocnorm, rRot[3][3],dloc[3],ax[3],delxy[3],velCurr;
  // assemble voxel index for all vState=2 over all processes and give to all
  ncomp=0;
  for (int j1=0;j1<Ntot;++j1){if (vState[j1]==2){ncomp+=1; } }
  int ncompv[_part->nprocs];
  MPI_Allgather(&ncomp,1,MPI_INT,&ncompv[0],1,MPI_INT,MPI_COMM_WORLD);
  ncompT=0;
  for (int j1=0;j1<_part->nprocs;++j1){ncompT+=ncompv[j1];}
  int icid0[ncomp],*is=icid0, icidT[ncompT],gid0[ncomp],*gs=gid0,gidT[ncompT];
  double vel0[ncomp],*vs=vel0,velT[ncompT];
  for (int j1=0;j1<Ntot;++j1){
    if (vState[j1]==2){
      *(is++) = _part->icellidLoc[j1];
      *(gs++) = gID[j1];
      _xyz->tL > _temp->TempCurr[j1] ?
	*(vs++)=(5.51*pow(M_PI,2.0)*pow((- _xyz->mL)*(1-_xyz->kP),1.5)*
       (_xyz->Gamma))*( pow((_xyz->tL - _temp->TempCurr[j1]),2.5)/
			pow(_xyz->c0,1.5)): *(vs++)=0.0;
    }
  }
  cc=0;
  for (int j1=0;j1<_part->nprocs;++j1){
    if (ncompv[j1]>0){
      if (_part->myid==j1){
	for (int j2=0;j2<ncompv[j1];++j2){
	  icidT[cc+j2] = icid0[j2];
	  gidT[cc+j2] = gid0[j2];
	  velT[cc+j2] = vel0[j2];
	}
      }
      MPI_Bcast(&icidT[cc],ncompv[j1],MPI_INT,j1,MPI_COMM_WORLD);
      MPI_Bcast(&gidT[cc],ncompv[j1],MPI_INT,j1,MPI_COMM_WORLD);
      MPI_Bcast(&velT[cc],ncompv[j1],MPI_DOUBLE,j1,MPI_COMM_WORLD);
      cc+=ncompv[j1];
    }
  }
  // end assemble voxel index for all vState=2 over all processes and give to all
  for (int j=0;j<Ntot;++j){
    if (vState[j]==1 && _temp->TempCurr[j]< _xyz->tL){
      jx3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      jx2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*jx3)/_xyz->nX[0]);
      jx1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*jx3 - _xyz->nX[0]*jx2;
      velCurr=(5.51*pow(M_PI,2.0)*pow((- _xyz->mL)*(1-_xyz->kP),1.5)*
       (_xyz->Gamma))*( pow((_xyz->tL - _temp->TempCurr[j]),2.5)/
			pow(_xyz->c0,1.5));
      double A[4][ncompT],p[ncompT],pg[nGrain],
	*ps=pg,dminT=_xyz->dX[0]*1000; // container for distribution of  Delta t (Dt_min,Dt_max, A,B)
      std::fill(p,p+ncompT,1); // initialize all to value 1
      std::fill(pg,pg+nGrain,0.0); // initialize all to value 0
      for (int j1=0;j1<ncompT;++j1){
	jy3 = floor(icidT[j1]/(_xyz->nX[0]*_xyz->nX[1]));
	jy2 = floor( (icidT[j1]- _xyz->nX[0]*_xyz->nX[1]*jy3)/_xyz->nX[0]);
	jy1 = icidT[j1] - _xyz->nX[0]*_xyz->nX[1]*jy3 - _xyz->nX[0]*jy2;
	delxy[0]=(jy1-jx1)*_xyz->dX[0]; delxy[1]=(jy2-jx2)*_xyz->dX[1]; delxy[2]=(jy3-jx3)*_xyz->dX[2];
	dmin = pow(pow(delxy[0],2.0)+ pow(delxy[1],2.0)+pow(delxy[2],2.0),.5);
	dminT = std::min(dmin,dminT);
	dmax = dmin*M_PI;
	omega = cTheta[4*(gidT[j1]-1)];
	ax[0]=cTheta[4*(gidT[j1]-1)+1];ax[1]=cTheta[4*(gidT[j1]-1)+2];ax[2]=cTheta[4*(gidT[j1]-1)+3];
	// matrix is local->global; need to multiply by transpose for global->local
	rRot[0][0] = cos(omega) + pow(ax[0],2.0)*(1-cos(omega));
	rRot[0][1] = ax[0]*ax[1]*(1-cos(omega)) - ax[2]*sin(omega);
	rRot[0][2] = ax[0]*ax[2]*(1-cos(omega)) + ax[1]*sin(omega);
	rRot[1][0] = ax[0]*ax[1]*(1-cos(omega)) + ax[2]*sin(omega);
	rRot[1][1] = cos(omega) + pow(ax[1],2.0)*(1-cos(omega));
	rRot[1][2] = ax[1]*ax[2]*(1-cos(omega)) - ax[0]*sin(omega);
	rRot[2][0] = ax[2]*ax[0]*(1-cos(omega)) - ax[1]*sin(omega);
	rRot[2][1] = ax[2]*ax[1]*(1-cos(omega)) + ax[0]*sin(omega);
	rRot[2][2] = cos(omega) + pow(ax[2],2.0)*(1-cos(omega));
	// put into 1st quadrant b/c of symmetry
	dloc[0] = std::fabs(rRot[0][0]*delxy[0]+rRot[1][0]*delxy[1]+rRot[2][0]*delxy[2]);
	dloc[1] = std::fabs(rRot[0][1]*delxy[0]+rRot[1][1]*delxy[1]+rRot[2][1]*delxy[2]);
	dloc[2] = std::fabs(rRot[0][2]*delxy[0]+rRot[1][2]*delxy[1]+rRot[2][2]*delxy[2]);
	dlocnorm = pow(pow(dloc[0],2.0)+pow(dloc[1],2.0)+pow(dloc[2],2.0),.5);
	vhat = .5*(velT[j1]+velCurr)*dlocnorm/(dloc[0]+dloc[1]+dloc[2]);
	A[0][j1] = dmin/vhat * 1e6;
	A[1][j1] = dmax/vhat * 1e6;
	A[2][j1] = log(_xyz->ethresh)*pow(_xyz->ethresh,A[0][j1]/(A[0][j1]-A[1][j1]))/(A[0][j1]-A[1][j1]);
	A[3][j1] = log(_xyz->ethresh)/(A[0][j1]-A[1][j1]);
      } // for (int j1=0...)
      for (int j1=0;j1<ncompT;++j1){
	for (int j2=0;j2<ncompT;++j2){
	  if (j1!=j2){
	    p[j1]*=A[2][j1]*A[2][j2]/A[3][j1]/A[3][j2]*
	      (exp(-A[3][j1]*A[0][j1]-A[3][j2]*std::max(A[0][j1],A[0][j2]))) - 
	      A[2][j1]*A[2][j2]/A[3][j1]/(A[3][j1]+A[3][j2])*
	      exp(-(A[3][j1]+A[3][j2])*std::max(A[0][j1],A[0][j2]));
	  } // if (j1!=j2
	} // for (int j2...
      } // for (int j1
      for (int j1=0;j1<ncompT;++j1){pg[gidT[j1]-1]+= p[j1];}
      std::default_random_engine gen1(seed0+10*_part->icellidLoc[j]+15*_xyz->tInd);
      std::discrete_distribution<int> pind(ps,ps+nGrain);
      std::vector<double> pmf=pind.probabilities();
      double maxP = *std::max_element(pmf.begin(),pmf.end());
      if (_temp->tInd==3){std::cout << std::distance(pmf.begin(),std::max_element(pmf.begin(),pmf.end()))<<","<<maxP<<std::endl;}
      if (dminT <= 2*_xyz->dX[0]){
	//int ii=pind(gen1)+1;
	gID[j] = pind(gen1)+1;
	//std::cout << "1,"<< ii<<std::endl;
	vState[j] = 2;
      } else if (dminT > 2*_xyz->dX[0] && maxP > _xyz->deltaThresh*5){
	//int ii=pind(gen1)+1;
	gID[j] = pind(gen1)+1;
	vState[j] = 2;
	//std::cout << "2,"<< ii<<std::endl;
      }
    } // if (vState[j]==1...
  } // for (int j ...

}; // end ComputeVoxelCapture2

void VoxelsCA::SetLiquid(){
  // makes cell liquid if temperature exceeds liquidus
  int Ntot = _part->ncellLoc + _part->nGhost; 
  int n1 = NzhBP*_xyz->nX[0]*_xyz->nX[1];
  int j;
  for (int j=0;j<Ntot;++j){
    if (_part->icellidLoc[j] <n1){continue;}
    if (_temp->TempCurr[j] >= _xyz->tL ) { 
      vState[j] = 1;
      gID[j] = 0; // flag that it loses its grain
      extents[j] = extentsInitialValue;
    } 
  } // end for j
};

void VoxelsCA::SetLiquid2(){
  // makes cell liquid if temperature exceeds liquidus
  int Ntot = _part->ncellLoc + _part->nGhost; 
  int n1 = NzhBP*_xyz->nX[0]*_xyz->nX[1];
  int j,j3;
  for (int j=0;j<Ntot;++j){
    if (_part->icellidLoc[j] <n1){
      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      if (j3==NzhBP-1 && _temp->TempCurr[j]> _xyz->tL*1.0){vState[j] = 2;}
    } else{
      if (_temp->TempCurr[j] >= _xyz->tL ) {
	vState[j] = 1;
	gID[j] = 0; // flag that it loses its grain
	extents[j] = extentsInitialValue;
      }
    } 
  } // end for j
};

void VoxelsCA::SetLiquid3(){
  // makes cell liquid if temperature exceeds liquidus
  int Ntot = _part->ncellLoc + _part->nGhost,n1; 
  int j;
  n1 = _xyz->nX[0]*_xyz->nX[1];
  for (int j=0;j<Ntot;++j){
    if (_temp->TempCurr[j] >= _xyz->tL ) { 
      if (_part->icellidLoc[j]<n1){
	vState[j]=2;
      } else{
	vState[j] = 1;
	gID[j] = 0; // flag that it loses its grain
	//extents[j] = extentsInitialValue;
      } // if (j<n1
    } 
  } // end for j
};

void VoxelsCA::ExtentsInitialize(){
  int Ntot = _part->ncellLoc + _part->nGhost;
  int j,j1,j2,j3,j1n,j2n,j3n;
  double xj,yj,zj,xn,yn,zn;
  for (int j=0;j<Ntot;++j){
    if (vState[j]==2){
      j3n = floor(gNucleus[gID[j]-1]/(_xyz->nX[0]*_xyz->nX[1]));
      j2n = floor( (gNucleus[gID[j]-1]-_xyz->nX[0]*_xyz->nX[1]*j3n)  /_xyz->nX[0]);
      j1n = gNucleus[gID[j]-1] - _xyz->nX[0]*_xyz->nX[1]*j3n - _xyz->nX[0]*j2n;
      xn = (double(j1n)+.5)*(_xyz->dX[0]);
      yn = (double(j2n)+.5)*(_xyz->dX[1]);
      zn = (double(j3n)+.5)*(_xyz->dX[2]);

      j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      xj = (double(j1)+.5)*(_xyz->dX[0]);
      yj = (double(j2)+.5)*(_xyz->dX[1]);
      zj = (double(j3)+.5)*(_xyz->dX[2]);
      extents[j] = pow( pow(xj-xn,2.0)+pow(yj-yn,2.0)+pow(zj-zn,2.0),.5);
    } else {
	extents[j] = extentsInitialValue;
    } // if/else (vState[j]==2...
  } // for (int j=0...
}; // end ExtentsInitialize

void VoxelsCA::ZeroVoxels(){
  // sets voxels to zero where laser hasn't reached yet
  int j,j1,j2,j3;
  int Ntot = _part->ncellLoc + _part->nGhost,n1; 
  double x,y;
  if (_temp->patternID==0){
    x = fmod(_temp->tInd,(_temp->nTTemp[0]))*_temp->bmDX[0] - _temp->offset[0];
    y = floor(fmod(_temp->tInd,(_temp->nTTemp[0]*_temp->nTTemp[1]))/_temp->nTTemp[0])*_temp->bmDX[1]-_temp->offset[1];
    j1 = int(floor(std::max(x,0.0)/_xyz->dX[0]));
    j2 = int(floor(std::max(y,0.0)/_xyz->dX[1]));
    j3 = _temp->ilaserLoc;
    n1 = _xyz->nX[0]*_xyz->nX[1]*j3 + _xyz->nX[0]*j2 + j1;
    for (int j=0;j<Ntot;++j){
      if (_part->icellidLoc[j] <n1){continue;}
      vState[j] = 0;
      gID[j] = 0;
    } // for j
  } // if (_temp->patternID==0
};
void VoxelsCA::ZeroVoxels1(){
  // sets voxels to zero where laser hasn't reached yet
  int j,j1b,j2b,j3b,j2,j3;
  int Ntot = _part->ncellLoc + _part->nGhost,n1; 
  double x,y;
  if (_temp->patternID==0){
    n1 = _xyz->nX[0]*_xyz->nX[1]*_temp->ilaserLoc;
    for (int j=0;j<Ntot;++j){
      if (_part->icellidLoc[j] >=n1){
	vState[j] = 0;
	gID[j] = 0;
      } /*else {
	j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
	j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
	x = fmod(_temp->tInd,(_temp->nTTemp[0]))*_temp->bmDX[0] - _temp->offset[0];
	y = floor(fmod(_temp->tInd,(_temp->nTTemp[0]*_temp->nTTemp[1]))/_temp->nTTemp[0])*_temp->bmDX[1]-_temp->offset[1];
	j2b = int(floor(std::max(y,0.0)/_xyz->dX[1]));
	j3b = NzhBP + floor( (floor( (_temp->tInd)/(_temp->nTTemp[0]*_temp->nTTemp[1])))*_xyz->layerT/_xyz->dX[2]);
	if (j3 >= j3b && j2 > j2b){
	  vState[j] = 0;
	  gID[j] = 0;
	} // (if j3 > ...
      } // if (_part->icellidLoc ...
      */
    } // for j
  } // if (_temp->patternID==0
}; // ZeroVoxels1

void VoxelsCA::CheckTimeSkip(){
  // checks if any voxel with vState=1 has temperature < melting; if none, then 
  // skip ahead to next time increment for temperature
  int Ntot = _part->ncellLoc + _part->nGhost,icheck=0,icheckT;
  for (int j=0;j<Ntot;++j){
    if (vState[j]==1 && _temp->TempCurr[j]<_xyz->tL){icheck=1;}
  }
  MPI_Allreduce(&icheck,&icheckT,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  if (icheckT==0){_xyz->SkipTime(_temp->DelT); if (_part->myid==0){std::cout << _temp->tInd<<","<<_xyz->tInd<<std::endl;}}

}; // end CheckTimeSkip

void VoxelsCA::ComputeExtents(){
  // compute velocity and update extent of [1,0] direction if vState=2
  int j,Ntot,jneigh;
  double etaScale;
  Ntot = _part->ncellLoc;
  std::vector<double> velocity(Ntot,0);
  std::vector<int> neigh,vneigh;
  std::string strNeigh;
  vmax = 0.0;
  strNeigh ="first";
  for (int j=0;j<Ntot;++j){
    if (vState[j] == 2){
      _xyz->ComputeNeighborhood(_part->icellidLoc[j],strNeigh,neigh);
      vneigh.assign(neigh.size(),0);
      for (int jj=0;jj<neigh.size();++jj){
	jneigh = std::distance(_part->icellidLoc.begin(),
			     std::find(_part->icellidLoc.begin(),_part->icellidLoc.end(),neigh[jj]));
	vneigh[jj]=vState[jneigh];
      } // for (int jj ...
      if (_temp->TempCurr[j]>=_xyz->tL){
	velocity[j] = 0.0;
	//extents[j] = extentsInitialValue;
      } else if (!std::any_of(vneigh.begin(),vneigh.end(),[](int n){return n==1;})) {
	velocity[j] = 0.0;
      } else{
	velocity[j] = (_xyz->dL)/
	  (5.51*pow(M_PI,2.0)*pow((- _xyz->mL)*(1-_xyz->kP),1.5)*
	   (_xyz->Gamma))*( pow((_xyz->tL - _temp->TempCurr[j]),2.5)/
			  pow(_xyz->c0,1.5));
       }	
      vmax = std::max(vmax,velocity[j]);
    } // if (vState[j]==2)
  } // for j  
  double vmax2 = vmax;
  MPI_Allreduce(&vmax2,&vmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  _xyz->UpdateTimeTest1(vmax);

  /*
  etaScale = .5*_xyz->deltaXmin/_temp->DelT/vmax*20;
  for (int j=0;j<Ntot;++j){
    if (vState[j] == 2){velocity[j]*=etaScale;}
  }
  _xyz->UpdateTime2(_temp->DelT/20.0);
  */
  
  for (int j=0;j<Ntot;++j){
    if (vState[j] == 2){
      extents[j] +=velocity[j]*(_xyz->dt);
      
    }
  }
};

void VoxelsCA::ComputeNucleation1(){
  int Ntot = _part->ncellLoc, nmushy,itmp,ngtmp=0,ngnew,cc,seed,cc1;
  std::vector<double> nrate(Ntot,0),gidtmp,cthtmp,tmpall1,tmp1;
  std::vector<std::vector<double>> aa;
  std::vector<int> ngvec(_part->nprocs,0),ind,tmp2,tmpall2,gnuctmp;
  unsigned int sdloc;
  std::srand(std::time(nullptr));
  std::default_random_engine g(seed1+21*_xyz->tInd);
  std::uniform_real_distribution<double> xrand(0,1);

  for (int j=0;j<Ntot;++j){
    nrate[j] = -2*_xyz->muN *(_xyz->tL - _temp->TempCurr[j])*_temp->DDtTemp[j];
  } // for j
  itmp = 0;
  for (int j=0;j<Ntot;++j){
    if (vState[j] == 2){itmp +=1;}
  }
  MPI_Allreduce(&itmp,&nmushy,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  sdloc = (unsigned( double(g())/double(g.max())*pow(2.0,32.0)));
  for (int j=0;j<Ntot;++j){
    if (vState[j] == 2){
      std::default_random_engine g1(seed1+11*_part->icellidLoc[j] + 21*_xyz->tInd);
      if (xrand(g1) < nrate[j]/nmushy*_xyz->dt){
	ind.push_back(j);
	ngtmp += 1;
	gnuctmp.push_back(_part->icellidLoc[j]);
	//	sdloc.push_back(int(g1()/1000 + 1234));
      } // if (xrand(g) < ...
    } // if (vState[j]==2
  } // for (int j ...
  aa.assign(ngtmp,std::vector<double>(4));
  SampleOrientation sa;
  sa.GenerateSamples(ngtmp,sdloc,aa);
  for (int j1=0;j1<ngtmp;++j1){
    cthtmp.push_back(aa[j1][0]);
    cthtmp.push_back(aa[j1][1]);
    cthtmp.push_back(aa[j1][2]);
    cthtmp.push_back(aa[j1][3]);
  } // end for j1   

  MPI_Allgather(&ngtmp,1,MPI_INT,&ngvec[0],1,MPI_INT,MPI_COMM_WORLD);
  ngnew = std::accumulate(ngvec.begin(),ngvec.end(),0);
  if (ngnew > 0){
    tmpall1.assign(4*ngnew,0);
    tmpall2.assign(ngnew,0);
    cc=0;
    cc1=0;
    for (int j=0;j<_part->nprocs;++j){
      tmp1.assign(4*ngvec[j]+1,0);
      tmp2.assign(ngvec[j]+1,0);
      if (_part->myid==j){
	for (int j1=0;j1<ind.size();++j1){
	  tmp1[4*j1+1] = cthtmp[4*j1];
	  tmp1[4*j1+1+1] = cthtmp[4*j1+1];
	  tmp1[4*j1+2+1] = cthtmp[4*j1+2];
	  tmp1[4*j1+3+1] = cthtmp[4*j1+3];
	  tmp2[j1+1] = gnuctmp[j1];
	}
      }
      MPI_Bcast(&tmp1[0],4*ngvec[j]+1,MPI_DOUBLE,j,MPI_COMM_WORLD);
      MPI_Bcast(&tmp2[0],ngvec[j]+1,MPI_INT,j,MPI_COMM_WORLD);
      for (int j1=0;j1<ngvec[j];++j1){
	tmpall1[cc1+4*j1] = tmp1[4*j1+1];
	tmpall1[cc1+4*j1+1] = tmp1[4*j1+1+1];
	tmpall1[cc1+4*j1+2] = tmp1[4*j1+2+1];
	tmpall1[cc1+4*j1+3] = tmp1[4*j1+3+1];
	tmpall2[cc+j1] = tmp2[j1+1];
      } // for (int j1...
      cc +=ngvec[j];
      cc1 +=4*ngvec[j];
    } // for (int j ...
    for (int j=0; j<4*ngnew;++j){cTheta.push_back(tmpall1[j]);}
    for (int j=0; j<ngnew;++j){gNucleus.push_back(tmpall2[j]);}
    cc = 1;
    for (int j=0;j<_part->nprocs;++j){    
      if (_part->myid ==j){
	for (int j1=0;j1<ind.size();++j1){
	  gID[ind[j1]] = nGrain + cc + j1;
	  extents[ind[j1]] = extentsInitialValue2;
	} // for (int j1 ...
      } // if (_part->myid ...
      cc +=ngvec[j];
    } // for (int j ...
    nGrain +=ngnew;
    _part->PassInformation(gID);
    _part->PassInformation(extents);

  } // if (ngnew > 0)  
};

void VoxelsCA::WriteToVTU1(const std::string &filename)
{
  // writes gID, vState, cTheta per voxel
  std::string vtkFilename = filename + ".vtu";
  std::ofstream fp;
  int Offset=0, nC=_part->ncellLoc, nP=_part->npointLoc,j1,j2,j3,cell_offsets;
  unsigned char cell_type;
  std::vector< float> TempOut(nC,0);
  float IPFmapBD[3*nC];
  double vBD[3]={0.0,0.0,1.0},omega,ax[3],vCD[3],mxAng,blue,green,red,rRot[3][3];
  for (int j=0;j<nC;++j){
    TempOut[j] = _temp->TempCurr[j];
    if (gID[j]<1){
      IPFmapBD[3*j] = 0.0;
      IPFmapBD[3*j+1] = 0.0;
      IPFmapBD[3*j+2] = 0.0;
    } else {
      omega = cTheta[4*(gID[j]-1)];
      ax[0]= cTheta[4*(gID[j]-1)+1];
      ax[1]= cTheta[4*(gID[j]-1)+2];
      ax[2]= cTheta[4*(gID[j]-1)+3];
      // matrix is local->global; need to multiply by transpose for global->local            
      rRot[0][0] = cos(omega) + pow(ax[0],2.0)*(1-cos(omega));
      rRot[0][1] = ax[0]*ax[1]*(1-cos(omega)) - ax[2]*sin(omega);
      rRot[0][2] = ax[0]*ax[2]*(1-cos(omega)) + ax[1]*sin(omega);
      rRot[1][0] = ax[0]*ax[1]*(1-cos(omega)) + ax[2]*sin(omega);
      rRot[1][1] = cos(omega) + pow(ax[1],2.0)*(1-cos(omega));
      rRot[1][2] = ax[1]*ax[2]*(1-cos(omega)) - ax[0]*sin(omega);
      rRot[2][0] = ax[2]*ax[0]*(1-cos(omega)) - ax[1]*sin(omega);
      rRot[2][1] = ax[2]*ax[1]*(1-cos(omega)) + ax[0]*sin(omega);
      rRot[2][2] = cos(omega) + pow(ax[2],2.0)*(1-cos(omega));
      vCD[0] = std::fabs(rRot[0][0]*vBD[0]+rRot[1][0]*vBD[1]+rRot[2][0]*vBD[2]);
      vCD[1] = std::fabs(rRot[0][1]*vBD[0]+rRot[1][1]*vBD[1]+rRot[2][1]*vBD[2]);
      vCD[2] = std::fabs(rRot[0][2]*vBD[0]+rRot[1][2]*vBD[1]+rRot[2][2]*vBD[2]);
      std::sort(vCD,vCD+3);
      std::swap(vCD[0],vCD[1]);
      mxAng = M_PI/4.0;
      red = std::fabs( (mxAng - acos(vCD[2]))/mxAng );
      blue = atan2(vCD[1],vCD[0]);
      green = mxAng - blue;
      blue *= (1-red)/mxAng;
      green *= (1-red)/mxAng;
      IPFmapBD[3*j] = red;
      IPFmapBD[3*j+1] = green;
      IPFmapBD[3*j+2] = blue;
    }
  }
  if (_xyz->nnodePerCell==4){cell_type=9;}
  if (_xyz->nnodePerCell==8){cell_type=12;}
  if (_part->myid == 0) {
    fp.open(vtkFilename.c_str());
    if (!fp) throw std::runtime_error("Cannot create VTU file: " + filename);
    fp << "<?xml version='1.0'?>" << std::endl;
    fp << "<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'>" << std::endl;
    fp << "  <UnstructuredGrid>" << std::endl;
    fp << "    <FieldData>" << std::endl;
    fp << "      <DataArray type='Int32' Name='time_step' NumberOfTuples='1' format='ascii'>" << std::endl;
    fp << "        " << _xyz->tInd << std::endl;
    fp << "      </DataArray>" << std::endl;
    fp << "    </FieldData>" << std::endl;
    fp.close();
  } // if (myid==0)
  for (int rank=0;rank<_part->nprocs;++rank){
    int OffsetP=Offset;
    MPI_Allreduce(&OffsetP,&Offset,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (_part->myid != rank) { Offset=0; continue; }
    fp.open(vtkFilename.c_str(), std::fstream::app);
    fp << "    <Piece NumberOfPoints='" << nP << "' NumberOfCells='" << nC << "'>" << std::endl;
    fp << "      <Points>" << std::endl;
    fp << "        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='appended' offset='" << Offset << "'/>" << std::endl;
    fp << "      </Points>" << std::endl;
    Offset += 3 * nP * sizeof(float) + sizeof(int);
    fp << "      <Cells>" << std::endl;
    fp << "        <DataArray type='Int32' Name='connectivity' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += _part->iconnectivityLoc.size() * sizeof (int) + sizeof(int);
    fp << "        <DataArray type='Int32' Name='offsets' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += nC * sizeof (int) + sizeof(int);
    fp << "        <DataArray type='UInt8' Name='types' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += nC * sizeof (unsigned char) + sizeof(int);
    fp << "      </Cells>" << std::endl;
    fp << "      <PointData>" << std::endl;
    fp << "      </PointData>" << std::endl;
    fp << "      <CellData>" << std::endl;
    fp << "        <DataArray type='Int32' Name='vState' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += nC * sizeof (int) + sizeof (int);
    fp << "        <DataArray type='Int32' Name='gID' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += nC * sizeof (int) + sizeof (int);
    fp << "        <DataArray type='Float32' Name='IPFz' NumberOfComponents='3' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += 3*nC * sizeof (float) + sizeof (int);
    fp << "        <DataArray type='Float32' Name='Temperature' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += nC * sizeof (float) + sizeof (int);
    fp << "      </CellData>" << std::endl;
    fp << "    </Piece>" << std::endl;
    fp.close();
  } // for (int rank ...
  MPI_Barrier(MPI_COMM_WORLD);
  if (_part->myid == 0) {
    fp.open(vtkFilename.c_str(), std::fstream::app);
    fp << "  </UnstructuredGrid>" << std::endl;
    fp << "  <AppendedData encoding='raw'>" << std::endl;
    fp << "_";
    fp.close();
  }
  for (int rank = 0; rank < _part->nprocs; rank++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (_part->myid != rank) continue;
    fp.open(vtkFilename.c_str(), std::fstream::app);
    int Scalar = nP * sizeof (float), Vector = 3 * Scalar, Cells = nC * sizeof(int), 
      CellsTH=4*nC*sizeof(float), CellsTemp=nC*sizeof(float);
    int CellChars = nC * sizeof(unsigned char), Conn = _part->iconnectivityLoc.size() * sizeof(int);
    fp.write(reinterpret_cast<const char *>(&Vector), 4);
    for (int j=0;j<nP;j++) {
      j3 = floor(_part->ipointidLoc[j]/( (_xyz->nX[0]+1)*(_xyz->nX[1]+1)));
      j2 = floor( (_part->ipointidLoc[j] - (_xyz->nX[0]+1)*(_xyz->nX[1]+1)*j3)/(_xyz->nX[0]+1));
      j1 = _part->ipointidLoc[j] - (_xyz->nX[0]+1)*(_xyz->nX[1]+1)*j3 - (_xyz->nX[0]+1)*j2;
      float x = j1*_xyz->dX[0], y = j2*_xyz->dX[1], z = j3*_xyz->dX[2];
      fp.write(reinterpret_cast<const char *>(&x), sizeof(float));
      fp.write(reinterpret_cast<const char *>(&y), sizeof(float));
      fp.write(reinterpret_cast<const char *>(&z), sizeof(float));
    }
    fp.write(reinterpret_cast<const char *>(&Conn), 4);
    for (int i=0;i<nC;++i) {
      for (int j=0;j<_xyz->nnodePerCell;++j) {
	j1 = _part->iconnectivityLoc[_xyz->nnodePerCell*i+j];
	fp.write(reinterpret_cast<const char *>(&j1), sizeof(int));
      }
    }
    fp.write(reinterpret_cast<const char *>(&Cells), 4);
    cell_offsets = 0;
    for (int i=0;i<nC;i++) {
      cell_offsets += _xyz->nnodePerCell;
      fp.write(reinterpret_cast<const char *>(&cell_offsets), sizeof(int));
    }
    fp.write(reinterpret_cast<const char *>(&CellChars), 4);    
    for (int i=0;i<nC;i++) fp.write(reinterpret_cast<const char *>(&cell_type), sizeof (unsigned char));
    fp.write(reinterpret_cast<const char *>(&Cells), 4);
    for (int i=0;i<nC;i++) fp.write(reinterpret_cast<const char *>(&vState[i]), sizeof (int));
    fp.write(reinterpret_cast<const char *>(&Cells), 4);
    for (int i=0;i<nC;i++) fp.write(reinterpret_cast<const char *>(&gID[i]), sizeof (int));
    fp.write(reinterpret_cast<const char *>(&CellsTH), 4);
    for (int i=0;i<3*nC;i++) fp.write(reinterpret_cast<const char *>(&IPFmapBD[i]), sizeof (float));
    fp.write(reinterpret_cast<const char *>(&CellsTemp), 4);
    for (int i=0;i<nC;i++) fp.write(reinterpret_cast<const char *>(&TempOut[i]), sizeof (float));
    fp.close();
  } // for (int rank ...)
  MPI_Barrier(MPI_COMM_WORLD);
  if (_part->myid == 0) {
    fp.open(vtkFilename.c_str(), std::fstream::app);
    fp << "  </AppendedData>" << std::endl;
    fp << "</VTKFile>" << std::endl;
    fp.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);  
} // end WriteToVTU1

void VoxelsCA::WriteToVTU0(const std::string &filename)
{
  // writes gID, vState, cTheta per voxel
  std::string vtkFilename = filename + ".vtu";
  std::ofstream fp;
  int Offset=0, nC=_part->ncellLoc, nP=_part->npointLoc,j1,j2,j3,cell_offsets;
  unsigned char cell_type;
  std::vector< float> cThetaOut(4*nC,0),TempOut(nC,0);
  for (int j=0;j<nC;++j){
    TempOut[j] = _temp->TempCurr[j];
    if (gID[j]<1){
      cThetaOut[4*j] = 0;
      cThetaOut[4*j+1] = 0;
      cThetaOut[4*j+2] = 0;
      cThetaOut[4*j+3] = 0;
    } else {
      cThetaOut[4*j] = cTheta[4*(gID[j]-1)];
      cThetaOut[4*j+1] = cTheta[4*(gID[j]-1)+1];
      cThetaOut[4*j+2] = cTheta[4*(gID[j]-1)+2];
      cThetaOut[4*j+3] = cTheta[4*(gID[j]-1)+3];
    }
  }
  if (_xyz->nnodePerCell==4){cell_type=9;}
  if (_xyz->nnodePerCell==8){cell_type=12;}
  if (_part->myid == 0) {
    fp.open(vtkFilename.c_str());
    if (!fp) throw std::runtime_error("Cannot create VTU file: " + filename);
    fp << "<?xml version='1.0'?>" << std::endl;
    fp << "<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'>" << std::endl;
    fp << "  <UnstructuredGrid>" << std::endl;
    fp << "    <FieldData>" << std::endl;
    fp << "      <DataArray type='Int32' Name='time_step' NumberOfTuples='1' format='ascii'>" << std::endl;
    fp << "        " << _xyz->tInd << std::endl;
    fp << "      </DataArray>" << std::endl;
    fp << "    </FieldData>" << std::endl;
    fp.close();
  } // if (myid==0)
  for (int rank=0;rank<_part->nprocs;++rank){
    int OffsetP=Offset;
    MPI_Allreduce(&OffsetP,&Offset,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (_part->myid != rank) { Offset=0; continue; }
    fp.open(vtkFilename.c_str(), std::fstream::app);
    fp << "    <Piece NumberOfPoints='" << nP << "' NumberOfCells='" << nC << "'>" << std::endl;
    fp << "      <Points>" << std::endl;
    fp << "        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='appended' offset='" << Offset << "'/>" << std::endl;
    fp << "      </Points>" << std::endl;
    Offset += 3 * nP * sizeof(float) + sizeof(int);
    fp << "      <Cells>" << std::endl;
    fp << "        <DataArray type='Int32' Name='connectivity' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += _part->iconnectivityLoc.size() * sizeof (int) + sizeof(int);
    fp << "        <DataArray type='Int32' Name='offsets' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += nC * sizeof (int) + sizeof(int);
    fp << "        <DataArray type='UInt8' Name='types' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += nC * sizeof (unsigned char) + sizeof(int);
    fp << "      </Cells>" << std::endl;
    fp << "      <PointData>" << std::endl;
    fp << "      </PointData>" << std::endl;
    fp << "      <CellData>" << std::endl;
    fp << "        <DataArray type='Int32' Name='vState' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += nC * sizeof (int) + sizeof (int);
    fp << "        <DataArray type='Int32' Name='gID' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += nC * sizeof (int) + sizeof (int);
    fp << "        <DataArray type='Float32' Name='cTheta' NumberOfComponents='4' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += 4*nC * sizeof (float) + sizeof (int);
    fp << "        <DataArray type='Float32' Name='Temperature' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += nC * sizeof (float) + sizeof (int);
    fp << "      </CellData>" << std::endl;
    fp << "    </Piece>" << std::endl;
    fp.close();
  } // for (int rank ...
  MPI_Barrier(MPI_COMM_WORLD);
  if (_part->myid == 0) {
    fp.open(vtkFilename.c_str(), std::fstream::app);
    fp << "  </UnstructuredGrid>" << std::endl;
    fp << "  <AppendedData encoding='raw'>" << std::endl;
    fp << "_";
    fp.close();
  }
  for (int rank = 0; rank < _part->nprocs; rank++) {
    MPI_Barrier(MPI_COMM_WORLD);
    if (_part->myid != rank) continue;
    fp.open(vtkFilename.c_str(), std::fstream::app);
    int Scalar = nP * sizeof (float), Vector = 3 * Scalar, Cells = nC * sizeof(int), 
      CellsTH=4*nC*sizeof(float), CellsTemp=nC*sizeof(float);
    int CellChars = nC * sizeof(unsigned char), Conn = _part->iconnectivityLoc.size() * sizeof(int);
    fp.write(reinterpret_cast<const char *>(&Vector), 4);
    for (int j=0;j<nP;j++) {
      j3 = floor(_part->ipointidLoc[j]/( (_xyz->nX[0]+1)*(_xyz->nX[1]+1)));
      j2 = floor( (_part->ipointidLoc[j] - (_xyz->nX[0]+1)*(_xyz->nX[1]+1)*j3)/(_xyz->nX[0]+1));
      j1 = _part->ipointidLoc[j] - (_xyz->nX[0]+1)*(_xyz->nX[1]+1)*j3 - (_xyz->nX[0]+1)*j2;
      float x = j1*_xyz->dX[0], y = j2*_xyz->dX[1], z = j3*_xyz->dX[2];
      fp.write(reinterpret_cast<const char *>(&x), sizeof(float));
      fp.write(reinterpret_cast<const char *>(&y), sizeof(float));
      fp.write(reinterpret_cast<const char *>(&z), sizeof(float));
    }
    fp.write(reinterpret_cast<const char *>(&Conn), 4);
    for (int i=0;i<nC;++i) {
      for (int j=0;j<_xyz->nnodePerCell;++j) {
	j1 = _part->iconnectivityLoc[_xyz->nnodePerCell*i+j];
	fp.write(reinterpret_cast<const char *>(&j1), sizeof(int));
      }
    }
    fp.write(reinterpret_cast<const char *>(&Cells), 4);
    cell_offsets = 0;
    for (int i=0;i<nC;i++) {
      cell_offsets += _xyz->nnodePerCell;
      fp.write(reinterpret_cast<const char *>(&cell_offsets), sizeof(int));
    }
    fp.write(reinterpret_cast<const char *>(&CellChars), 4);    
    for (int i=0;i<nC;i++) fp.write(reinterpret_cast<const char *>(&cell_type), sizeof (unsigned char));
    fp.write(reinterpret_cast<const char *>(&Cells), 4);
    for (int i=0;i<nC;i++) fp.write(reinterpret_cast<const char *>(&vState[i]), sizeof (int));
    fp.write(reinterpret_cast<const char *>(&Cells), 4);
    for (int i=0;i<nC;i++) fp.write(reinterpret_cast<const char *>(&gID[i]), sizeof (int));
    fp.write(reinterpret_cast<const char *>(&CellsTH), 4);
    for (int i=0;i<4*nC;i++) fp.write(reinterpret_cast<const char *>(&cThetaOut[i]), sizeof (float));
    fp.write(reinterpret_cast<const char *>(&CellsTemp), 4);
    for (int i=0;i<nC;i++) fp.write(reinterpret_cast<const char *>(&TempOut[i]), sizeof (float));
    fp.close();
  } // for (int rank ...)
  MPI_Barrier(MPI_COMM_WORLD);
  if (_part->myid == 0) {
    fp.open(vtkFilename.c_str(), std::fstream::app);
    fp << "  </AppendedData>" << std::endl;
    fp << "</VTKFile>" << std::endl;
    fp.close();
  }
  MPI_Barrier(MPI_COMM_WORLD);  
} // end WriteToVTU0

void VoxelsCA::WriteCSVData(const std::string &filename)
{
  // write out csv file with grain id, axis-angle, vstate
  int j3,j2,j1,Ntot;
  double x,y,z,vvol;
  std::vector<double> gVol(nGrain,0.0),gVolT(nGrain,0);
  std::ofstream fp;
  Ntot = _part->ncellLoc;
  vvol = _xyz->dX[0]*_xyz->dX[1]*_xyz->dX[2];
  for (int j=0;j<Ntot;++j){
    if (gID[j]>0){gVol[gID[j]-1]+=vvol;}
  }
  MPI_Allreduce(&gVol[0],&gVolT[0],nGrain,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  if (_part->myid==0){
    fp.open(filename.c_str());
    fp << "grain nucleation (x,y,z), axis-angle (omega,n), grain Volume" << std::endl;
    for (int j=0;j<nGrain;++j){
      j3 = floor(gNucleus[j]/(_xyz->nX[0]*_xyz->nX[1]));
      j2 = floor( (gNucleus[j]-_xyz->nX[0]*_xyz->nX[1]*j3)  /_xyz->nX[0]);
      j1 = gNucleus[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      x = (double(j1)+.5)*(_xyz->dX[0]);
      y = (double(j2)+.5)*(_xyz->dX[1]);
      z = (double(j3)+.5)*(_xyz->dX[2]);
      fp << x <<"," << y << "," << z << "," << cTheta[4*j] << "," <<
	cTheta[4*j+1] << "," << cTheta[4*j+2] << "," << cTheta[4*j+3]<<","<<gVolT[j] << std::endl;
    } // end for (int j...
    fp.close();
  } // if (_part->myid ..)
  MPI_Barrier(MPI_COMM_WORLD);
} // WriteCSVData

void VoxelsCA::NucleateGrains(std::vector<int> &nIn, std::vector<double> &tIn)
{
  int Ntot=_part->ncellLoc,Nnew1,Nq=0,NqA,ind1[Ntot],Ntot2=_part->ncellLoc+_part->nGhost,cc,cc1;
  double tNp1= _xyz->time + _temp->DelT,rNuc[Ntot],x1,x2,rmax=0.0,rateX=0.0;
  std::fill(rNuc,rNuc+Ntot,0.0);
  std::vector<double> TempNp1(Ntot2,0.0);
  std::default_random_engine g1(30*_xyz->tInd+seed1);
  _temp->SchwalbachTempCurr(tNp1,TempNp1);
  for (int j=0;j<Ntot;++j){
    if (vState[j]==1 && _temp->TempCurr[j]< _xyz->tL){
      x2 = ( (_xyz->tL -TempNp1[j]) - _xyz->dTempM)/_xyz->dTempS/pow(2.0,.5);
      x1 = ( (_xyz->tL - _temp->TempCurr[j]) - _xyz->dTempM)/_xyz->dTempS/pow(2.0,.5);
      rNuc[j] = _xyz->rNmax*std::max(0.0, .5*(erf(x2) - erf(x1)));
      rmax = std::max(rmax,rNuc[j]);
      ind1[Nq] = j;
      Nq+=1;
    } // if (vState[j]==1...
  } // for (int j=0...
  MPI_Allreduce(&Nq,&NqA,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  rateX= (double)NqA *_xyz->dX[0]*_xyz->dX[1]*_xyz->dX[2] * rmax;
  // sample Poisson distribution from rmax and then do thinning
  std::poisson_distribution<int> Np(rateX);
  std::uniform_int_distribution<int> uind(0,Nq-1);
  std::uniform_real_distribution<double> xrand1(0.0,1.0);
  Nnew1 = (NqA==0)? 0:  (int)floor( ((double)Nq/(double)NqA)* Np(g1) );
  Nnew1 = std::min(Nnew1,Nq);
  int itmp1[Nnew1],itmp2[Nq];
  double tmp1[Nnew1];
  std::fill(itmp2,itmp2+Nq,0);
  cc1=0;
  while (cc1<Nnew1){
    cc = uind(g1);
    itmp2[cc]+=1;
    if (itmp2[cc]==1){
      itmp1[cc1] = ind1[cc];
      tmp1[cc1] = _temp->DelT*xrand1(g1); // what time to be 0 at current time
      cc1+=1;
    }
  } // while (cc1
  nIn.resize(Nnew1,0);
  tIn.resize(Nnew1,0);
  cc=0;
  //MPI_Barrier(MPI_COMM_WORLD);
  for (int j=0;j<Nnew1;++j){
    if (xrand1(g1) > (rNuc[itmp1[j]]/rmax)){continue;}
    nIn[cc] =itmp1[j];
    tIn[cc] = tmp1[j];
    cc+=1;
  } // for (int j ...
  nIn.resize(cc);
  tIn.resize(cc);
} // end NucleateGrains

void VoxelsCA::WriteToPVD(const std::string &filename, const std::vector<int> & filinds,
			  const std::vector<double> & filtime)
{
  std::string vtkFilename = filename + ".pvd";
  std::ofstream fp;
  std::string filvtu;
  if (_part->myid == 0) {
    fp.open(vtkFilename.c_str());
    fp << "<?xml version=\"1.0\"?>" << std::endl;
    fp << "<VTKFile type=\"Collection\" version=\"0.1\"" << std::endl;
    fp << "      byte_order=\"LittleEndian\"" << std::endl;
    fp << "      compressor=\"vtkZLibDataCompressor\">" << std::endl;
    fp << "  <Collection>" << std::endl;
    for (int j=0;j<filinds.size();++j){
      filvtu = filename+std::to_string(filinds[j])+".vtu";
      fp << "<DataSet timestep=\" "<< filtime[j] <<"\" group=\"\" part=\"0\"" << std::endl;
      fp << "   file=\"" << filvtu <<"\"/>" << std::endl;
    } // for (int j ...
    fp << "  </Collection>" << std::endl;
    fp <<"</VTKFile>" << std::endl;
  } // if (_part->myid
} // WriteToPVD
