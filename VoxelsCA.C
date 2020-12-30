// define member functions of VoxelCA

#include "Grid.h"
#include "BasePlate.h"
#include "VoxelsCA.h"
#include "Partition.h"
#include "iostream"
#include "fstream"
#include "sstream"
#include <math.h>
#include <algorithm>
#include <random>
#include <ctime>
#include "mpi.h"
#include "SampleOrientation.h"
#include <chrono>
#include <adios2.h>


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
  extentsInitialValue = 0.0; //(xyz->dx)/2.0
  extents.assign(Ntot1,extentsInitialValue);
  centroidOct.assign(3*Ntot1,0.0);
  seed0= 2132512;
  seed1=2912351;
  double velY=(5.51*pow(M_PI,2.0)*pow((- _xyz->mL)*(1-_xyz->kP),1.5)*
		 (_xyz->Gamma))*( pow((_xyz->tL - _xyz->tS),2.5)/pow(_xyz->c0,1.5));
  
  vXi = 3*_temp->bmV/velY;
  // establishes ineighID and ineighptr for convertSolid1 
  int cc=0;
  std::vector<int> neigh;
  ineighptr.assign(Ntot1+1,0);
  for (int j=0; j < Ntot1;++j){
    _xyz->ComputeNeighborhood(_part->icellidLoc[j],_xyz->neighOrder,neigh);
    cc+=neigh.size();
    ineighptr[j+1]=cc;
  }
  ineighID.assign(cc,0);
  cc=0;
  for (int j=0;j<Ntot1;++j){
    _xyz->ComputeNeighborhood(_part->icellidLoc[j],_xyz->neighOrder,neigh);
    for (int j1=0;j1<neigh.size();++j1){
      ineighID[cc] = std::distance(_part->icellidLoc.begin(),
			 std::find(_part->icellidLoc.begin(),_part->icellidLoc.end(),neigh[j1]));
      cc+=1;
    }
  }

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

void VoxelsCA::InitializeTest2()
{
  int Ntot = _part->ncellLoc + _part->nGhost,j2,jn[3],i2,i3;
  i3 = _xyz->nX[0]*_xyz->nX[1];
  i2 = _xyz->nX[0];
  nGrain = 8;
  std::uniform_real_distribution<> ir(0.0,1.0);
  std::default_random_engine gen1(1234125);
  double ax[3],anorm,xi;

  gID.assign(Ntot,0);
  gNucleus.resize(nGrain,0);
  //gNucleus[0] = _xyz->nX[0]*_xyz->nX[1]*floor(_xyz->nX[2]/2.0)+_xyz->nX[0]*floor(_xyz->nX[1]/2.0)+floor(_xyz->nX[0]/2.0);
  xi = .8;
  for (int j=0;j<nGrain;++j){
    jn[2] = int(floor( _xyz->nX[2]*(xi*ir(gen1)) + xi/2.0) );
    jn[1] = int(floor( _xyz->nX[1]*(xi*ir(gen1)) + xi/2.0) );
    jn[0] = int(floor( _xyz->nX[0]*(xi*ir(gen1)) + xi/2.0) );
    gNucleus[j] = i3*jn[2]+i2*jn[1]+jn[0];
  }
  //cTheta = {M_PI*60.0/180.0,1.0/pow(3.0,.5),1.0/pow(3.0,.5),1.0/pow(3.0,.5)};
  unsigned int sdloc;
  sdloc = unsigned(double(gen1())/double(gen1.max())*pow(2.0,32.0));
  SampleOrientation sa1;
  // randomly generate crystallographic orientations
  std::vector<double> aa;
  sa1.GenerateSamples(nGrain,sdloc,aa);
  cTheta.assign(nGrain*4,0);
  for (int j1=0;j1<nGrain;++j1){
    cTheta[4*j1] = aa[4*j1];
    cTheta[4*j1+1] =aa[4*j1+1];
    cTheta[4*j1+2] = aa[4*j1+2];
    cTheta[4*j1+3] = aa[4*j1+3];
  } // end for j1        
  //std::cout << cTheta[0]<<","<<cTheta[1]<<","<<cTheta[2]<<","<<cTheta[3]<<","<<std::endl;
  for (int j =0; j <Ntot; ++j){ vState[j] = 1;}
  for (int j=0;j<nGrain;++j){
    j2 = std::distance(_part->icellidLoc.begin(),
		     std::find(_part->icellidLoc.begin(),_part->icellidLoc.end(),gNucleus[j]));
    if (j2 != _part->icellidLoc.size()){
      vState[j2] = 2;
      gID[j2] = j+1;
      jn[2] = floor(_part->icellidLoc[j2]/i3);
      jn[1] =floor( (_part->icellidLoc[j2]- i3*jn[2])/i2);
      jn[0] = _part->icellidLoc[j2] -i3*jn[2] - i2*jn[1];
      centroidOct[3*j2]=(double(jn[0])+.5)*_xyz->dX[0];
      centroidOct[3*j2+1]=(double(jn[1])+.5)*_xyz->dX[1];
      centroidOct[3*j2+2]=(double(jn[2])+.5)*_xyz->dX[2];
    }
  } // for (int j...
} // end InitializeTest2()

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
  /*
    This is identical to UpdateVoxels5() except it has nucleation. You
    must make sure that modifications to the code are implemented in
    both files
   */
  // set liquid if temperature > TL 
  // and zero voxels above laser (no material there)
  //*******************  AFTER TEST, ELIMINATE THIS IF/ELSE AND TEST CASES
  // TEST CASE IS ICTRL=4
  if (_xyz->ictrl==4){
    SetLiquid4(); 
  } else{
    SetLiquid3();  
    //ZeroVoxels1(); // voxels are zero'ed in update layer
  }
  //*******************  AFTER TEST, ELIMINATE THIS IF/ELSE AND TEST CASES
  // solid (vState=3) to mushy (vState=2) if one neighbor liquid (vState=1)
  ConvertSolid1(1);
  //_part->PassInformation(vState);
  // assemble arrays to be used to compute grain growth
  int NlocA=0,Na=0,Ntot,cc,cc1,jn[3],nHaz=0,nHazT=0,
    NnlocA=0,i1,i2,i3,i4,i5,i6,Nna=0;
  Ntot = _part->ncellLoc;
  std::vector<int> neigh,j0(_part->nprocs,0),jn0(_part->nprocs,0),
    NvecA(_part->nprocs),NnvecA(_part->nprocs),indA(26*Ntot,Ntot+10);
  std::vector<std::vector<double>> sdiag0,sdiag(6,std::vector<double>(3));
  std::vector<std::vector<int>> sInd;
  loadS(sdiag0,sInd); // this is for the decentered octohedron algorithm
  for (int j=0;j<Ntot;++j){
    if (vState[j]==2 || (vState[j]==1 && _temp->TempCurr[j]<_xyz->tL) ){
      indA[NlocA]=j;
      NlocA+=1;
      for (int j1=ineighptr[j];j1<ineighptr[j+1];++j1){	
	if (_temp->TempCurr[ineighID[j1]]>=_xyz->tL && ineighID[j1]<Ntot){
	  indA[NlocA]=ineighID[j1];
	  NlocA+=1;
	}
      }
    }
  }
  std::partial_sort(indA.begin(),indA.begin()+NlocA,indA.end());
  NlocA=std::distance(indA.begin(),std::unique(indA.begin(),indA.begin()+NlocA));
  MPI_Allgather(&NlocA,1,MPI_INT,&NvecA[0],1,MPI_INT,MPI_COMM_WORLD);
  for (int j=0;j<_part->nprocs;++j){Na+=NvecA[j];}
  std::vector<double> T(Na,0),ExtA(Na,0.0),CentroidA(3*Na,0);
  std::vector<int> vS(Na,0),G(Na,0),vI(Na,0),iv(_part->nprocs+1,0);
  j0[0]=0;
  for (int j=1;j<_part->nprocs;++j){j0[j]= j0[j-1]+NvecA[j-1];}
  cc=0; 
  for (int j=0;j<NlocA;++j){
    i1 = indA[j];
    vS[j0[_part->myid]+cc] = vState[i1];
    vI[j0[_part->myid]+cc] = _part->icellidLoc[i1];
    G[j0[_part->myid]+cc] = gID[i1];
    T[j0[_part->myid]+cc] = _temp->TempCurr[i1];
    ExtA[j0[_part->myid]+cc] = extents[i1];
    CentroidA[3*j0[_part->myid]+3*cc]=centroidOct[3*i1];
    CentroidA[3*j0[_part->myid]+3*cc+1]=centroidOct[3*i1+1];
    CentroidA[3*j0[_part->myid]+3*cc+2]=centroidOct[3*i1+2];
    cc+=1;
  } // for (int j...
  for (int j=0;j<_part->nprocs;++j){
    if (NvecA[j]==0){continue;}
    MPI_Bcast(&vS[j0[j]],NvecA[j],MPI_INT,j,MPI_COMM_WORLD);
    MPI_Bcast(&vI[j0[j]],NvecA[j],MPI_INT,j,MPI_COMM_WORLD);
    MPI_Bcast(&G[j0[j]],NvecA[j],MPI_INT,j,MPI_COMM_WORLD);
    MPI_Bcast(&T[j0[j]],NvecA[j],MPI_DOUBLE,j,MPI_COMM_WORLD);
    MPI_Bcast(&ExtA[j0[j]],NvecA[j],MPI_DOUBLE,j,MPI_COMM_WORLD);
    MPI_Bcast(&CentroidA[3*j0[j]],3*NvecA[j],MPI_DOUBLE,j,MPI_COMM_WORLD);
  } // for (int j ..
  std::vector<int> inloc(ineighptr[Ntot],0),ineighAptr(Na+1,0);
  NnlocA=0;
  cc=0;
  for (int j=0;j<NlocA;++j){
    i3 = indA[j];
    ineighAptr[j0[_part->myid]+cc]=NnlocA;
    for (int j1=ineighptr[i3];j1<ineighptr[i3+1];++j1){
      i1 = _part->icellidLoc[ineighID[j1]];
      i2 = std::distance(vI.begin(),std::find(vI.begin(),vI.end(),i1));
      if (i2<Na){
	inloc[NnlocA] = i2;
	NnlocA+=1;
      } // if (i2<Na...
    } // for (int j1 ...
    cc+=1;
  } // for (int j...
  MPI_Allgather(&NnlocA,1,MPI_INT,&NnvecA[0],1,MPI_INT,MPI_COMM_WORLD);
  jn0[0]=0;
  for (int j=1;j<_part->nprocs;++j){jn0[j]= jn0[j-1]+NnvecA[j-1];}
  for (int j=0;j<_part->nprocs;++j){Nna+=NnvecA[j];}
  for (int j=0;j<NlocA;++j){ineighAptr[j0[_part->myid]+j]+=jn0[_part->myid];}
  for (int j=0;j<_part->nprocs;++j){
    if (NvecA[j]==0){continue;}
    MPI_Bcast(&ineighAptr[j0[j]],NvecA[j],MPI_INT,j,MPI_COMM_WORLD);
  } // for (int j ..
  ineighAptr[Na]=Nna;
  std::vector<int> ineighIDA(Nna,0);
  for (int j=0;j<NnlocA;++j){ineighIDA[jn0[_part->myid]+j]=inloc[j];}
  for (int j=0;j<_part->nprocs;++j){
    if (NnvecA[j]==0){continue;}
    MPI_Bcast(&ineighIDA[jn0[j]],NnvecA[j],MPI_INT,j,MPI_COMM_WORLD);
  }
  i1 = ceil( (double)Na / (double)_part->nprocs);
  i2 = floor( (double)Na / (double)i1);
  for (int j=0;j<(_part->nprocs+1);++j){
    if (j < (i2+1)){iv[j] = i1*j;}
    if (j>i2 && j<_part->nprocs){iv[j] = i2*i1 + 
	floor( (double)(Na-i2*i1)/(double)(_part->nprocs-i2));}
    if (j==_part->nprocs){iv[j]=Na;}
  } // for (int j...

  // nucleate grains
  /*
  int NnucA;
  std::vector<int> nucA;
  std::vector<double> tnucA,quatnuc;
  NucleateGrains(nucA,tnucA,quatnuc,iv,vI,vS,T);
  NnucA=nucA.size();

  }
  */
  unsigned int sdloc;
  SampleOrientation sa;
  sdloc= seed0 + 32*_xyz->tInd +64*nGrain;
  std::vector<double> quatnuc(4);
  std::default_random_engine g1(30*_xyz->tInd+seed1);
  std::uniform_real_distribution<double> xrand1(0.0,1.0);
  double rX;
  // set nucleation discrete
  rX = _xyz->rNmax*pow(_xyz->dX[0]*1e6,3.);
  // end nucleation grains

  // end assemble arrays to be used to compute grain growth
  // capture all undercooled liquid voxels by growing grains
  int js=0,j1s=0,jx[3],jy[3],countS=Na;
  double velX,velY,omega,ax[3],
    dlocX[3],locX[3],dr,vhat,tmelt=_xyz->tL,ph,th,
    timeUntil,dnx[3],tinc=0.0;
  std::vector<int> vSneigh(26,0),jTs(_part->nprocs,0),j1Ts(_part->nprocs,0);
  std::vector<double> Tneigh(26,0.0),DtT(_part->nprocs,0.0);
  std::vector<std::vector<double>> rRot(3,std::vector<double>(3,0));
  i2 = _xyz->nX[0]; i3 = _xyz->nX[0]*_xyz->nX[1];
  cc=0;
  struct{double DtMin; int rank;} xin,xout;
  xin.rank=_part->myid;
  
  while (std::count(vS.begin(),vS.end(),2)!=countS && tinc<_temp->DelT) {
  //while (std::count(vS.begin(),vS.end(),2)!=countS) {
    xin.DtMin=1e6;
    countS = std::count(vS.begin(),vS.end(),2);
    i1 = iv[_part->myid];
    std::vector<double> vhatvec(iv[_part->myid+1]-i1,0);
    for (int j=iv[_part->myid];j<iv[_part->myid+1];++j){
      if (vS[j]==2){
        i5 = ineighAptr[j];
        i6 = ineighAptr[j+1]-i5;
        for (int j1=0;j1<i6;++j1){
          i4 = ineighIDA[i5+j1];
          vSneigh[j1]=vS[i4];
          Tneigh[j1] = T[i4];
        } // for (int j1... 
	if (std::any_of(vSneigh.begin(),vSneigh.begin()+i6,[](int n){return n==1;}) &&
	    std::all_of(Tneigh.begin(),Tneigh.begin()+i6,[&tmelt](double tchk){return tchk < tmelt;})){
	  omega = cTheta[4*(G[j]-1)];
	  ax[0]=cTheta[4*(G[j]-1)+1];
	  ax[1]=cTheta[4*(G[j]-1)+2];
	  ax[2]=cTheta[4*(G[j]-1)+3];
	  loadRotMat(omega,ax,rRot);	    
	  T[j]>=_xyz->tL ? velY=0.0: velY=(5.51*pow(M_PI,2.0)*pow((- _xyz->mL)*(1-_xyz->kP),1.5)*
		 (_xyz->Gamma))*( pow((_xyz->tL - T[j]),2.5)/pow(_xyz->c0,1.5));
	  vhatvec[j-i1] = vXi*velY;
	  for (int j1=0;j1<i6;++j1){
	    if (vSneigh[j1] != 1 ){continue;}
            i4 = ineighIDA[i5+j1];
	    jx[2] = floor(vI[i4]/i3);
	    jx[1] = floor((vI[i4]- i3*jx[2])/i2);
	    jx[0] = vI[i4] - i3*jx[2] - i2*jx[1];
	    dnx[0] = (double(jx[0])+.5)*_xyz->dX[0] - CentroidA[3*j];
	    dnx[1] = (double(jx[1])+.5)*_xyz->dX[1] - CentroidA[3*j+1];
	    dnx[2] = (double(jx[2])+.5)*_xyz->dX[2] - CentroidA[3*j+2];
            th = atan2(std::fabs(dnx[1]),std::fabs(dnx[0]));
            th > M_PI/4.0 ? th= M_PI/2.0 - th: th ;
            ph = atan2(pow(pow(dnx[0],2.0)+pow(dnx[1],2.0),.5),std::fabs(dnx[2]));
            ph < M_PI/4.0 ? ph = M_PI/2.0 - ph: ph ;
	    // matrix is local->global; need to multiply by transpose for global->local            
	    // put into 1st quadrant b/c of symmetry     
	    dlocX[0] = std::fabs(rRot[0][0]*dnx[0]+rRot[1][0]*dnx[1]+rRot[2][0]*dnx[2]);
	    dlocX[1] = std::fabs(rRot[0][1]*dnx[0]+rRot[1][1]*dnx[1]+rRot[2][1]*dnx[2]);
	    dlocX[2] = std::fabs(rRot[0][2]*dnx[0]+rRot[1][2]*dnx[1]+rRot[2][2]*dnx[2]);
	    dr = pow(cos(th)*sin(ph),.5)*(dlocX[0]+dlocX[1]+dlocX[2]) - ExtA[j];
	    timeUntil = dr/vhatvec[j-i1];
	    if (timeUntil < xin.DtMin){
	      xin.DtMin = timeUntil;
	      js = j;  // js is capturing voxel index
	      j1s = j1;  // V[js][j1s] is index of captured voxel
	    } // if (timeUntil ...
	  } // for (int j1...
	} // if (std::any_of ...  
      } // if (vS[j]==2...      
    } // for (int j...
    MPI_Allreduce(&xin,&xout,1,MPI_DOUBLE_INT,MPI_MINLOC,MPI_COMM_WORLD);
    MPI_Bcast(&js,1,MPI_INT,xout.rank,MPI_COMM_WORLD);
    MPI_Bcast(&j1s,1,MPI_INT,xout.rank,MPI_COMM_WORLD);
    if (xout.DtMin>=1e6){break;}
    // test: make captured grain new grain based on rate
    //rX = _xyz->rNmax*exp( - 25*pow( (T[js]-_xyz->tS)/(_xyz->tL-_xyz->tS) ,2.0));
    if (xrand1(g1)< rX){      
      i1 = ineighIDA[ineighAptr[js]+j1s];
      sa.GenerateSamples(1,sdloc,quatnuc);
      vS[i1] = 2;
      nGrain+=1;
      G[i1] = nGrain;
      cTheta.insert(cTheta.end(),quatnuc.begin(),quatnuc.end());
      jx[2] = floor(vI[i1]/i3);
      jx[1] =floor( (vI[i1]- i3*jx[2])/i2);
      jx[0] = vI[i1] -i3*jx[2] - i2*jx[1];
      CentroidA[3*i1] = (double(jx[0])+.5)*_xyz->dX[0];
      CentroidA[3*i1+1] = (double(jx[1])+.5)*_xyz->dX[1];
      CentroidA[3*i1+2] = (double(jx[2])+.5)*_xyz->dX[2];	
      sdloc +=10;
      xout.DtMin=0.0;
      ExtA[i1] = extentsInitialValue;

    /*
    i4 = std::distance(tnucA.begin(),std::upper_bound(tnucA.begin(),tnucA.end(),tinc));
    i5 = std::distance(tnucA.begin(),std::upper_bound(tnucA.begin(),tnucA.end(),tinc+xout.DtMin));
    cc1=0;    
    while (i4<i5){
      if (vS[nucA[i4]]==1){
	cc1=1;
	break;
      }
      i4+=1;
    }
    if (cc1==1){

      
      xout.DtMin = tnucA[i4]-tinc;
      js = nucA[i4];
      i1=iv[_part->myid];
      for (int j=iv[_part->myid];j<iv[_part->myid+1];++j){
	if (vS[j]!=2){continue;}
	ExtA[j]+=vhatvec[j-i1]*xout.DtMin;
	ExtA[j] = std::max(ExtA[j],0.0);
      } // for (int j...
      vS[js] = 2;
      nGrain+=1;
      G[js] = nGrain;
      cTheta.insert(cTheta.end(),quatnuc.begin()+4*i4,quatnuc.begin()+4*i4+4);
      jx[2] = floor(vI[js]/i3);
      jx[1] = floor((vI[js]- i3*jx[2])/i2);
      jx[0] = vI[js] - i3*jx[2] - i2*jx[1];
      CentroidA[3*js] = (double(jx[0])+.5)*_xyz->dX[0];
      CentroidA[3*js+1] = (double(jx[1])+.5)*_xyz->dX[1];
      CentroidA[3*js+2] = (double(jx[2])+.5)*_xyz->dX[2];
      ExtA[js] = 0.0;    
      */
    } else {
      double xI[3],xJ[3],d1I,dI2,d1J,dJ3,L12,L13,l,Lmud,dnorm[3],xiL;
      int jInd,nvoxproc;
      // compute new decentered octohedron return centroid and extents
      // grow all extents that satify appropriate condition
      i1=iv[_part->myid];
      for (int j=iv[_part->myid];j<iv[_part->myid+1];++j){
	if (vS[j]!=2){continue;}
	ExtA[j]+=vhatvec[j-i1]*xout.DtMin;
	ExtA[j] = std::max(ExtA[j],0.0);
      } // for (int j... 
      MPI_Bcast(&ExtA[js],1,MPI_DOUBLE,xout.rank,MPI_COMM_WORLD);
      i1 = ineighIDA[ineighAptr[js]+j1s];
      vS[i1] = 2;
      G[i1] = G[js];
      jx[2] = floor(vI[i1]/i3);
      jx[1] =floor( (vI[i1]- i3*jx[2])/i2);
      jx[0] = vI[i1] -i3*jx[2] - i2*jx[1];
      jy[2] = floor(vI[js]/i3);
      jy[1] = floor((vI[js]- i3*jy[2])/i2);
      jy[0] = vI[js] - i3*jy[2] - i2*jy[1];
      l = pow( pow((jx[0]-jy[0])*_xyz->dX[0],2)+ pow((jx[1]-jy[1])*_xyz->dX[1],2)+
	       pow((jx[2]-jy[2])*_xyz->dX[2],2),.5);
      dnx[0] = (double(jx[0])+.5)*_xyz->dX[0] - CentroidA[3*js];
      dnx[1] = (double(jx[1])+.5)*_xyz->dX[1] - CentroidA[3*js+1];
      dnx[2] = (double(jx[2])+.5)*_xyz->dX[2] - CentroidA[3*js+2];
      omega = cTheta[4*(G[js]-1)];
      ax[0]=cTheta[4*(G[js]-1)+1];
      ax[1]=cTheta[4*(G[js]-1)+2];
      ax[2]=cTheta[4*(G[js]-1)+3];
      loadRotMat(omega,ax,rRot);	    
      // matrix is local->global; need to multiply by transpose for global->local            
      locX[0] = rRot[0][0]*dnx[0]+rRot[1][0]*dnx[1]+rRot[2][0]*dnx[2];
      locX[1] = rRot[0][1]*dnx[0]+rRot[1][1]*dnx[1]+rRot[2][1]*dnx[2];
      locX[2] = rRot[0][2]*dnx[0]+rRot[1][2]*dnx[1]+rRot[2][2]*dnx[2];
      th = atan2(std::fabs(dnx[1]),std::fabs(dnx[0]));
      th > M_PI/4.0 ? th= M_PI/2.0 - th: th;
      ph = atan2(pow(pow(dnx[0],2.0)+pow(dnx[1],2.0),.5),std::fabs(dnx[2]));
      ph < M_PI/4.0 ? ph = M_PI/2.0 - ph: ph ;
      // signbit returns 0 if positive and 1 if negative
      dr = std::fabs(locX[0])+ std::fabs(locX[1])+ std::fabs(locX[2]);
      for (int j1=0;j1<6;++j1){sdiag[j1]={sdiag0[j1][0]*dr,sdiag0[j1][1]*dr,sdiag0[j1][2]*dr};}
      jInd = 4*std::signbit(locX[2])+ 2*std::signbit(locX[1])+ std::signbit(locX[0]);
      for (int j1=0;j1<3;++j1){ 
	dnorm[j1] = pow(locX[0]-sdiag[sInd[jInd][j1]][0],2.0)+
	  pow(locX[1]-sdiag[sInd[jInd][j1]][1],2.0)+pow(locX[2]-sdiag[sInd[jInd][j1]][2],2.0);
      } // for (int j1...
      jy[0] = std::distance(dnorm,std::min_element(dnorm,dnorm+3));
      jy[2] = std::distance(dnorm,std::max_element(dnorm,dnorm+3));
      if (jy[0]==jy[2]){
	  jy[0]=0;
	  jy[1]=1;
	  jy[2]=2;
      } else{
	jy[1] = 3 - jy[0] - jy[2];
      }
      projectPointLine(locX,&sdiag[sInd[jInd][jy[0]]][0],&sdiag[sInd[jInd][jy[1]]][0],xI);
      projectPointLine(locX,&sdiag[sInd[jInd][jy[0]]][0],&sdiag[sInd[jInd][jy[2]]][0],xJ);
      d1I = pow(pow(sdiag[sInd[jInd][jy[0]]][0]-xI[0],2.0) + 
		pow(sdiag[sInd[jInd][jy[0]]][1]-xI[1],2.0) + 
		pow(sdiag[sInd[jInd][jy[0]]][2]-xI[2],2.0),.5);
      dI2 = pow(pow(sdiag[sInd[jInd][jy[1]]][0]-xI[0],2.0) + 
		pow(sdiag[sInd[jInd][jy[1]]][1]-xI[1],2.0) + 
		pow(sdiag[sInd[jInd][jy[1]]][2]-xI[2],2.0),.5);
      d1J = pow(pow(sdiag[sInd[jInd][jy[0]]][0]-xJ[0],2.0) + 
		pow(sdiag[sInd[jInd][jy[0]]][1]-xJ[1],2.0) + 
		pow(sdiag[sInd[jInd][jy[0]]][2]-xJ[2],2.0),.5);
      dJ3 = pow(pow(sdiag[sInd[jInd][jy[2]]][0]-xJ[0],2.0) + 
		pow(sdiag[sInd[jInd][jy[2]]][1]-xJ[1],2.0) + 
		pow(sdiag[sInd[jInd][jy[2]]][2]-xJ[2],2.0),.5);
      L12 = .5*(std::min(d1I,pow(3.0,.5)*l) + std::min(dI2,pow(3.0,.5)*l) );  
      L13 = .5*(std::min(d1J,pow(3.0,.5)*l) + std::min(dJ3,pow(3.0,.5)*l) );  
      Lmud =  pow(3.0,.5)* (pow(2.0/3.0,.5)*std::max(L12,L13));
      //xiL = -.4*(l/_xyz->dX[0] - 1.0)/(pow(3.0,.5)-1.0) + .9;
      //xiL = -.2*(l/_xyz->dX[0] - 1.0)/(pow(3.0,.5)-1.0) + .6;
      //xiL = .4 + .13*(l/_xyz->dX[0]-1.0);
      //xiL = .25 + .75*(l/_xyz->dX[0]-1.0);
      xiL = 1.0;
      ExtA[i1] = pow(cos(th)*sin(ph),.5)*Lmud*xiL;
      /*
	dnx[0] = sdiag[sInd[jInd][jy[0]]][0] - (2-xiL)*Lmud*sdiag0[sInd[jInd][jy[0]]][0];
	dnx[1] = sdiag[sInd[jInd][jy[0]]][1] - (2-xiL)*Lmud*sdiag0[sInd[jInd][jy[0]]][1];
	dnx[2] = sdiag[sInd[jInd][jy[0]]][2] - (2-xiL)*Lmud*sdiag0[sInd[jInd][jy[0]]][2];
      */
      dnx[0] = sdiag[sInd[jInd][jy[0]]][0] - Lmud*sdiag0[sInd[jInd][jy[0]]][0];
      dnx[1] = sdiag[sInd[jInd][jy[0]]][1] - Lmud*sdiag0[sInd[jInd][jy[0]]][1];
      dnx[2] = sdiag[sInd[jInd][jy[0]]][2] - Lmud*sdiag0[sInd[jInd][jy[0]]][2];
      locX[0] = rRot[0][0]*dnx[0]+rRot[0][1]*dnx[1]+rRot[0][2]*dnx[2];
      locX[1] = rRot[1][0]*dnx[0]+rRot[1][1]*dnx[1]+rRot[1][2]*dnx[2];
      locX[2] = rRot[2][0]*dnx[0]+rRot[2][1]*dnx[1]+rRot[2][2]*dnx[2];
      CentroidA[3*i1] = CentroidA[3*js] + locX[0];
      CentroidA[3*i1+1] = CentroidA[3*js+1] + locX[1];
      CentroidA[3*i1+2] = CentroidA[3*js+2] + locX[2];      
      // end compute new decentered octohedron return centroid and extents      
    } // if (tnucA[i4]<(tinc+xout.DtMin)) ...
    cc+=1;    
    tinc += xout.DtMin;
  } // while (std::count...
  for (int j=0;j<_part->nprocs;++j){
    i1 = iv[j+1]-iv[j];
    MPI_Bcast(&ExtA[iv[j]],i1,MPI_DOUBLE,j,MPI_COMM_WORLD);  
  } // for (int j ...			
  // end capture all undercooled liquid voxels by growing grains
  // bring global variables back to local variables
  cc1=0;
  for (int j=0;j<NlocA;++j){
    i1 = indA[j];
    vState[i1] = vS[j0[_part->myid]+cc1];
    gID[i1] = G[j0[_part->myid]+cc1];
    extents[i1] = ExtA[j0[_part->myid]+cc1];
    centroidOct[3*i1] = CentroidA[3*j0[_part->myid]+3*cc1];
    centroidOct[3*i1+1] = CentroidA[3*j0[_part->myid]+3*cc1+1];
    centroidOct[3*i1+2] = CentroidA[3*j0[_part->myid]+3*cc1+2];
    cc1+=1;
  } // for (int j...
  // pass information 
  for (int j=Ntot;j<(_part->nGhost+Ntot);++j){
    cc1=std::distance(vI.begin(),std::find(vI.begin(),vI.end(),_part->icellidLoc[j]));
    if (cc1<Na){
      vState[j] = vS[cc1];
      gID[j] = G[cc1];
    } // if (cc1<
  } // for (int j...
  // mushy (vState=2) to solid (vState=3) if all neighbors 2 
  ConvertSolid1(0);
  _part->PassInformation(vState);
}; // end UpdateVoxels
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
  int Ntot = (_part->ncellLoc),cc,j2,cc1;
  std::vector<int> vneigh,neigh,i1(_part->ncellLoc,0);
  if (iswitch==0){
    // this checks if all neighbors are >=2 and converts to 3 if currently 2
    /*
    cc=0;
    for (int j=0; j < Ntot;++j){
      if (vState[j]==2){
	_xyz->ComputeNeighborhood(_part->icellidLoc[j],_xyz->neighOrder,neigh);
	vneigh.assign(neigh.size(),0);
	for (int j1=0;j1<neigh.size();++j1){
	  j2 = std::distance(_part->icellidLoc.begin(),
			     std::find(_part->icellidLoc.begin(),_part->icellidLoc.end(),neigh[j1]));
	  vneigh[j1] = vState[j2];
	}
	if (std::all_of(vneigh.begin(),vneigh.end(),[](int n){return n>=2;})) {
	//if (!std::any_of(vneigh.begin(),vneigh.end(),[](int n){return n==1;})) {
	  i1[cc] = j;
	  cc +=1;
	} // if
      } // if (vState[j]==2)
    } // for j
    */
    cc=0;
    for (int j=0; j < Ntot;++j){
      if (vState[j]==2){      
	vneigh.assign(ineighptr[j+1]-ineighptr[j],0);
	cc1=0;
	for (int j1=ineighptr[j];j1<ineighptr[j+1];++j1){
	  vneigh[cc1]=vState[ineighID[j1]];
	  cc1+=1;
	}
	if (std::all_of(vneigh.begin(),vneigh.end(),[](int n){return n>=2;})) {
	//if (!std::any_of(vneigh.begin(),vneigh.end(),[](int n){return n==1;})) {
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
    /*
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
//	if (!std::all_of(vneigh.begin(),vneigh.end(),[](int n){return n>1;})) {
	if (std::any_of(vneigh.begin(),vneigh.end(),[](int n){return n==1;})) {
	  i1[cc] = j;
	  cc += 1;
	} // end if
      } // if (vState[j]==3)
    } // for j
    */
    cc=0;
    for (int j=0; j < Ntot;++j){
      if (vState[j]==3){
	vneigh.assign(ineighptr[j+1]-ineighptr[j],0);
	cc1=0;
	for (int j1=ineighptr[j];j1<ineighptr[j+1];++j1){
	  vneigh[cc1]=vState[ineighID[j1]];
	  cc1+=1;
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
void VoxelsCA::CleanLayer(){
 /*
   at beginning of each layer, function purges indices of cTheta associated
   with grains having no volume and renumbers gID to eliminate grains
   with no volume
 */

  // purge indices of cTheta associated with grains having no volume
  std::vector<double> gVol(nGrain,0.0),gVolT(nGrain,0),ctmp;
  std::vector<int> gtmp,itmp1(nGrain,0);
  int Ntot,Ntot2,nGtmp,i1t;
  double vvol;
  Ntot = _part->ncellLoc;
  Ntot2 = _part->ncellLoc + _part->nGhost;
  vvol = _xyz->dX[0]*_xyz->dX[1]*_xyz->dX[2];
  for (int j=0;j<Ntot;++j){
    if (gID[j]>0){gVol[gID[j]-1]+=vvol;}
  }
  MPI_Allreduce(&gVol[0],&gVolT[0],nGrain,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  nGtmp=0;
  ctmp.assign(4*nGrain,0);
  for (int j=0;j<nGrain;++j){
    if (gVolT[j]>0.0){
      ctmp[4*nGtmp] = cTheta[4*j];
      ctmp[4*nGtmp+1] = cTheta[4*j+1];
      ctmp[4*nGtmp+2] = cTheta[4*j+2];
      ctmp[4*nGtmp+3] = cTheta[4*j+3];
      itmp1[j] = nGtmp+1;
      nGtmp+=1;
    }
  }
  cTheta.assign(ctmp.begin(),ctmp.begin()+4*nGtmp);
  nGrain = nGtmp;
  for (int j=0;j<Ntot2;++j){
    if (gID[j]==0){continue;}
    i1t = itmp1[gID[j]-1];
    gID[j] =  i1t;
  }
 
}; // CleanLayer
void VoxelsCA::AddLayer(){
  /*
    function adds the layer of powder and is called at the beginning of new layer
   */
  int Ntot,Ntot2,i1t,nVlayer,j1,j2,j3,iz1,jvox0;
  Ntot = _part->ncellLoc;
  Ntot2 = _part->ncellLoc + _part->nGhost;
  nVlayer = _xyz->nX[0]*_xyz->nX[1]*_xyz->nZlayer;
  std::default_random_engine gen1(3134525);
  unsigned int sdloc;
  sdloc = unsigned(double(gen1())/double(gen1.max())*pow(2.0,32.0));
  SampleOrientation sa1;
  // randomly generate crystallographic orientations
  std::vector<double> aa;
  sa1.GenerateSamples(nVlayer,sdloc,aa);
  cTheta.insert(cTheta.end(), aa.begin(),aa.end());
  // ASSIGN GID NGRAIN1+1...NGRAIN TO PROPER VOXELS
  // LOOK INTO WHETHER GNUCLEUS IS IMPORTANT IN WHICH
  // CASE IT WILL BE THE VOXEL CENTERS
  iz1 = _temp->ilaserLoc - _xyz->nZlayer;
  jvox0 = _xyz->nX[0]*_xyz->nX[1]*iz1;
  for (int j=0;j<Ntot2;++j){
    j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
    if (j3>=iz1 && j3<_temp->ilaserLoc){
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      gID[j] = nGrain + _xyz->nX[0]*_xyz->nX[1]*j3+_xyz->nX[0]*j2+j1 - jvox0 + 1;
      vState[j] = 3;
      if (j<Ntot){
	centroidOct[3*j]=(double(j1)+.5)*_xyz->dX[0];
	centroidOct[3*j+1]=(double(j2)+.5)*_xyz->dX[1];
	centroidOct[3*j+2]=(double(j3)+.5)*_xyz->dX[2];
      } // if (j<Ntot)
    } // if (j3>=iz1 && j3<_temp->ilaserLoc ...
  } // for (int j...
  nGrain += _xyz->nX[0]*_xyz->nX[1]*_xyz->nZlayer;  
  // assign appropriate vState (including zeroing states above ilaserLoc
  iz1 = _xyz->nX[0]*_xyz->nX[1]*_temp->ilaserLoc;
  for (int j=0;j<Ntot2;++j){
    if (_part->icellidLoc[j] >=iz1){
      vState[j] = 0;
      gID[j] = 0;
    }    
  } // for j  
}; // AddLayer


void VoxelsCA::UpdateLayer(std::string &filCSV){
  /*
    Function adds powder to layer by making every voxel in layer a unique grain.
    It is called at the start of a new layer
    
  */

  // purge indices of cTheta associated with grains having no volume
  std::vector<double> gVol(nGrain,0.0),gVolT(nGrain,0),ctmp;
  std::vector<int> gtmp,itmp1(nGrain,0);
  int Ntot,Ntot2,nGtmp,i1t,nVlayer;
  double vvol;
  Ntot = _part->ncellLoc;
  Ntot2 = _part->ncellLoc + _part->nGhost;
  vvol = _xyz->dX[0]*_xyz->dX[1]*_xyz->dX[2];
  for (int j=0;j<Ntot;++j){
    if (gID[j]>0){gVol[gID[j]-1]+=vvol;}
  }
  MPI_Allreduce(&gVol[0],&gVolT[0],nGrain,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  nGtmp=0;
  ctmp.assign(4*nGrain,0);
  for (int j=0;j<nGrain;++j){
    if (gVolT[j]>0.0){
      ctmp[4*nGtmp] = cTheta[4*j];
      ctmp[4*nGtmp+1] = cTheta[4*j+1];
      ctmp[4*nGtmp+2] = cTheta[4*j+2];
      ctmp[4*nGtmp+3] = cTheta[4*j+3];
      itmp1[j] = nGtmp+1;
      nGtmp+=1;
    }
  }
  cTheta.assign(ctmp.begin(),ctmp.begin()+4*nGtmp);
  nGrain = nGtmp;
  for (int j=0;j<Ntot2;++j){
    if (gID[j]==0){continue;}
    i1t = itmp1[gID[j]-1];
    gID[j] =  i1t;
  }
  // create new grains for each voxel in layer
  WriteCSVData1(filCSV);
  nVlayer = _xyz->nX[0]*_xyz->nX[1]*_xyz->nZlayer;
  std::default_random_engine gen1(3134525);
  unsigned int sdloc;
  sdloc = unsigned(double(gen1())/double(gen1.max())*pow(2.0,32.0));
  SampleOrientation sa1;
  // randomly generate crystallographic orientations
  std::vector<double> aa;
  sa1.GenerateSamples(nVlayer,sdloc,aa);
  cTheta.insert(cTheta.end(), aa.begin(),aa.end());
  // ASSIGN GID NGRAIN1+1...NGRAIN TO PROPER VOXELS
  // LOOK INTO WHETHER GNUCLEUS IS IMPORTANT IN WHICH
  // CASE IT WILL BE THE VOXEL CENTERS
  int j1,j2,j3,iz1,jvox0;
  iz1 = _temp->ilaserLoc - _xyz->nZlayer;
  jvox0 = _xyz->nX[0]*_xyz->nX[1]*iz1;
  for (int j=0;j<Ntot2;++j){
    j3 = floor(_part->icellidLoc[j]/(_xyz->nX[0]*_xyz->nX[1]));
    if (j3>=iz1 && j3<_temp->ilaserLoc){
      j2 = floor( (_part->icellidLoc[j]- _xyz->nX[0]*_xyz->nX[1]*j3)/_xyz->nX[0]);
      j1 = _part->icellidLoc[j] - _xyz->nX[0]*_xyz->nX[1]*j3 - _xyz->nX[0]*j2;
      gID[j] = nGrain + _xyz->nX[0]*_xyz->nX[1]*j3+_xyz->nX[0]*j2+j1 - jvox0 + 1;
      vState[j] = 3;
      if (j<Ntot){
	centroidOct[3*j]=(double(j1)+.5)*_xyz->dX[0];
	centroidOct[3*j+1]=(double(j2)+.5)*_xyz->dX[1];
	centroidOct[3*j+2]=(double(j3)+.5)*_xyz->dX[2];
      } // if (j<Ntot)
    } // if (j3>=iz1 && j3<_temp->ilaserLoc ...
  } // for (int j...
  nGrain += _xyz->nX[0]*_xyz->nX[1]*_xyz->nZlayer;
  
  // assign appropriate vState (including zeroing states above ilaserLoc
  iz1 = _xyz->nX[0]*_xyz->nX[1]*_temp->ilaserLoc;
  for (int j=0;j<Ntot2;++j){
    if (_part->icellidLoc[j] >=iz1){
      vState[j] = 0;
      gID[j] = 0;
    }    
  } // for j  
}; // UpdateLayer

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
	if (j<_part->ncellLoc){extents[j] = extentsInitialValue;}
      }
    } 
  } // end for j
};

void VoxelsCA::SetLiquid3(){
  // makes cell liquid if temperature exceeds liquidus
  int Ntot = _part->ncellLoc + _part->nGhost,n1,j,nZlayer,iz1; 
  int Ntot1 = _part->ncellLoc;
  n1 = _xyz->nX[0]*_xyz->nX[1];
  nZlayer = round(_xyz->layerT/_xyz->dX[2]);
  iz1 = _temp->ilaserLoc*_xyz->nX[0]*_xyz->nX[1];
  for (int j=0;j<Ntot;++j){
    if (_temp->TempCurr[j] >= _xyz->tL ) { 
      if (_part->icellidLoc[j]<n1){
	vState[j]=2;
      } else if ( _part->icellidLoc[j]<iz1){
	vState[j] = 1;
	gID[j] = 0; // flag that it loses its grain
      } // if (_part->icellidLoc[j]<n1..
      if (j<Ntot1){extents[j] = extentsInitialValue;}
    } // if (_temp->TempCurr[j]...
  } // end for j
};
void VoxelsCA::SetLiquid4(){
  // makes cell liquid if temperature exceeds liquidus
  // this is for test case
  int Ntot = _part->ncellLoc + _part->nGhost,n1; 
  int Ntot1=_part->ncellLoc;
  int j;
  for (int j=0;j<Ntot;++j){
    n1 = std::distance(gNucleus.begin(),std::find(
       gNucleus.begin(),gNucleus.end(),_part->icellidLoc[j]));
    if (_temp->TempCurr[j] >= _xyz->tL && n1==nGrain) { 
	vState[j] = 1;
	gID[j] = 0; // flag that it loses its grain
	// if (j<Ntot1){extents[j] = extentsInitialValue;}
    } // if (
  } // end for j
};
void VoxelsCA::SetLiquid5(){
  // makes cell liquid if temperature exceeds liquidus
  int Ntot = _part->ncellLoc + _part->nGhost,n1,j,nZlayer,iz1;
  int Ntot1 = _part->ncellLoc;
  double tcheck = .99*_temp->tBeg0;
  n1 = _xyz->nX[0]*_xyz->nX[1];
  nZlayer = round(_xyz->layerT/_xyz->dX[2]);
  iz1 = _temp->ilaserLoc*_xyz->nX[0]*_xyz->nX[1];
  for (int j=0;j<Ntot;++j){
    if (_temp->tBeg[j] < tcheck) {
      if (_part->icellidLoc[j]<n1){
        vState[j]=2;
      } else if ( _part->icellidLoc[j]<iz1){
        vState[j] = 1;
        gID[j] = 0; // flag that it loses its grain
      } // if (_part->icellidLoc[j]<n1..
      if (j<Ntot1){extents[j] = extentsInitialValue;}
    } // if (_temp->tBeg[j]...
  } // end for j
};
void VoxelsCA::SetLiquid6(){
  // makes cell liquid if temperature exceeds liquidus
  int Ntot = _part->ncellLoc + _part->nGhost,n1,j,nZlayer,iz1;
  int Ntot1 = _part->ncellLoc;
  double tcheck1=1.01*_temp->tBeg0;
  n1 = _xyz->nX[0]*_xyz->nX[1];
  nZlayer = round(_xyz->layerT/_xyz->dX[2]);
  iz1 = _temp->ilaserLoc*_xyz->nX[0]*_xyz->nX[1];
  for (int j=0;j<Ntot;++j){
    if (_temp->tBeg[j] > tcheck1) {
      if (_part->icellidLoc[j]<n1){
        vState[j]=2;
      } else if ( _part->icellidLoc[j]<iz1){
        vState[j] = 1;
        //gID[j] = 0; // flag that it loses its grain
      } // if (_part->icellidLoc[j]<n1..
      if (j<Ntot1){extents[j] = extentsInitialValue;}
    } // if (_temp->tBeg[j]...
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
  n1 = _xyz->nX[0]*_xyz->nX[1]*_temp->ilaserLoc;
  for (int j=0;j<Ntot;++j){
    if (_part->icellidLoc[j] >=n1){
      vState[j] = 0;
      gID[j] = 0;
    } 
  } // for j
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
  std::vector<double> nrate(Ntot,0),gidtmp,cthtmp,tmpall1,tmp1,aa;
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
  aa.assign(ngtmp*4,0.0);
  SampleOrientation sa;
  sa.GenerateSamples(ngtmp,sdloc,aa);
  for (int j1=0;j1<ngtmp;++j1){
    cthtmp.push_back(aa[4*j1]);
    cthtmp.push_back(aa[4*j1+1]);
    cthtmp.push_back(aa[4*j1+2]);
    cthtmp.push_back(aa[4*j1+3]);
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
  float IPFmapBD[3*nC], IPFmapx[3*nC], IPFmapy[3*nC];
  double vBD[3]={0.0,0.0,1.0},omega,ax[3],vCD[3],mxAng,blue,green,red,rRot[3][3],mscale,
    vX[3]={1.0,0.0,0.0},vY[3]={0.0,1.0,0.0};
  for (int j=0;j<nC;++j){
    TempOut[j] = _temp->TempCurr[j];
    if (gID[j]<1){
      IPFmapBD[3*j] = 0.0;
      IPFmapBD[3*j+1] = 0.0;
      IPFmapBD[3*j+2] = 0.0;
      IPFmapx[3*j] = 0.0;
      IPFmapx[3*j+1] = 0.0;
      IPFmapx[3*j+2] = 0.0;
      IPFmapy[3*j] = 0.0;
      IPFmapy[3*j+1] = 0.0;
      IPFmapy[3*j+2] = 0.0;
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
      vCD[2] = std::min(vCD[2],1.0);
      mxAng = M_PI/4.0;
      red = std::fabs( (mxAng - acos(vCD[2]))/mxAng );
      blue = atan2(vCD[1],vCD[0]);
      green = mxAng - blue;
      blue *= (1-red)/mxAng;
      green *= (1-red)/mxAng;
      mscale = std::max(red,std::max(green,blue));
      IPFmapBD[3*j] = red/mscale;
      IPFmapBD[3*j+1] = green/mscale;
      IPFmapBD[3*j+2] = blue/mscale;
      // x dir
      vCD[0] = std::fabs(rRot[0][0]*vX[0]+rRot[1][0]*vX[1]+rRot[2][0]*vX[2]);
      vCD[1] = std::fabs(rRot[0][1]*vX[0]+rRot[1][1]*vX[1]+rRot[2][1]*vX[2]);
      vCD[2] = std::fabs(rRot[0][2]*vX[0]+rRot[1][2]*vX[1]+rRot[2][2]*vX[2]);
      std::sort(vCD,vCD+3);
      std::swap(vCD[0],vCD[1]);
      vCD[2]=std::min(vCD[2],1.0);
      mxAng = M_PI/4.0;
      red = std::fabs( (mxAng - acos(vCD[2]))/mxAng );
      blue = atan2(vCD[1],vCD[0]);
      green = mxAng - blue;
      blue *= (1-red)/mxAng;
      green *= (1-red)/mxAng;
      mscale = std::max(red,std::max(green,blue));
      IPFmapx[3*j] = red/mscale;
      IPFmapx[3*j+1] = green/mscale;
      IPFmapx[3*j+2] = blue/mscale;
      // y dir 
      vCD[0] = std::fabs(rRot[0][0]*vY[0]+rRot[1][0]*vY[1]+rRot[2][0]*vY[2]);
      vCD[1] = std::fabs(rRot[0][1]*vY[0]+rRot[1][1]*vY[1]+rRot[2][1]*vY[2]);
      vCD[2] = std::fabs(rRot[0][2]*vY[0]+rRot[1][2]*vY[1]+rRot[2][2]*vY[2]);
      std::sort(vCD,vCD+3);
      std::swap(vCD[0],vCD[1]);
      vCD[2]=std::min(vCD[2],1.0);
      mxAng = M_PI/4.0;
      red = std::fabs( (mxAng - acos(vCD[2]))/mxAng );
      blue = atan2(vCD[1],vCD[0]);
      green = mxAng - blue;
      blue *= (1-red)/mxAng;
      green *= (1-red)/mxAng;
      mscale = std::max(red,std::max(green,blue));
      IPFmapy[3*j] = red/mscale;
      IPFmapy[3*j+1] = green/mscale;
      IPFmapy[3*j+2] = blue/mscale;
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
    fp << "        <DataArray type='Float32' Name='IPFx' NumberOfComponents='3' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += 3*nC * sizeof (float) + sizeof (int);
    fp << "        <DataArray type='Float32' Name='IPFy' NumberOfComponents='3' format='appended' offset='" << Offset << "'/>" << std::endl;
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
      CellsTH=3*nC*sizeof(float), CellsTemp=nC*sizeof(float);
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
    fp.write(reinterpret_cast<const char *>(&CellsTH), 4);
    for (int i=0;i<3*nC;i++) fp.write(reinterpret_cast<const char *>(&IPFmapx[i]), sizeof (float));
    fp.write(reinterpret_cast<const char *>(&CellsTH), 4);
    for (int i=0;i<3*nC;i++) fp.write(reinterpret_cast<const char *>(&IPFmapy[i]), sizeof (float));
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
void VoxelsCA::WriteToHDF1(const std::string &filename)
{
  // writes gID, vState, cTheta per voxel
  int Ntot=_part->ncellLoc, i1,i2,i3,icase;
  std::string hdf5Filename = filename + ".h5";
  std::vector< float> TempOut(Ntot,0),IPFmapBD(3*Ntot,0), IPFmapx(3*Ntot,0), IPFmapy(3*Ntot,0),cth(4*nGrain,0);
  double vBD[3]={0.0,0.0,1.0},omega,ax[3],vCD[3],mxAng,blue,green,red,rRot[3][3],mscale,
    vX[3]={1.0,0.0,0.0},vY[3]={0.0,1.0,0.0},xp,yp,x0,y0,m,a,b,c,H,S,V,sMax,ff,p,q,t;
  std::vector<std::vector<double>> triPts(2,std::vector<double>(3,0));
  triPts[0][0]=0.0;
  triPts[0][1]=2./pow(2,.5)/(1.+1./pow(2,.5));
  triPts[0][2]=2./pow(3,.5)/(1.+1./pow(3,.5));
  triPts[1][0]=0.0;
  triPts[1][1]=0.0;
  triPts[1][2]=2./pow(3,.5)/(1.+1./pow(3,.5));
  m=tan(1./2.*atan2(triPts[1][2],triPts[0][2]));
  a=pow(pow(triPts[1][2]-triPts[1][1],2.)+pow(triPts[0][2]-triPts[0][1],2.),.5);
  b=pow(pow(triPts[0][1],2.)+pow(triPts[1][1],2.) ,0.5);
  c=pow(pow(triPts[0][2],2.)+pow(triPts[1][2],2.),0.5);
  y0=1./2.*pow((b+c-a)*(c+a-b)*(a+b-c)/(a+b+c),.5);
  x0=y0/m;
  sMax=pow(pow(x0,2.)+pow(y0,2.),.5);
  for (int j=0;j<Ntot;++j){
    TempOut[j] = _temp->TempCurr[j];
    if (gID[j]<1){
      IPFmapBD[3*j] = 0.0;
      IPFmapBD[3*j+1] = 0.0;
      IPFmapBD[3*j+2] = 0.0;
      IPFmapx[3*j] = 0.0;
      IPFmapx[3*j+1] = 0.0;
      IPFmapx[3*j+2] = 0.0;
      IPFmapy[3*j] = 0.0;
      IPFmapy[3*j+1] = 0.0;
      IPFmapy[3*j+2] = 0.0;
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
      vCD[2] = std::min(vCD[2],1.0);
      xp = 2.*vCD[0]/(1.+vCD[2]);
      yp = 2.*vCD[1]/(1.+vCD[2]);
      H=atan( (yp-y0)/(xp-x0))*180./M_PI;
      xp < x0 ? H+=180: H;
      H = H+240-atan((triPts[1][2]-y0)/(triPts[0][2]-x0))*180/M_PI;
      V=1.;
      S=pow(pow(xp-x0,2.)+pow(yp-y0,2.),.5);
      S=S/sMax*0.8 + 0.2;
      H>=360.0 ? H=0.0 : H;
      icase = floor(H/60.0);
      ff= H/6.0 - icase;
      p=V*(1.-S);
      q=V*(1.-S*ff);
      t=V*(1.-(S*(1.-ff)));
      if (S<=0.0){
	IPFmapBD[3*j]=0.0;
	IPFmapBD[3*j+1]=0.0;
	IPFmapBD[3*j+2]=0.0;
      } else {
	if (icase==0){
	  IPFmapBD[3*j]=V;
	  IPFmapBD[3*j+1]=t;
	  IPFmapBD[3*j+2]=p;
	}
	if (icase==1){
	  IPFmapBD[3*j]=q;
	  IPFmapBD[3*j+1]=V;
	  IPFmapBD[3*j+2]=p;
	}
	if (icase==2){
	  IPFmapBD[3*j]=p;
	  IPFmapBD[3*j+1]=V;
	  IPFmapBD[3*j+2]=t;
	}
	if (icase==3){
	  IPFmapBD[3*j]=p;
	  IPFmapBD[3*j+1]=q;
	  IPFmapBD[3*j+2]=V;
	}
	if (icase==4){
	  IPFmapBD[3*j]=t;
	  IPFmapBD[3*j+1]=p;
	  IPFmapBD[3*j+2]=V;
	}
	if (icase==5){
	  IPFmapBD[3*j]=V;
	  IPFmapBD[3*j+1]=p;
	  IPFmapBD[3*j+2]=q;
	}
      }
      // x dir
      vCD[0] = std::fabs(rRot[0][0]*vX[0]+rRot[1][0]*vX[1]+rRot[2][0]*vX[2]);
      vCD[1] = std::fabs(rRot[0][1]*vX[0]+rRot[1][1]*vX[1]+rRot[2][1]*vX[2]);
      vCD[2] = std::fabs(rRot[0][2]*vX[0]+rRot[1][2]*vX[1]+rRot[2][2]*vX[2]);
      std::sort(vCD,vCD+3);
      std::swap(vCD[0],vCD[1]);
      vCD[2]=std::min(vCD[2],1.0);
      xp = 2.*vCD[0]/(1.+vCD[2]);
      yp = 2.*vCD[1]/(1.+vCD[2]);
      H=atan( (yp-y0)/(xp-x0))*180./M_PI;
      xp < x0 ? H+=180: H;
      H = H+240-atan((triPts[1][2]-y0)/(triPts[0][2]-x0))*180/M_PI;
      V=1.;
      S=pow(pow(xp-x0,2.)+pow(yp-y0,2.),.5);
      S=S/sMax*0.8 + 0.2;
      H>=360.0 ? H=0.0 : H;
      icase = floor(H/60.0);
      ff= H/6.0 - icase;
      p=V*(1.-S);
      q=V*(1.-S*ff);
      t=V*(1.-(S*(1.-ff)));
      if (S<=0.0){
	IPFmapx[3*j]=0.0;
	IPFmapx[3*j+1]=0.0;
	IPFmapx[3*j+2]=0.0;
      } else {
	if (icase==0){
	  IPFmapx[3*j]=V;
	  IPFmapx[3*j+1]=t;
	  IPFmapx[3*j+2]=p;
	}
	if (icase==1){
	  IPFmapx[3*j]=q;
	  IPFmapx[3*j+1]=V;
	  IPFmapx[3*j+2]=p;
	}
	if (icase==2){
	  IPFmapx[3*j]=p;
	  IPFmapx[3*j+1]=V;
	  IPFmapx[3*j+2]=t;
	}
	if (icase==3){
	  IPFmapx[3*j]=p;
	  IPFmapx[3*j+1]=q;
	  IPFmapx[3*j+2]=V;
	}
	if (icase==4){
	  IPFmapx[3*j]=t;
	  IPFmapx[3*j+1]=p;
	  IPFmapx[3*j+2]=V;
	}
	if (icase==5){
	  IPFmapx[3*j]=V;
	  IPFmapx[3*j+1]=p;
	  IPFmapx[3*j+2]=q;
	}
      }
      // y dir 
      vCD[0] = std::fabs(rRot[0][0]*vY[0]+rRot[1][0]*vY[1]+rRot[2][0]*vY[2]);
      vCD[1] = std::fabs(rRot[0][1]*vY[0]+rRot[1][1]*vY[1]+rRot[2][1]*vY[2]);
      vCD[2] = std::fabs(rRot[0][2]*vY[0]+rRot[1][2]*vY[1]+rRot[2][2]*vY[2]);
      std::sort(vCD,vCD+3);
      std::swap(vCD[0],vCD[1]);
      vCD[2]=std::min(vCD[2],1.0);
      xp = 2.*vCD[0]/(1.+vCD[2]);
      yp = 2.*vCD[1]/(1.+vCD[2]);
      H = H+240-atan((triPts[1][2]-y0)/(triPts[0][2]-x0))*180/M_PI;
      V=1.;
      S=pow(pow(xp-x0,2.)+pow(yp-y0,2.),.5);
      S=S/sMax*0.8 + 0.2;
      H>=360.0 ? H=0.0 : H;
      icase = floor(H/60.0);
      ff= H/6.0 - icase;
      p=V*(1.-S);
      q=V*(1.-S*ff);
      t=V*(1.-(S*(1.-ff)));
      if (S<=0.0){
	IPFmapy[3*j]=0.0;
	IPFmapy[3*j+1]=0.0;
	IPFmapy[3*j+2]=0.0;
      } else {
	if (icase==0){
	  IPFmapy[3*j]=V;
	  IPFmapy[3*j+1]=t;
	  IPFmapy[3*j+2]=p;
	}
	if (icase==1){
	  IPFmapy[3*j]=q;
	  IPFmapy[3*j+1]=V;
	  IPFmapy[3*j+2]=p;
	}
	if (icase==2){
	  IPFmapy[3*j]=p;
	  IPFmapy[3*j+1]=V;
	  IPFmapy[3*j+2]=t;
	}
	if (icase==3){
	  IPFmapy[3*j]=p;
	  IPFmapy[3*j+1]=q;
	  IPFmapy[3*j+2]=V;
	}
	if (icase==4){
	  IPFmapy[3*j]=t;
	  IPFmapy[3*j+1]=p;
	  IPFmapy[3*j+2]=V;
	}
	if (icase==5){
	  IPFmapy[3*j]=V;
	  IPFmapy[3*j+1]=p;
	  IPFmapy[3*j+2]=q;
	}
      } // if (S<0.0...
    } //     if (gID[j]<1){
  } // for (int j..
  for (int j=0;j<4*nGrain;++j){cth[j]=cTheta[j];}
  unsigned int nVoxT, js,jc,ncth;
  nVoxT = _xyz->nX[0]*_xyz->nX[1]*_xyz->nX[2];
  ncth=4*nGrain;
  js = _part->icellidLoc[0];
  jc = _part->icellidLoc[Ntot-1]-js + 1;
  std::vector<float> dX(_xyz->dX.begin(),_xyz->dX.end());
  adios2::ADIOS adios(MPI_COMM_WORLD);
  adios2::IO hdf5IO = adios.DeclareIO("HDFFileIO");
  hdf5IO.SetEngine("HDF5");
  // global array : name, { shape (total) }, { start (local) }, { count (local) }  all are constant dimensions
  adios2::Variable<int> dimsa = hdf5IO.DefineVariable<int>(
	      "dims", {3}, {0}, {3});
  adios2::Variable<float> dxa = hdf5IO.DefineVariable<float>(
	      "VoxelDX", {3}, {0}, {3});
  adios2::Variable<int> gida = hdf5IO.DefineVariable<int>(
	      "gID", {nVoxT}, {js}, {jc});
  adios2::Variable<int> vStatea = hdf5IO.DefineVariable<int>(
	      "vState", {nVoxT}, {js}, {jc});
  adios2::Variable<float> TempOuta = hdf5IO.DefineVariable<float>(
	      "Temperature", {nVoxT}, {js}, {jc});
  adios2::Variable<float> IPFmapBDa = hdf5IO.DefineVariable<float>(
	      "IPFz", {3*nVoxT}, {3*js}, {3*jc});
  adios2::Variable<float> IPFmapxa = hdf5IO.DefineVariable<float>(
	      "IPFx", {3*nVoxT}, {3*js}, {3*jc});
  adios2::Variable<float> IPFmapya = hdf5IO.DefineVariable<float>(
	      "IPFy", {3*nVoxT}, {3*js}, {3*jc});
  adios2::Variable<float> angAx = hdf5IO.DefineVariable<float>(
	      "angleAxis", {ncth}, {0}, {ncth});
  adios2::Engine hdf5Writer =
      hdf5IO.Open(hdf5Filename, adios2::Mode::Write);
  hdf5Writer.Put<int>(dimsa, _xyz->nX.data());
  hdf5Writer.Put<float>(dxa, dX.data());
  hdf5Writer.Put<int>(gida, gID.data());
  hdf5Writer.Put<int>(vStatea, vState.data());
  hdf5Writer.Put<float>(TempOuta, TempOut.data());
  hdf5Writer.Put<float>(IPFmapBDa, IPFmapBD.data());
  hdf5Writer.Put<float>(IPFmapxa, IPFmapx.data());
  hdf5Writer.Put<float>(IPFmapya, IPFmapy.data());
  hdf5Writer.Put<float>(angAx, cth.data());
  hdf5Writer.Close();
  MPI_Barrier(MPI_COMM_WORLD);  
} // end WriteToHDF1

void VoxelsCA::WriteToVTU2(const std::string &filename)
{
  // writes gID, vState, cTheta per voxel
  std::string vtkFilename = filename + ".vtu";
  std::ofstream fp;
  int Offset=0, nC=_part->ncellLoc, nP=_part->npointLoc,j1,j2,j3,cell_offsets;
  unsigned char cell_type;
  std::vector< float> TempOut(nC,0);
  float IPFmapBD[3*nC];
  double vBD[3]={0.0,0.0,1.0},omega,ax[3],vCD[3],mxAng,blue,green,red,rRot[3][3],mscale;
  for (int j=0;j<nC;++j){
    TempOut[j] = _temp->tBeg[j];
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
      vCD[2] = std::min(vCD[2],1.0);
      mxAng = M_PI/4.0;
      red = std::fabs( (mxAng - acos(vCD[2]))/mxAng );
      blue = atan2(vCD[1],vCD[0]);
      green = mxAng - blue;
      blue *= (1-red)/mxAng;
      green *= (1-red)/mxAng;
      mscale = std::max(red,std::max(green,blue));
      IPFmapBD[3*j] = red/mscale;
      IPFmapBD[3*j+1] = green/mscale;
      IPFmapBD[3*j+2] = blue/mscale;
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
} // end WriteToVTU2

void VoxelsCA::WriteToPVTU1(const std::string &filename, int &nFile)
{
  /* 
     writes gID, vState, IPFx,IPFy,IPFz per voxel
     - if vtu file > 2G then need to break up into smaller files
       for paraview to read. This function is for that case. it
       produces the same data output as writeToVTU1.
   */
  int Offset=0, nC=_part->ncellLoc, nP=_part->npointLoc,j1,j2,j3,cell_offsets,rank;
  unsigned char cell_type;
  std::vector< float> TempOut(nC,0);
  float IPFmapBD[3*nC], IPFmapx[3*nC], IPFmapy[3*nC];
  double vBD[3]={0.0,0.0,1.0},omega,ax[3],vCD[3],mxAng,blue,green,red,rRot[3][3],mscale,
    vX[3]={1.0,0.0,0.0},vY[3]={0.0,1.0,0.0};
  for (int j=0;j<nC;++j){
    TempOut[j] = _temp->TempCurr[j];
    if (gID[j]<1){
      IPFmapBD[3*j] = 0.0;
      IPFmapBD[3*j+1] = 0.0;
      IPFmapBD[3*j+2] = 0.0;
      IPFmapx[3*j] = 0.0;
      IPFmapx[3*j+1] = 0.0;
      IPFmapx[3*j+2] = 0.0;
      IPFmapy[3*j] = 0.0;
      IPFmapy[3*j+1] = 0.0;
      IPFmapy[3*j+2] = 0.0;
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
      vCD[2] = std::min(vCD[2],1.0);
      mxAng = M_PI/4.0;
      red = std::fabs( (mxAng - acos(vCD[2]))/mxAng );
      blue = atan2(vCD[1],vCD[0]);
      green = mxAng - blue;
      blue *= (1-red)/mxAng;
      green *= (1-red)/mxAng;
      mscale = std::max(red,std::max(green,blue));
      IPFmapBD[3*j] = red/mscale;
      IPFmapBD[3*j+1] = green/mscale;
      IPFmapBD[3*j+2] = blue/mscale;
      // x dir
      vCD[0] = std::fabs(rRot[0][0]*vX[0]+rRot[1][0]*vX[1]+rRot[2][0]*vX[2]);
      vCD[1] = std::fabs(rRot[0][1]*vX[0]+rRot[1][1]*vX[1]+rRot[2][1]*vX[2]);
      vCD[2] = std::fabs(rRot[0][2]*vX[0]+rRot[1][2]*vX[1]+rRot[2][2]*vX[2]);
      std::sort(vCD,vCD+3);
      std::swap(vCD[0],vCD[1]);
      vCD[2]=std::min(vCD[2],1.0);
      mxAng = M_PI/4.0;
      red = std::fabs( (mxAng - acos(vCD[2]))/mxAng );
      blue = atan2(vCD[1],vCD[0]);
      green = mxAng - blue;
      blue *= (1-red)/mxAng;
      green *= (1-red)/mxAng;
      mscale = std::max(red,std::max(green,blue));
      IPFmapx[3*j] = red/mscale;
      IPFmapx[3*j+1] = green/mscale;
      IPFmapx[3*j+2] = blue/mscale;
      // y dir 
      vCD[0] = std::fabs(rRot[0][0]*vY[0]+rRot[1][0]*vY[1]+rRot[2][0]*vY[2]);
      vCD[1] = std::fabs(rRot[0][1]*vY[0]+rRot[1][1]*vY[1]+rRot[2][1]*vY[2]);
      vCD[2] = std::fabs(rRot[0][2]*vY[0]+rRot[1][2]*vY[1]+rRot[2][2]*vY[2]);
      std::sort(vCD,vCD+3);
      std::swap(vCD[0],vCD[1]);
      vCD[2]=std::min(vCD[2],1.0);
      mxAng = M_PI/4.0;
      red = std::fabs( (mxAng - acos(vCD[2]))/mxAng );
      blue = atan2(vCD[1],vCD[0]);
      green = mxAng - blue;
      blue *= (1-red)/mxAng;
      green *= (1-red)/mxAng;
      mscale = std::max(red,std::max(green,blue));
      IPFmapy[3*j] = red/mscale;
      IPFmapy[3*j+1] = green/mscale;
      IPFmapy[3*j+2] = blue/mscale;
    }
  }
  if (_xyz->nnodePerCell==4){cell_type=9;}
  if (_xyz->nnodePerCell==8){cell_type=12;}
  std::string vtuf,pvtuFilename = filename + ".pvtu";
  std::ofstream fp[nFile],f0;
  std::vector<std::string> vtuFilenames(nFile);
  for (int f = 0; f < nFile; f++) {
    std::stringstream ftmp;
    ftmp << filename << "_f" << f << ".vtu";
    ftmp >> vtuf;
    vtuFilenames[f] = vtuf;
  }  
  
  
  if (_part->myid == 0 ) {
    f0.open(pvtuFilename.c_str());
    f0 << "<?xml version='1.0'?>" << std::endl;
    f0 << "<VTKFile type='PUnstructuredGrid' version='0.1' byte_order='LittleEndian'>" << std::endl;
    f0 << "  <PUnstructuredGrid GhostLevel='0'>" << std::endl;
    f0 << "      <PPoints>" << std::endl;
    f0 << "        <PDataArray type='Float32' Name='Position' NumberOfComponents='3'/>" << std::endl;
    f0 << "      </PPoints>" << std::endl;
    f0 << "      <PCells>" << std::endl;
    f0 << "        <PDataArray type='Int32' Name='connectivity' NumberOfComponents='1'/>" << std::endl;
    f0 << "        <PDataArray type='Int32' Name='offsets' NumberOfComponents='1'/>" << std::endl;
    f0 << "        <PDataArray type='UInt8' Name='types' NumberOfComponents='1'/>" << std::endl;
    f0 << "      </PCells>" << std::endl;
    f0 << "      <PPointData>" << std::endl;
    f0 << "      </PPointData>" << std::endl;
    f0 << "      <PCellData>" << std::endl;
    f0 << "        <PDataArray type='Int32' Name='vState' NumberOfComponents='1'/>" << std::endl;
    f0 << "        <PDataArray type='Int32' Name='gID' NumberOfComponents='1'/>" << std::endl;
    f0 << "        <PDataArray type='Float32' Name='IPFz' NumberOfComponents='3'/>" << std::endl;
    f0 << "        <PDataArray type='Float32' Name='IPFx' NumberOfComponents='3'/>" << std::endl;
    f0 << "        <PDataArray type='Float32' Name='IPFy' NumberOfComponents='3'/>" << std::endl;
    f0 << "        <PDataArray type='Float32' Name='Temperature' NumberOfComponents='1'/>" << std::endl;
    f0 << "      </PCellData>" << std::endl;
    for (int f=0; f<nFile; ++f) {
      f0 << "    <Piece Source='" << vtuFilenames[f] << "'/>" << std::endl;
    }
    f0 << "  </PUnstructuredGrid>" << std::endl;
    f0 << "</VTKFile>" << std::endl;
    f0.close();
  } // if (_part->myid==0 ...

  if (_part->myid == 0) {
    for (int f = 0; f < nFile; f++) {
      vtuf = vtuFilenames[f];
      fp[f].open(vtuf.c_str());
      fp[f] << "<?xml version='1.0'?>" << std::endl;
      fp[f] << "<VTKFile type='UnstructuredGrid' version='0.1' byte_order='LittleEndian'>" << std::endl;
      fp[f] << "  <UnstructuredGrid>" << std::endl;
      fp[f] << "    <FieldData>" << std::endl;
      fp[f] << "      <DataArray type='Int32' Name='time_step' NumberOfTuples='1' format='ascii'>" << std::endl;
      fp[f] << "        " << _xyz->tInd << std::endl;
      fp[f] << "      </DataArray>" << std::endl;
      fp[f] << "    </FieldData>" << std::endl;
      fp[f].close();
    } // for (int f ...
  } // if (_part->myid...


  int nPIncr = int(_part->nprocs/nFile);
  for (int f = 0; f < nFile; f++) {
    Offset = 0;
    vtuf = vtuFilenames[f];
    for (int j = 0; j < nPIncr; j++) {
      rank = f*nPIncr + j;
      if (rank>=_part->nprocs) break;
      int OffsetP=Offset;
      MPI_Allreduce(&OffsetP,&Offset,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    if (_part->myid != rank) { Offset=0; continue; }
    fp[f].open(vtuf.c_str(), std::fstream::app);
    fp[f] << "    <Piece NumberOfPoints='" << nP << "' NumberOfCells='" << nC << "'>" << std::endl;
    fp[f] << "      <Points>" << std::endl;
    fp[f] << "        <DataArray type='Float32' Name='Position' NumberOfComponents='3' format='appended' offset='" << Offset << "'/>" << std::endl;
    fp[f] << "      </Points>" << std::endl;
    Offset += 3 * nP * sizeof(float) + sizeof(int);
    fp[f] << "      <Cells>" << std::endl;
    fp[f] << "        <DataArray type='Int32' Name='connectivity' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += _part->iconnectivityLoc.size() * sizeof (int) + sizeof(int);
    fp[f] << "        <DataArray type='Int32' Name='offsets' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += nC * sizeof (int) + sizeof(int);
    fp[f] << "        <DataArray type='UInt8' Name='types' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += nC * sizeof (unsigned char) + sizeof(int);
    fp[f] << "      </Cells>" << std::endl;
    fp[f] << "      <PointData>" << std::endl;
    fp[f] << "      </PointData>" << std::endl;
    fp[f] << "      <CellData>" << std::endl;
    fp[f] << "        <DataArray type='Int32' Name='vState' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += nC * sizeof (int) + sizeof (int);
    fp[f] << "        <DataArray type='Int32' Name='gID' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += nC * sizeof (int) + sizeof (int);
    fp[f] << "        <DataArray type='Float32' Name='IPFz' NumberOfComponents='3' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += 3*nC * sizeof (float) + sizeof (int);
    fp[f] << "        <DataArray type='Float32' Name='IPFx' NumberOfComponents='3' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += 3*nC * sizeof (float) + sizeof (int);
    fp[f] << "        <DataArray type='Float32' Name='IPFy' NumberOfComponents='3' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += 3*nC * sizeof (float) + sizeof (int);
    fp[f] << "        <DataArray type='Float32' Name='Temperature' NumberOfComponents='1' format='appended' offset='" << Offset << "'/>" << std::endl;
    Offset += nC * sizeof (float) + sizeof (int);
    fp[f] << "      </CellData>" << std::endl;
    fp[f] << "    </Piece>" << std::endl;
    fp[f].close();
    } // for (int j ..
  } // for (int f ...
  MPI_Barrier(MPI_COMM_WORLD);
  if (_part->myid == 0) {
    for (int f=0;f<nFile;++f){
      vtuf = vtuFilenames[f];
      fp[f].open(vtuf.c_str(), std::fstream::app);
      fp[f] << "  </UnstructuredGrid>" << std::endl;
      fp[f] << "  <AppendedData encoding='raw'>" << std::endl;
      fp[f] << "_";
      fp[f].close();
    }
  }
  for (int f = 0; f < nFile; f++) {
    Offset = 0;
    vtuf = vtuFilenames[f];
    for (int jf = 0; jf < nPIncr; jf++) {
      rank = f*nPIncr + jf;
      if (rank>=_part->nprocs) break;
      MPI_Barrier(MPI_COMM_WORLD);
      if (_part->myid != rank) continue;
      fp[f].open(vtuf.c_str(), std::fstream::app);
      int Scalar = nP * sizeof (float), Vector = 3 * Scalar, Cells = nC * sizeof(int), 
	CellsTH=3*nC*sizeof(float), CellsTemp=nC*sizeof(float);
      int CellChars = nC * sizeof(unsigned char), Conn = _part->iconnectivityLoc.size() * sizeof(int);
      fp[f].write(reinterpret_cast<const char *>(&Vector), 4);
      for (int j=0;j<nP;j++) {
	j3 = floor(_part->ipointidLoc[j]/( (_xyz->nX[0]+1)*(_xyz->nX[1]+1)));
	j2 = floor( (_part->ipointidLoc[j] - (_xyz->nX[0]+1)*(_xyz->nX[1]+1)*j3)/(_xyz->nX[0]+1));
	j1 = _part->ipointidLoc[j] - (_xyz->nX[0]+1)*(_xyz->nX[1]+1)*j3 - (_xyz->nX[0]+1)*j2;
	float x = j1*_xyz->dX[0], y = j2*_xyz->dX[1], z = j3*_xyz->dX[2];
	fp[f].write(reinterpret_cast<const char *>(&x), sizeof(float));
	fp[f].write(reinterpret_cast<const char *>(&y), sizeof(float));
	fp[f].write(reinterpret_cast<const char *>(&z), sizeof(float));
      }
      fp[f].write(reinterpret_cast<const char *>(&Conn), 4);
      for (int i=0;i<nC;++i) {
	for (int j=0;j<_xyz->nnodePerCell;++j) {
	  j1 = _part->iconnectivityLoc[_xyz->nnodePerCell*i+j];
	  fp[f].write(reinterpret_cast<const char *>(&j1), sizeof(int));
	}
      }
      fp[f].write(reinterpret_cast<const char *>(&Cells), 4);
      cell_offsets = 0;
      for (int i=0;i<nC;i++) {
	cell_offsets += _xyz->nnodePerCell;
	fp[f].write(reinterpret_cast<const char *>(&cell_offsets), sizeof(int));
      }
      fp[f].write(reinterpret_cast<const char *>(&CellChars), 4);    
      for (int i=0;i<nC;i++) fp[f].write(reinterpret_cast<const char *>(&cell_type), sizeof (unsigned char));
      fp[f].write(reinterpret_cast<const char *>(&Cells), 4);
      for (int i=0;i<nC;i++) fp[f].write(reinterpret_cast<const char *>(&vState[i]), sizeof (int));
      fp[f].write(reinterpret_cast<const char *>(&Cells), 4);
      for (int i=0;i<nC;i++) fp[f].write(reinterpret_cast<const char *>(&gID[i]), sizeof (int));
      fp[f].write(reinterpret_cast<const char *>(&CellsTH), 4);
      for (int i=0;i<3*nC;i++) fp[f].write(reinterpret_cast<const char *>(&IPFmapBD[i]), sizeof (float));
      fp[f].write(reinterpret_cast<const char *>(&CellsTH), 4);
      for (int i=0;i<3*nC;i++) fp[f].write(reinterpret_cast<const char *>(&IPFmapx[i]), sizeof (float));
      fp[f].write(reinterpret_cast<const char *>(&CellsTH), 4);
      for (int i=0;i<3*nC;i++) fp[f].write(reinterpret_cast<const char *>(&IPFmapy[i]), sizeof (float));
      fp[f].write(reinterpret_cast<const char *>(&CellsTemp), 4);
      for (int i=0;i<nC;i++) fp[f].write(reinterpret_cast<const char *>(&TempOut[i]), sizeof (float));
      fp[f].close();
    } // for (int jf ...
  } // for (int f ...
  MPI_Barrier(MPI_COMM_WORLD);
  if (_part->myid == 0) {
    for (int f=0;f<nFile;++f){
      vtuf = vtuFilenames[f];
      fp[f].open(vtuf.c_str(), std::fstream::app);
      fp[f] << "  </AppendedData>" << std::endl;
      fp[f] << "</VTKFile>" << std::endl;
      fp[f].close();
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);  
} // end WriteToPVTU1


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
    fp << "axis-angle (omega,n), grain Volume" << std::endl;
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

void VoxelsCA::WriteCSVData1(const std::string &filename)
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
      fp << cTheta[4*j] << "," << cTheta[4*j+1] << "," << cTheta[4*j+2] << "," << 
	cTheta[4*j+3]<<","<<gVolT[j] << std::endl;
    } // end for (int j...
    fp.close();
  } // if (_part->myid ..)
  MPI_Barrier(MPI_COMM_WORLD);
} // WriteCSVData1

void VoxelsCA::NucleateGrains(std::vector<int> &nIn, std::vector<double> &tIn)
{
  int Ntot=_part->ncellLoc,Nnew1,Nq=0,NqA,Ntot2=_part->ncellLoc+_part->nGhost,cc=0,cc1;
  std::vector<int> ind1(Ntot);
  double tNp1= _xyz->time + _temp->DelT,x1,x2,rmax=0.0,rateX=0.0;
  std::vector<double> TempNp1,rNuc(Ntot,0);
  std::default_random_engine g1(30*_xyz->tInd+seed1);
  _temp->AnalyticTempCurr(tNp1,TempNp1,_part->icellidLoc,Ntot2);
  for (int j=0;j<Ntot;++j){
    if (vState[j]==1 && _temp->TempCurr[j]< _xyz->tL){
      x2 = ( (_xyz->tL -TempNp1[j]) - _xyz->dTempM)/_xyz->dTempS/pow(2.0,.5);
      x1 = ( (_xyz->tL - _temp->TempCurr[j]) - _xyz->dTempM)/_xyz->dTempS/pow(2.0,.5);
      rNuc[j] = _xyz->rNmax*std::max(0.0, .5*(erf(x2) - erf(x1)));
      rmax = std::max(rmax,rNuc[j]);
      ind1[cc] = j;
      cc+=1;
    } // if (vState[j]==1...
  } // for (int j=0...
  Nq=cc;
  MPI_Allreduce(&cc,&NqA,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  //rateX= std::min(double(NqA),(double)NqA *_xyz->dX[0]*_xyz->dX[1]*_xyz->dX[2] * rmax);
  rateX= (double)NqA *_xyz->dX[0]*_xyz->dX[1]*_xyz->dX[2] * rmax;
  // sample Poisson distribution from rmax and then do thinning
  std::poisson_distribution<int> Np(rateX);
  std::uniform_int_distribution<int> uind(0,std::max(Nq-1,0));
  std::uniform_real_distribution<double> xrand1(0.0,1.0);
  Nnew1 = (Nq==0)? 0:  int(floor( (double(Nq)/ double(NqA))* Np(g1)));
  Nnew1 = std::min(Nnew1,Nq);
  std::vector<int> itmp1(Nnew1),itmp2(Nq,0);
  std::vector<double> tmp1(Nnew1);
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
  // perform thinning process
  cc=0;
  for (int j=0;j<Nnew1;++j){
    if (xrand1(g1) > (rNuc[itmp1[j]]/rmax)){continue;}
    nIn[cc] =itmp1[j];
    tIn[cc] = tmp1[j];
    cc+=1;
  } // for (int j ...
  nIn.resize(cc);
  tIn.resize(cc);
} // end NucleateGrains

void VoxelsCA::NucleateGrains(std::vector<int> &nucA, std::vector<double> &tnucA,
			      std::vector<double> &quatnuc,std::vector<int> &iv,
			      std::vector<int> &vI,std::vector<int> &vS,std::vector<double> &tA)
{

  int Ntot=iv[_part->myid+1]-iv[_part->myid],Nnew1,Nq=0,NqA,cc=0,cc1;
  std::vector<int> ind1(Ntot),nIn;
  double tNp1= _xyz->time + _temp->DelT,x1,x2,rmaxloc=0.0,rmax,rateX=0.0;
  std::vector<double> rNuc(Ntot,0),tIn;
  std::default_random_engine g1(30*_xyz->tInd+seed1);
  cc1=iv[_part->myid];
  cc1=0;
  for (int j=iv[_part->myid];j<iv[_part->myid+1];++j){
    if (vS[j]==1 && tA[j]<_xyz->tL){
      //x2 = ( (_xyz->tL -TempNp1[cc1]) - _xyz->dTempM)/_xyz->dTempS/pow(2.0,.5);
      x1 = ( (_xyz->tL - tA[j]) - _xyz->dTempM)/_xyz->dTempS/pow(2.0,.5);
      //rNuc[cc1] = _xyz->rNmax*std::max(0.0, .5*(erf(x2) - erf(x1)));
      //rNuc[cc1] = std::max(0.0, .5*erf(x1)+.5 );
      rNuc[cc1] = exp(-.5*pow(x1,2.0));
      rmaxloc = std::max(rmaxloc,rNuc[cc1]);
      ind1[cc] = j;
      cc+=1;
    } // if (vS[j]...
    cc1+=1;
  } // for (int j=0...
  Nq=cc;


  MPI_Allreduce(&Nq,&NqA,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&rmaxloc,&rmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  //rateX= std::min(double(NqA),(double)NqA *_xyz->dX[0]*_xyz->dX[1]*_xyz->dX[2] * rmax);
  rateX= double(NqA) * _xyz->rNmax* rmax;
  // sample Poisson distribution from rmax and then do thinning
  std::poisson_distribution<int> Np(rateX);
  std::uniform_int_distribution<int> uind(0,std::max(Nq-1,0));
  std::uniform_real_distribution<double> xrand1(0.0,1.0);
  Nnew1 = (Nq==0)? 0:  int(floor( (double(Nq)/ double(NqA))* Np(g1)));
  Nnew1 = std::min(Nnew1,Nq);
  std::vector<int> itmp1(Nnew1),itmp2(Nq,0);
  std::vector<double> tmp1(Nnew1);
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
  // perform thinning process
  cc=0;
  for (int j=0;j<Nnew1;++j){
    if (xrand1(g1) > (rNuc[itmp1[j]]/rmax)){continue;}
    nIn[cc] =itmp1[j];
    tIn[cc] = tmp1[j];
    cc+=1;
  } // for (int j ...
  nIn.resize(cc);
  tIn.resize(cc);
  Nq=cc;
  int NnucA=0;
  std::vector<int> Nnucvec(_part->nprocs,0),jnuc0(_part->nprocs,0);
  MPI_Allgather(&Nq,1,MPI_INT,&Nnucvec[0],1,MPI_INT,MPI_COMM_WORLD);
  for (int j=0;j<_part->nprocs;++j){NnucA+=Nnucvec[j];}
  jnuc0[0]=0;
  for (int j=1;j<_part->nprocs;++j){jnuc0[j] = jnuc0[j-1]+Nnucvec[j-1];}
  nucA.resize(NnucA);
  tnucA.resize(NnucA);
  std::vector<int> ivec(NnucA),itmp(NnucA);
  std::iota(ivec.begin(),ivec.end(),0);
  std::vector<double> tmp(NnucA);
  cc=jnuc0[_part->myid];
  for (int j=0;j<Nq;++j){
    itmp[cc+j] = nIn[j];
    tmp[cc+j] = tIn[j];
  }
  for (int j=0;j<_part->nprocs;++j){
    if (Nnucvec[j]==0){continue;}
    MPI_Bcast(&itmp[jnuc0[j]],Nnucvec[j],MPI_INT,j,MPI_COMM_WORLD);
    MPI_Bcast(&tmp[jnuc0[j]],Nnucvec[j],MPI_DOUBLE,j,MPI_COMM_WORLD);
  } // for (int j ..    
  std::sort(ivec.begin(),ivec.end(),[&tmp](int j1,int j2){return tmp[j1]<tmp[j2];});
  for (int j=0;j<NnucA;++j){
    cc=ivec[j];
    tnucA[j]=tmp[cc];
    nucA[j]=itmp[cc];
  }
  unsigned int sdloc;
  SampleOrientation sa;
  sdloc= seed0 + 32*_xyz->tInd +64*nGrain;
  quatnuc.resize(4*NnucA);
  sa.GenerateSamples(NnucA,sdloc,quatnuc);

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
