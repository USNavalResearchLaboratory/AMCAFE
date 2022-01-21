// define member functions of VoxelCA

#include "Grid.cuh"
#include "VoxelsCA.cuh"
#include "iostream"
#include "fstream"
#include "sstream"
#include <math.h>
#include <algorithm>
#include <random>
#include <ctime>
#include "SampleOrientation.cuh"
#include <chrono>
#include <adios2.h>
#include <curand_kernel.h>
#include <math_constants.h>

static void HandleError(cudaError_t err) {
    if (err != cudaSuccess) {
        printf("%s \n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

__global__ void getSites(Grid *g,VoxelsCA *vx,double *xs, int numPG)
{
// calculates sites for Voronoi (done in own global function to avoid 
// inter-block race condition)
  int tid=threadIdx.x+blockDim.x*blockIdx.x, subsq=0,iz1,js,
    stride=blockDim.x*gridDim.x;
  unsigned int seedL = (vx->seed0+tid*100+g->tInd*10) % 4294967295;
  curandState_t s1;
  curand_init(seedL,subsq,0,&s1);
  double Lx=g->nX[0]*g->dX[0], Ly=g->nX[1]*g->dX[1],
    lt=g->layerT,zloc0;
  iz1 = g->ilaserLoc - g->nZlayer;
  zloc0 = iz1*g->dX[0];
  js=tid;
  while (js <numPG){
    xs[3*js]= curand_uniform(&s1)*Lx;
    xs[3*js+1]= curand_uniform(&s1)*Ly;
    xs[3*js+2]= curand_uniform(&s1)*lt+zloc0;
    js+=stride;
  }
}
__global__ void addLayer1Part1(Grid *g,VoxelsCA *vx,double *xs, double *troids,
			  int *gD, int *vs, int *itmp, int numPG)
{
  int tid=threadIdx.x+blockDim.x*blockIdx.x,iz1,js,
    stride=blockDim.x*gridDim.x,j1,j2,j3,nX1=g->nX[0],
    iz2=g->ilaserLoc,i3=g->nX[0]*g->nX[1],i2,ng1=vx->nGrain;
  double dx=g->dX[0],dsqc,dsq,xc,yc,zc;
  iz1 = g->ilaserLoc - g->nZlayer;
  js=tid+i3*iz1;
  while (js < i3*iz2){
    j3 = js/i3; // note that int/int is equivalent to floor                                                                         
    j2 = (js - i3*j3)/(nX1);
    j1 = js - i3*j3 - nX1*j2;
    zc = (j3+.5)*dx;;
    yc = (j2+.5)*dx;
    xc = (j1+.5)*dx;
    dsqc=1e6;
    for (int jz=0;jz<numPG;++jz){
      dsq = pow(xc - xs[3*jz],2.)+pow(yc-xs[3*jz+1],2.)+pow(zc-xs[3*jz+2],2.);
      if (dsq<dsqc){i2=jz; dsqc=dsq;}
    } // for (int jz
    gD[js] = i2+1+ng1;
    itmp[i2]=1;
    vs[js] = 3;
    troids[3*js] = xc;
    troids[3*js+1] = yc;
    troids[3*js+2] = zc;
    js+=stride;
  } // while (js < i3*...
} // end addlayer1part1

__global__ void addLayer1Part2(Grid *g,VoxelsCA *vx, double *dctheta,
			  int *gD, int *itmp, int numPG)
{

  // !!***************************************************
  // FUNCTION CAN ONLY RUN WITH 1 BLOCK. OTHERWISE, BUG
  // !!***************************************************

  int tid=threadIdx.x+blockDim.x*blockIdx.x, subsq=0,iz1,js,
    stride=blockDim.x*gridDim.x,ng1=vx->nGrain;
  unsigned int seedL = (vx->seed0+tid*100+g->tInd*10) % 4294967295;
  curandState_t s1;
  curand_init(seedL,subsq,0,&s1);
  double lt=g->layerT;
  iz1 = g->ilaserLoc - g->nZlayer;
  int i3=g->nX[0]*g->nX[1],i2,i1;
  // ensures a continuous numbering of grain ids
  __shared__ int i2s;
  if (tid==0){
    i2=0;
    for (int j=0;j<numPG;++j){
      if (itmp[j]==1){
	itmp[j]=i2;
	i2+=1;
      }
    }
    i2s=i2;
    vx->nGrain += i2+1;
  } // if (tid==0
  __syncthreads();
  if (i2s != numPG){
    js=tid+i3*iz1;
    while (js < i3*(lt+iz1)){
      i1=gD[js]-1;
      gD[js]=itmp[i1]+1+ng1;
      js+=stride;
    }
  }
  // randomly generate crystallographic orientations 
  double axAng[4];
  i2s+=1;
  js=tid;
  while (js<i2s){
    GenerateSamples(1,seedL,subsq,s1, axAng);
    dctheta[4*(js+ng1)]=axAng[0];
    dctheta[4*(js+ng1)+1]=axAng[1];
    dctheta[4*(js+ng1)+2]=axAng[2];
    dctheta[4*(js+ng1)+3]=axAng[3];
    js+=stride;
  }
} // __global__ void addlayer1part2

__global__ void addLayer1Part3(const Grid *g,int *gD, int *vs) 
{
  int tid=threadIdx.x+blockDim.x*blockIdx.x,js,i1,ntot,stride;
  stride=blockDim.x*gridDim.x;
  i1=g->nX[0]*g->nX[1]*g->ilaserLoc;
  js=tid+i1;
  ntot=g->nX[0]*g->nX[1]*g->nX[2];
  while (js<ntot){
    vs[js]=0;
    gD[js]=0;
    js+=stride;
  }
} // __global__ void addlayer1par3


__global__ void copyGlobal(double *x1,double *x0, int n)
{
  int tid=threadIdx.x + blockDim.x*blockIdx.x,js,stride;
  js=tid;
  stride=blockDim.x*gridDim.x;
  while (js < n){
    x1[js] = x0[js];
    js +=stride;
  }
}

__global__ void cleanLayerPart1(VoxelsCA *dvx,int *dgid,int *gvolflg, int *itmp,int Ntot)
{
  int tid=threadIdx.x+blockIdx.x*blockDim.x,js,i1,stride;
  stride = blockDim.x*gridDim.x;
  js=tid;
  while (js < Ntot){
    i1=dgid[js]-1;
    gvolflg[i1]=1;
    js+=stride;
  }
}
__global__ void cleanLayerPart2(VoxelsCA *dvx,int *gvolflg, int *itmp)
{
  int ng=dvx->nGrain,tid=threadIdx.x+blockIdx.x*blockDim.x,i2;
  if (tid==0){
    i2=0;
    for (int j=0;j<ng;++j){
      if (gvolflg[j]>0){
	itmp[j] = i2+1;
	i2+=1;
      }
    }
    gvolflg[ng]=i2;
  } 
}

__global__ void cleanLayerPart3(VoxelsCA *dvx,int *dgid,int *gvolflg, int *itmp,
				double *ctmp, double *dcth,int Ntot)
{
  int ng,tid=threadIdx.x+blockIdx.x*blockDim.x,js,stride,i1;
  stride=blockDim.x*gridDim.x;
  js=tid;
  ng=dvx->nGrain;
  while(js<ng){
    if (gvolflg[js]>0){
      ctmp[4*(itmp[js]-1)] = dcth[4*js];
      ctmp[4*(itmp[js]-1)+1] = dcth[4*js+1];
      ctmp[4*(itmp[js]-1)+2] = dcth[4*js+2];
      ctmp[4*(itmp[js]-1)+3] = dcth[4*js+3];
    }
    js+=stride;
  }
  js=tid;
  while (js<Ntot){
    if (dgid[js]!=0){
      i1 = itmp[dgid[js]-1];
      dgid[js]=i1;
    }
    js+=stride;
  }
}
__global__ void cleanLayerPart4(VoxelsCA *dvx, int *gvolflg)
{
  int tid=threadIdx.x+blockDim.x*blockIdx.x;
  int ng = dvx->nGrain;
  if (tid==0){dvx->nGrain = gvolflg[ng];}
}
__global__ void convertSolid1Part1(Grid *dg,int *vst,double *exts,int ntot)
{
  // ONLY RUN WITH 1 BLOCK
  int inei[27],vnei[26],js,nb,
    tid=threadIdx.x,jz,stride=blockDim.x;
  bool issol;
  extern __shared__ bool i1[];
  nb=ntot/stride + 1;
  for (int j=0;j<nb;++j){
    js = tid+j*stride;
    __syncthreads();
    i1[tid]=0;
    if (js<ntot){
      if (vst[js]==2){
	dg->GetNeighbors(js,inei);
	for (int j1=0;j1<inei[26];++j1){vnei[j1]=vst[inei[j1]];}
	jz=0;
	issol=1;
	while (issol && (jz<inei[26])){
	  issol = issol && (vnei[jz]>=2);
	  jz+=1;
	}
	if (issol){
	  i1[threadIdx.x]=1;
	  exts[js]=0.0;
	}
      }
    }
    __syncthreads();
    if (i1[tid]){vst[js]=3;}
  }
}
__global__ void convertSolid1Part2(Grid *dg,int *vst,int ntot)
{
  // ONLY RUN WITH 1 BLOCK 
  int inei[27],vnei[26],js,nb,
    tid=threadIdx.x,stride=blockDim.x;
  extern __shared__ bool i1[];
  nb=ntot/stride+1;
  for (int j=0;j<nb;++j){
    js=tid+j*stride;
    __syncthreads();
    i1[threadIdx.x]=0;
    inei[26]=0;
    if (js<ntot){
      if (vst[js]==1){
	dg->GetNeighbors(js,inei);
	for (int j1=0;j1<inei[26];++j1){vnei[j1]=vst[inei[j1]];}
      }
    }
    __syncthreads();
    for (int j1=0;j1<inei[26];++j1){
      if (vnei[j1]==3){vst[inei[j1]]=2;}
    }
  }
}
__global__ void setLiquid3(Grid *dg,int *dgid,int *dvstate,double *dtempval, double *dexts)
{
  // makes cell liquid if temperature exceeds liquidus
  int n1=dg->nX[0]*dg->nX[1],ntot,iz1,js,
    tid=threadIdx.x+blockIdx.x*blockDim.x,stride=blockDim.x*gridDim.x,
    tL=dg->tL;
  ntot=n1*dg->nX[2];
  iz1=dg->ilaserLoc*n1;
  js=tid;
  while (js<ntot){
    if (dtempval[js]>=tL){
      if (js < n1){
	dvstate[js]=2;
      } else if (js < iz1){
	dvstate[js]=1;
	dgid[js]=0; 
	dexts[js] = 0.0;
      }
    }
    js+=stride;
  }
}
__global__ void updateVoxelsPart1(int *vstate,int *vs2cc,int ntot)
{
  int tid=threadIdx.x+blockIdx.x*blockDim.x,js,stride=blockDim.x*gridDim.x,s,
    tidL=threadIdx.x,bidL=blockIdx.x,nthread=blockDim.x;
  extern __shared__ volatile int sh[];
  js=tid;
  sh[tidL]=0;
  while (js<ntot){
    if (vstate[js]==2){sh[tidL]+=1;}
    js+=stride;
  }
  __syncthreads();
  s=nthread;
  while (s>=128){
    if (nthread >= s){
      if (tidL<s/2){sh[tidL] += sh[tidL+s/2];}
      __syncthreads();
    }
    s/=2;
  }
  if (tidL<32) {
    if (nthread >=64) {sh[tidL] += sh[tidL+32];}
    if (nthread >=32) {sh[tidL] += sh[tidL+16];}
    if (nthread >=16) {sh[tidL] += sh[tidL+8];}
    if (nthread >= 8) {sh[tidL] += sh[tidL+4];}
    if (nthread >= 4) {sh[tidL] += sh[tidL+2];}
    if (nthread >= 2) {sh[tidL] += sh[tidL+1];}
  }
  if (tidL==0) {
    vs2cc[bidL] = sh[tidL];
  }
}

__global__ void updateVoxelsPart2(int *vs2cc, int *vs2, bool *disf,double *dtinc)
{
  vs2[0]=vs2cc[0];
  vs2[1]=0;
  dtinc[0]=0.0;
  disf[0]=1; // initialize to true (false exits voxel growth loop)
}

__global__ void updateVoxelsPart3(Grid *dg,int *dgid,int *dvstate, double *dctheta,double *dtempval,
				  double *dctroid,double *dexts,double *dtinc,int *vs2,bool *disf,
				  float *dtminG,int *jindG,int *j1indG,int ntot)
{
  int tid=threadIdx.x,gid=threadIdx.x+blockIdx.x*blockDim.x,inei[27],js,
    stride=gridDim.x*blockDim.x,jx[3],timeUntil,dr,i3=dg->nX[0]*dg->nX[1],
    i2=dg->nX[0],s,nthread=blockDim.x;
  double vhat,dnx[3],dlocX[3],omega,ax[3],rRot[3][3],th,ph,tmelt=dg->tL,
    Avel=dg->Avel,nvel=dg->nvel,dx=dg->dX[0];
  bool isaxt,isaxv;
  extern __shared__ volatile int sh[];
  volatile int *jvec=sh;
  volatile int *j1vec=&jvec[nthread];
  volatile float *dtminv=(float*)&j1vec[nthread];
  disf[0]=(vs2[0]!=vs2[1] && dtinc[0]<dg->bmDelT);
  if (disf[0]){
    dtminv[tid]=1e6;
    if (tid==0){vs2[1]=vs2[0];}
    js=gid;
    while (js<ntot){
      if (dvstate[js]==2){
	dg->GetNeighbors(js,inei);
	isaxt=1;
	isaxv=1;
	for (int j=0;j<inei[26];++j){
	  isaxt=isaxt && dtempval[inei[j]]<tmelt;
	  isaxv=isaxv && dvstate[inei[j]]!=1;
	}
      	if ( (!isaxv) && isaxt){
	  omega=dctheta[4*(dgid[js]-1)];
	  ax[0]=dctheta[4*(dgid[js]-1)+1];
	  ax[1]=dctheta[4*(dgid[js]-1)+2];
	  ax[2]=dctheta[4*(dgid[js]-1)+3];
	  loadRotMat(omega,ax,rRot);
	  dtempval[js]>=tmelt ? vhat=0. : vhat = Avel*pow(tmelt-dtempval[js],nvel);
          for (int j=0;j<inei[26];++j){
            if (dvstate[inei[j]] != 1 ){continue;}
            jx[2] = inei[j]/i3;
            jx[1] = (inei[j]- i3*jx[2])/i2;
            jx[0] = inei[j] - i3*jx[2] - i2*jx[1];      
            dnx[0] = (double(jx[0])+.5)*dx - dctroid[3*js];
            dnx[1] = (double(jx[1])+.5)*dx - dctroid[3*js+1];
            dnx[2] = (double(jx[2])+.5)*dx - dctroid[3*js+2];
            th = atan2(fabs(dnx[1]),fabs(dnx[0]));
            th > CUDART_PI/4.0 ? th= CUDART_PI/2.0 - th: th ;
            ph = atan2(pow(pow(dnx[0],2.0)+pow(dnx[1],2.0),.5),fabs(dnx[2]));
            ph < CUDART_PI/4.0 ? ph = CUDART_PI/2.0 - ph: ph ;
            // matrix is local->global; need to multiply by transpose for global->local
            // put into 1st quadrant b/c of symmetry
            dlocX[0] = fabs(rRot[0][0]*dnx[0]+rRot[1][0]*dnx[1]+rRot[2][0]*dnx[2]);
            dlocX[1] = fabs(rRot[0][1]*dnx[0]+rRot[1][1]*dnx[1]+rRot[2][1]*dnx[2]);
            dlocX[2] = fabs(rRot[0][2]*dnx[0]+rRot[1][2]*dnx[1]+rRot[2][2]*dnx[2]);
            dr = pow(cos(th)*sin(ph),.5)*(dlocX[0]+dlocX[1]+dlocX[2]) - dexts[js];
            timeUntil = dr/vhat;
	    if (timeUntil<dtminv[tid]){
	      dtminv[tid] = timeUntil;
	      jvec[tid] = js;
	      j1vec[tid] = inei[j];
	    }
	  } // for (int j=0
	} // if ((!isaxv) && isaxt...
      }	// if (dvstate[js]==2...
      js+=stride;
    } // while (js<ntot
    // reduce within block    
    __syncthreads();
    s=nthread;
    while (s>=128){
      if (nthread >= s){
	if (tid<s/2){
	  if (dtminv[tid]>dtminv[tid+s/2]){
	    dtminv[tid] = dtminv[tid+s/2];
	    jvec[tid] = jvec[tid+s/2];
	    j1vec[tid] = j1vec[tid+s/2];
	  }
	}
	__syncthreads();
      }
      s/=2;
    }
    if (tid<32) {
      if (nthread >=64) {
	if (dtminv[tid]> dtminv[tid+32]){
	  dtminv[tid] = dtminv[tid+32];
	  jvec[tid] = jvec[tid+32];
	  j1vec[tid] = j1vec[tid+32];
	}
      }
      if (nthread >=32) {
	if (dtminv[tid]> dtminv[tid+16]){
	  dtminv[tid] = dtminv[tid+16];
	  jvec[tid] = jvec[tid+16];
	  j1vec[tid] = j1vec[tid+16];
	}
      }
      if (nthread >=16) {
	if (dtminv[tid]> dtminv[tid+8]){
	  dtminv[tid] = dtminv[tid+8];
	  jvec[tid] = jvec[tid+8];
	  j1vec[tid] = j1vec[tid+8];
	}
      }
      if (nthread >= 8) {
	if (dtminv[tid]> dtminv[tid+4]){
	  dtminv[tid] = dtminv[tid+4];
	  jvec[tid] = jvec[tid+4];
	  j1vec[tid] = j1vec[tid+4];
	}
      }
      if (nthread >= 4) {
	if (dtminv[tid]> dtminv[tid+2]){
	  dtminv[tid] = dtminv[tid+2];
	  jvec[tid] = jvec[tid+2];
	  j1vec[tid] = j1vec[tid+2];
	}
      }
      if (nthread >= 2) {
	if (dtminv[tid]> dtminv[tid+1]){
	  dtminv[tid] = dtminv[tid+1];
	  jvec[tid] = jvec[tid+1];
	  j1vec[tid] = j1vec[tid+1];
	}
      }
    }  // if (tid<32) ...
    if (tid==0) {
      dtminG[blockIdx.x] = dtminv[tid];
      jindG[blockIdx.x] = jvec[tid];
      j1indG[blockIdx.x] = j1vec[tid];
    }
  } // if (disf[0])
}

__global__ void reduceVoxelCapture(float *dtminG,int *jindG, int *j1indG, int n)
{
  int js,tid=threadIdx.x+blockDim.x*blockIdx.x,stride=blockDim.x*gridDim.x,s,
    nthread=blockDim.x,tidL=threadIdx.x,bidL=blockIdx.x;
  extern __shared__ volatile int sh[];
  volatile int *jvec=sh;
  volatile int *j1vec=&jvec[nthread];
  volatile float *dtminv=(float*)&j1vec[nthread];
  js=tid;
  dtminv[tidL]=1e6;
  while (js<n){
    if (dtminv[tidL]>dtminG[js]){
      dtminv[tidL] = dtminG[js];
      jvec[tidL] = jindG[js];
      j1vec[tidL] = j1indG[js];
    }      
    js+=stride;
  }
  __syncthreads();
  s=nthread;
  while (s>=128){
    if (nthread >= s){
      if (tidL<s/2){
	if (dtminv[tidL]>dtminv[tidL+s/2]){
	  dtminv[tidL] = dtminv[tidL+s/2];
	  jvec[tidL] = jvec[tidL+s/2];
	  j1vec[tidL] = j1vec[tidL+s/2];
	}
      }
      __syncthreads();
    }
    s/=2;
  }
  if (tidL<32) {
    if (nthread >=64) {
      if (dtminv[tidL]>dtminv[tidL+32]){
	dtminv[tidL] = dtminv[tidL+32];
	jvec[tidL] = jvec[tidL+32];
	j1vec[tidL] = j1vec[tidL+32];
      }
    }
    if (nthread >=32) {
      if (dtminv[tidL]>dtminv[tidL+16]){
	dtminv[tidL] = dtminv[tidL+16];
	jvec[tidL] = jvec[tidL+16];
	j1vec[tidL] = j1vec[tidL+16];
      }
    }
    if (nthread >=16) {
      if (dtminv[tidL]>dtminv[tidL+8]){
	dtminv[tidL] = dtminv[tidL+8];
	jvec[tidL] = jvec[tidL+8];
	j1vec[tidL] = j1vec[tidL+8];
      }
    }
    if (nthread >= 8) {
      if (dtminv[tidL]>dtminv[tidL+4]){
	dtminv[tidL] = dtminv[tidL+4];
	jvec[tidL] = jvec[tidL+4];
	j1vec[tidL] = j1vec[tidL+4];
      }
    }
    if (nthread >= 4) {
      if (dtminv[tidL]>dtminv[tidL+2]){
	dtminv[tidL] = dtminv[tidL+2];
	jvec[tidL] = jvec[tidL+2];
	j1vec[tidL] = j1vec[tidL+2];
      }
    }
    if (nthread >= 2) {
      if (dtminv[tidL]>dtminv[tidL+1]){
	dtminv[tidL] = dtminv[tidL+1];
	jvec[tidL] = jvec[tidL+1];
	j1vec[tidL] = j1vec[tidL+1];
      }
    }
  } // if (tid<32 ...
  if (tidL==0) {
    dtminG[bidL] = dtminv[tidL];
    jindG[bidL] = jvec[tidL];
    j1indG[bidL] = j1vec[tidL];
  }
}

__global__ void updateVoxelsPart4(Grid *dg,VoxelsCA *dvox,int *dgid,int *dvstate, double *dctheta,
				  double *dctroid,double *dexts,int *vs2,
				  float *dtminG,int *j1indG)
{
  // ONLY RUN WITH 1 THREAD: SERIAL PROCESS
  // this is for nucleation
  int subsq=0,js1,nsamp=1,ng=dvox->nGrain,jx[3],
    i3=dg->nX[0]*dg->nX[1],i2=dg->nX[0];
  unsigned int seedL= (dg->tInd*30+dvox->seed0+64*ng) % 4294967295;
  curandState_t s1;
  double dtmin,axAng[4],dx=dg->dX[0];
  dtmin=dtminG[0]; 
  __syncthreads();
  if (dtmin<1e6){
    curand_init(seedL,subsq,0,&s1);
    js1=j1indG[0];
    GenerateSamples(nsamp,seedL,subsq,s1,axAng);
    dvstate[js1] = 2;
    vs2[1]+=1;
    dvox->nGrain+=1;
    dctheta[4*ng]=axAng[0];
    dctheta[4*ng+1]=axAng[1];
    dctheta[4*ng+2]=axAng[2];
    dctheta[4*ng+3]=axAng[3];
    dgid[js1]=ng+1;
    jx[2]=js1/i3;
    jx[1]=(js1-i3*jx[2])/i2;
    jx[0]=js1-i3*jx[2]-jx[1];
    dctroid[3*js1]=(double(jx[0])+.5)*dx;
    dctroid[3*js1+1]=(double(jx[1])+.5)*dx;
    dctroid[3*js1+2]=(double(jx[2])+.5)*dx;
    dexts[js1]=0.0;
  }
}
__global__ void updateVoxelsPart5(Grid *dg,int *dvstate, double *dtempval,
				  double *dexts,float *dtminG,int ntot)
{
  // UPDATES EXTENTS 
  int tid=threadIdx.x+blockIdx.x*blockDim.x,js,
    stride=blockDim.x*gridDim.x,inei[27];
  double dtmin,Avel=dg->Avel,nvel=dg->nvel,vhat,tmelt=dg->tL;
  bool isaxt,isaxv;
  dtmin=dtminG[0];
  __syncthreads();
  if (dtmin<1e6){
    js=tid;
    while (js<ntot){
      vhat=0.;
      if (dvstate[js]==2){
	dg->GetNeighbors(js,inei);
	isaxt=1;
	isaxv=1;
	for (int j=0;j<inei[26];++j){
	  isaxt=isaxt && dtempval[inei[j]]<tmelt;
	  isaxv=isaxv && dvstate[inei[j]]!=1;
	}
      	if ( (!isaxv) && isaxt){
	  dtempval[js]>=tmelt ? vhat=0. : vhat = Avel*pow(tmelt-dtempval[js],nvel);
	}	
      }
      dexts[js]+=vhat*dtmin;
      js+=stride;
    }
  }
}

__global__ void updateVoxelsPart6(Grid *dg,VoxelsCA *dvox,int *dgid,int *dvstate, double *dctheta,
				  double *dctroid,double *dexts,int *vs2,double *dtinc,
				  float *dtminG,int *jindG,int *j1indG)
{
  // capture voxel **ONLY USE 1 THREAD - this is serial operation**
  int js,js1,jInd,i3=dg->nX[0]*dg->nX[1],
    i2=dg->nX[0],jx[3],sInd[8][3],jy[3];
  double xI[3],xJ[3],d1I,dI2,d1J,dJ3,L12,L13,l,Lmud,dnorm[3],xiL,
    dx=dg->dX[0],rRot[3][3],dr,th,ph,sdiag0[6][3],sdiag[6][3],dtmin,
    dnx[3],omega,ax[3],locX[3];
      
  dtmin=dtminG[0];
  if (dtmin<1e6){
    js=jindG[0];
    js1=j1indG[0];
    dvstate[js1]=2;
    vs2[1]+=1;
    dgid[js1]=dgid[js];
    jx[2] = js1/i3;
    jx[1] = (js1- i3*jx[2])/i2;
    jx[0] = js1 -i3*jx[2] - i2*jx[1];
    jy[2] = js/i3;
    jy[1] = (js- i3*jy[2])/i2;
    jy[0] = js - i3*jy[2] - i2*jy[1];
    l = pow( pow((jx[0]-jy[0])*dx,2)+ pow((jx[1]-jy[1])*dx,2)+
	     pow((jx[2]-jy[2])*dx,2),.5);
      dnx[0] = (double(jx[0])+.5)*dx - dctroid[3*js];
      dnx[1] = (double(jx[1])+.5)*dx - dctroid[3*js+1];
      dnx[2] = (double(jx[2])+.5)*dx - dctroid[3*js+2];
      omega = dctheta[4*(dgid[js]-1)];
      ax[0] = dctheta[4*(dgid[js]-1)+1];
      ax[1] = dctheta[4*(dgid[js]-1)+2];
      ax[2] = dctheta[4*(dgid[js]-1)+3];
      loadRotMat(omega,ax,rRot);    
      // matrix is local->global; need to multiply by transpose for global->local
      locX[0] = rRot[0][0]*dnx[0]+rRot[1][0]*dnx[1]+rRot[2][0]*dnx[2];
      locX[1] = rRot[0][1]*dnx[0]+rRot[1][1]*dnx[1]+rRot[2][1]*dnx[2];
      locX[2] = rRot[0][2]*dnx[0]+rRot[1][2]*dnx[1]+rRot[2][2]*dnx[2];
      th = atan2(fabs(dnx[1]),fabs(dnx[0]));
      th > CUDART_PI/4.0 ? th= CUDART_PI/2.0 - th: th;
      ph = atan2(pow(pow(dnx[0],2.0)+pow(dnx[1],2.0),.5),fabs(dnx[2]));
      ph < CUDART_PI/4.0 ? ph = CUDART_PI/2.0 - ph: ph;
      // signbit returns 0 if positive and 1 if negative
      dr = fabs(locX[0])+fabs(locX[1])+fabs(locX[2]);
      loadS(sdiag0,sInd); 
      for (int j1=0;j1<6;++j1){
	sdiag[j1][0]=sdiag0[j1][0]*dr;
	sdiag[j1][1]=sdiag0[j1][1]*dr;
	sdiag[j1][2]=sdiag0[j1][2]*dr;
      }
      if (locX[2]<0){jx[2]=1;} else{jx[2]=0;}
      if (locX[1]<0){jx[1]=1;} else{jx[1]=0;}
      if (locX[0]<0){jx[0]=1;} else{jx[0]=0;}
      jInd = 4*jx[2]+ 2*jx[1]+ jx[0];
      for (int j1=0;j1<3;++j1){
        dnorm[j1] = pow(locX[0]-sdiag[sInd[jInd][j1]][0],2.0)+
          pow(locX[1]-sdiag[sInd[jInd][j1]][1],2.0)+pow(locX[2]-sdiag[sInd[jInd][j1]][2],2.0);
      } // for (int j1...                                                                                                                                     
      jy[0]=0; jy[2]=0;
      if (dnorm[0]<=dnorm[1] && dnorm[0]<=dnorm[2]){jy[0]=0;}
      if (dnorm[1]<=dnorm[0] && dnorm[1]<=dnorm[2]){jy[0]=1;}
      if (dnorm[2]<=dnorm[0] && dnorm[2]<=dnorm[1]){jy[0]=2;}
      if (dnorm[0]>=dnorm[1] && dnorm[0]>=dnorm[2]){jy[2]=0;}
      if (dnorm[1]>=dnorm[0] && dnorm[1]>=dnorm[2]){jy[2]=1;}
      if (dnorm[2]>=dnorm[0] && dnorm[2]>=dnorm[1]){jy[2]=2;}
      if (jy[0]==jy[2]){
          jy[0]=0; jy[1]=1; jy[2]=2;
      } else {
        jy[1] = 3 - jy[0] - jy[2];
      }
      projectPointLine(locX,&(sdiag[sInd[jInd][jy[0]]][0]),&(sdiag[sInd[jInd][jy[1]]][0]),xI);
      projectPointLine(locX,&(sdiag[sInd[jInd][jy[0]]][0]),&(sdiag[sInd[jInd][jy[2]]][0]),xJ);
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
      L12 = .5*(fmin(d1I,pow(3.0,.5)*l) + fmin(dI2,pow(3.0,.5)*l) );
      L13 = .5*(fmin(d1J,pow(3.0,.5)*l) + fmin(dJ3,pow(3.0,.5)*l) );
      Lmud =  pow(3.0,.5)* (pow(2.0/3.0,.5)*fmax(L12,L13));
      xiL = 1.0;
      dexts[js1] = pow(cos(th)*sin(ph),.5)*Lmud*xiL;
      dnx[0] = sdiag[sInd[jInd][jy[0]]][0] - Lmud*sdiag0[sInd[jInd][jy[0]]][0];
      dnx[1] = sdiag[sInd[jInd][jy[0]]][1] - Lmud*sdiag0[sInd[jInd][jy[0]]][1];
      dnx[2] = sdiag[sInd[jInd][jy[0]]][2] - Lmud*sdiag0[sInd[jInd][jy[0]]][2];
      locX[0] = rRot[0][0]*dnx[0]+rRot[0][1]*dnx[1]+rRot[0][2]*dnx[2];
      locX[1] = rRot[1][0]*dnx[0]+rRot[1][1]*dnx[1]+rRot[1][2]*dnx[2];
      locX[2] = rRot[2][0]*dnx[0]+rRot[2][1]*dnx[1]+rRot[2][2]*dnx[2];
      dctroid[3*js1] = dctroid[3*js] + locX[0];
      dctroid[3*js1+1] = dctroid[3*js+1] + locX[1];
      dctroid[3*js1+2] = dctroid[3*js+2] + locX[2];
      dtinc[0]+=dtmin;
  }
}

__device__ void loadRotMat(double omega, double *ax, double  rRot[][3])
{
  // loads the rotation matrix from (omega,ax), note that
  // ax is a 3x1  and rRot is a 3x3 static arrays
  rRot[0][0] = cos(omega) + pow(ax[0],2.0)*(1-cos(omega));
  rRot[0][1] = ax[0]*ax[1]*(1-cos(omega)) - ax[2]*sin(omega);
  rRot[0][2] = ax[0]*ax[2]*(1-cos(omega)) + ax[1]*sin(omega);
  rRot[1][0] = ax[0]*ax[1]*(1-cos(omega)) + ax[2]*sin(omega);
  rRot[1][1] = cos(omega) + pow(ax[1],2.0)*(1-cos(omega));
  rRot[1][2] = ax[1]*ax[2]*(1-cos(omega)) - ax[0]*sin(omega);
  rRot[2][0] = ax[2]*ax[0]*(1-cos(omega)) - ax[1]*sin(omega);
  rRot[2][1] = ax[2]*ax[1]*(1-cos(omega)) + ax[0]*sin(omega);
  rRot[2][2] = cos(omega) + pow(ax[2],2.0)*(1-cos(omega));
} 

__device__ void projectPointLine(double *A, double *x0, double *x1, double *xproj)
{
  // computes the coordinate of projecting a point A on line connecting x0 and x1 = xproj
  // note that A,x0,x1,xproj are all length 3 arrays
  double n[3],snorm,d,t;
  snorm = pow(pow(x1[0]-x0[0],2)+pow(x1[1]-x0[1],2)+pow(x1[2]-x0[2],2),.5);
  n[0] = (x1[0]-x0[0])/snorm;
  n[1] = (x1[1]-x0[1])/snorm;
  n[2] = (x1[2]-x0[2])/snorm;
  d = -(n[0]*A[0]+n[1]*A[1]+n[2]*A[2]);
  t = (-d - n[0]*x0[0] - n[1]*x0[1] - n[2]*x0[2])/
    (n[0]*(x1[0]-x0[0])+n[1]*(x1[1]-x0[1])+n[2]*(x1[2]-x0[2]));
  xproj[0] = x0[0] + (x1[0]-x0[0])*t;
  xproj[1] = x0[1] + (x1[1]-x0[1])*t;
  xproj[2] = x0[2] + (x1[2]-x0[2])*t;
} //end projectPointLine...

__device__  void loadS(double S[][3],int sInd[][3])
{
  // this is for decentered octahedron method:
  // S is 6 corners of octahedron in local coor and sInd gives the 3 corner
  // indices for a given octant
  S[0][0]=1;  S[0][1]=0;  S[0][2]=0;
  S[1][0]=0;  S[1][1]=1;  S[1][2]=0;
  S[2][0]=0;  S[2][1]=0;  S[2][2]=1;
  S[3][0]=-1;  S[3][1]=0;  S[3][2]=0;
  S[4][0]=0;  S[4][1]=-1;  S[4][2]=0;
  S[5][0]=0;  S[5][1]=0;  S[5][2]=-1;
  sInd[0][0]=0;  sInd[0][1]=1;  sInd[0][2]=2;
  sInd[1][0]=1;  sInd[1][1]=2;  sInd[1][2]=3;
  sInd[2][0]=0;  sInd[2][1]=2;  sInd[2][2]=4;
  sInd[3][0]=2;  sInd[3][1]=3;  sInd[3][2]=4;
  sInd[4][0]=0;  sInd[4][1]=1;  sInd[4][2]=5;
  sInd[5][0]=1;  sInd[5][1]=3;  sInd[5][2]=5;
  sInd[6][0]=0;  sInd[6][1]=4;  sInd[6][2]=5;
  sInd[7][0]=3;  sInd[7][1]=4;  sInd[7][2]=5;
}// end inline void loadS

__global__ void reduceGlobalArray(int *ig, int n,int isw)
{
  int js,tid=threadIdx.x+blockDim.x*blockIdx.x,stride=blockDim.x*gridDim.x,s,
    nthread=blockDim.x,tidL=threadIdx.x,bidL=blockIdx.x;
  extern __shared__ volatile int sh[];
  if (isw==0){
    // isw=0 is for summation
    js=tid;
    sh[tidL]=0;
    while (js<n){
      sh[tidL]+=ig[js];
      js+=stride;
    }
    __syncthreads();
    s=nthread;
    while (s>=128){
      if (nthread >= s){
	if (tidL<s/2){sh[tidL] += sh[tidL+s/2];}
	__syncthreads();
      }
      s/=2;
    }
    if (tidL<32) {
      if (nthread >=64) {sh[tidL] += sh[tidL+32];}
      if (nthread >=32) {sh[tidL] += sh[tidL+16];}
      if (nthread >=16) {sh[tidL] += sh[tidL+8];}
      if (nthread >= 8) {sh[tidL] += sh[tidL+4];}
      if (nthread >= 4) {sh[tidL] += sh[tidL+2];}
      if (nthread >= 2) {sh[tidL] += sh[tidL+1];}
    }
    if (tidL==0) {
      ig[bidL] = sh[tidL];
    }
  }
}

void resizeGlobalArray(double **y, int &n0, int &n1)
{
  double *d_ctmp;
  HandleError(cudaMallocManaged((void**)&d_ctmp,n0*sizeof(double)));
  int nThreads=512;
  int nBlocks=n0/nThreads;
  copyGlobal<<<nBlocks,nThreads>>>(d_ctmp,*y, n0);
  cudaFree(*y);
  HandleError(cudaMallocManaged((void**)y,n1*sizeof(double)));
  copyGlobal<<<nBlocks,nThreads>>>(*y,d_ctmp, n0);
  cudaFree(d_ctmp);
}

void resizeArray(double **y, int &n)
{
  double *tmp;
  tmp=(double*)malloc(n*sizeof(double));
  free(*y);
  *y=tmp;
  tmp=NULL;
  free(tmp);

}

// constructor
VoxelsCA::VoxelsCA(Grid &g)
{ 
  int Ntot = g.nX[0]*g.nX[1]*g.nX[2];
  gID = (int*)malloc(Ntot*sizeof(int));
  vState = (int*)malloc(Ntot*sizeof(int));
  extents = (double*)malloc(Ntot*sizeof(double));
  centroidOct = (double*)malloc(3*Ntot*sizeof(double));
  memset(gID,0,Ntot*sizeof(int));
  memset(vState,0,Ntot*sizeof(int));
  memset(extents,0,Ntot*sizeof(double));
  memset(centroidOct,0,3*Ntot*sizeof(double));
  seed0= 2132512;
  seed1=2912351;
  genlayer.seed(seed1);
} // end constructor
void VoxelsCA::AddLayer1Macro(VoxelsCA *d_vx,Grid &g,Grid *d_g,double **d_cthptr,
		    double *d_troids, int *d_gid, int *d_vst,int &nbuf2)
{
  int npg, nThreads,nBlocks,nbuf1, *d_itmp,Ntot;
  double *d_Sites;
  Ntot=g.nX[0]*g.nX[1]*g.nX[2];
  getNumPowderGrains(g,npg);
  HandleError(cudaMallocManaged((void**)&d_Sites,npg*3*sizeof(double)));
  HandleError(cudaMallocManaged((void**)&d_itmp,npg*sizeof(int)));
  // below is buffer for size of cTheta to account for nucleation: 100*expected # of new grains in 3 layers
  nbuf2 = 4*(nGrain+ npg + int(ceil(g.nX[0]*g.nX[1]*(g.layerT/g.dX[2])*3*g.rNmax*pow(g.dX[0]*1e6,3.)*100)));
  nbuf1=4*(nGrain);
  resizeGlobalArray(d_cthptr,nbuf1,nbuf2);
  cudaDeviceSynchronize();
  nThreads=512;
  nBlocks = npg/nThreads;
  getSites<<<nBlocks,nThreads>>>(d_g,d_vx,d_Sites,npg);
  nBlocks=(g.nZlayer*g.nX[0]*g.nX[1])/nThreads;
  addLayer1Part1<<<nBlocks,nThreads>>>(d_g,d_vx,d_Sites,d_troids,d_gid,d_vst,d_itmp,npg);
  HandleError( cudaPeekAtLastError() );
  nBlocks=1;
  nThreads=256;
  addLayer1Part2<<<1,nThreads>>>(d_g,d_vx,*d_cthptr,d_gid,d_itmp,npg);
  HandleError( cudaPeekAtLastError() );
  nThreads=512;
  nBlocks=Ntot/nThreads;
  addLayer1Part3<<<nBlocks,nThreads>>>(d_g,d_gid,d_vst);
  HandleError( cudaPeekAtLastError() );
  HandleError(cudaMemcpy(&(nGrain), &(d_vx->nGrain), sizeof(int), cudaMemcpyDeviceToHost));
  cudaFree(d_Sites);
  cudaFree(d_itmp);
}
void VoxelsCA::CleanLayerMacro(VoxelsCA *dvx,int *dgid,double **dcthetaptr, int &nTot)
{
  int nThreads,nBlocks, *dgidtmp, *dgvflg;
  double *dcthtmp;
  HandleError(cudaMalloc((void**)&dgidtmp,nGrain*sizeof(int)));
  HandleError(cudaMalloc((void**)&dgvflg,(nGrain+1)*sizeof(int)));
  HandleError(cudaMemset(dgidtmp,0,nGrain*sizeof(int)));
  HandleError(cudaMemset(dgvflg,0,(nGrain+1)*sizeof(int)));
  HandleError(cudaMalloc((void**)&dcthtmp,4*(nGrain)*sizeof(double)));
  nThreads=512;
  nBlocks=nGrain/nThreads+1;
  cleanLayerPart1<<<nBlocks,nThreads>>>(dvx,dgid,dgvflg,dgidtmp,nTot);
  cleanLayerPart2<<<nBlocks,nThreads>>>(dvx,dgvflg,dgidtmp);
  cleanLayerPart3<<<nBlocks,nThreads>>>(dvx,dgid,dgvflg,dgidtmp,
					dcthtmp, *dcthetaptr,nTot);
  cleanLayerPart4<<<nBlocks,nThreads>>>(dvx,dgvflg);  
  cudaFree(*dcthetaptr);
  HandleError(cudaMemcpy(&(nGrain), &(dvx->nGrain), sizeof(int), cudaMemcpyDeviceToHost));
  HandleError(cudaMalloc((void**)dcthetaptr,4*(nGrain)*sizeof(double)));
  copyGlobal<<<nBlocks,nThreads>>>(*dcthetaptr,dcthtmp, 4*nGrain);
  cudaFree(dgidtmp);
  cudaFree(dgvflg);
  cudaFree(dcthtmp);
}
void VoxelsCA::ConvertSolid1Macro(Grid *dg,int *vstate,double *dextents,const int &iswitch,int nThreads,const int &ntot)
{
    if (iswitch==0){
    // converts 2 to 3 if all neighbors >=2 (ie mushy to solid)
    // can only use 1 block
    convertSolid1Part1<<<1,nThreads,nThreads*sizeof(int)>>>(dg,vstate,dextents,ntot);
  } else {
    // converts 3 to 2 if any neighbor is 1 (ie solid to mushy)
    // can only use 1 block
    convertSolid1Part2<<<1,nThreads>>>(dg,vstate,ntot);
  }
}
void VoxelsCA::SetLiquid3Macro(Grid *dg,int *dgid, int *dvstate,double *dtempval, double *dexts,int &nThreads, int &nBlocks)
{
  setLiquid3<<<nBlocks,nThreads>>>(dg,dgid,dvstate,dtempval,dexts);
}
void VoxelsCA::UpdateVoxelsMacro(Grid *dg, Grid &gg,VoxelsCA *dvox,int *dgid, int *dvstate,double *dtempval,double *dexts, 
			    double *troids, double *dctheta,int &nThreads, int &nBlocks, int &ntot)
{
  nThreads=1024;
  nBlocks=ntot/nThreads+1;
  SetLiquid3Macro(dg,dgid,dvstate,dtempval,dexts,nThreads,nBlocks);
  int isw=1, *dvs2c, *djindG,*dj1indG,n1;
  ConvertSolid1Macro(dg,dvstate,dexts,isw,nThreads,ntot);
  // calculate # of mushy voxels (vstate=2)
  bool isf=1,*d_isf;
  double *dtinc;
  float *ddtminG;
  nBlocks= ntot/nThreads+1;
  n1=nBlocks;
  HandleError(cudaMalloc((void**)&dtinc,sizeof(double)));
  HandleError(cudaMalloc((void**)&djindG,n1*sizeof(int)));
  HandleError(cudaMalloc((void**)&dvs2c,2*sizeof(int)));
  HandleError(cudaMalloc((void**)&d_isf,sizeof(bool)));
  updateVoxelsPart1<<<nBlocks,nThreads,nThreads*sizeof(int)>>>(dvstate,djindG,ntot);
  HandleError( cudaPeekAtLastError() );
  while (nBlocks>3){
    n1=nBlocks;
    nBlocks/=2;
    reduceGlobalArray<<<nBlocks,nThreads,nThreads*sizeof(int)>>>(djindG,n1,0);
  }
  n1=nBlocks;
  reduceGlobalArray<<<1,nThreads,nThreads*sizeof(int)>>>(djindG,n1,0);
  // initialize while loop
  updateVoxelsPart2<<<1,1>>>(djindG,dvs2c,d_isf,dtinc);
  HandleError( cudaPeekAtLastError() );
  HandleError(cudaFree(djindG));
  nThreads=512;
  nBlocks=ntot/nThreads+1;
  HandleError(cudaMalloc((void**)&ddtminG,nBlocks*sizeof(float)));
  HandleError(cudaMalloc((void**)&djindG,nBlocks*sizeof(int)));
  HandleError(cudaMalloc((void**)&dj1indG,nBlocks*sizeof(int)));
  double rX=gg.rNmax*pow(gg.dX[0]*1e6,3.);
  std::default_random_engine g1(30*gg.tInd+seed1);
  std::uniform_real_distribution<double> xrand1(0.0,1.0);
  while (isf){
    nBlocks=ntot/nThreads+1;  
    updateVoxelsPart3<<<nBlocks,nThreads,2*nThreads*sizeof(int)+nThreads*sizeof(float)>>>
      (dg,dgid,dvstate,dctheta,dtempval,troids,dexts,dtinc,dvs2c,d_isf,ddtminG,djindG,dj1indG,ntot);
    HandleError(cudaMemcpy(&isf,d_isf,sizeof(bool),cudaMemcpyDeviceToHost));
    HandleError( cudaPeekAtLastError() );
    while (nBlocks>3){
      n1=nBlocks;
      nBlocks/=2;
      reduceVoxelCapture<<<nBlocks,nThreads,2*nThreads*sizeof(int)+nThreads*sizeof(float)>>>
	(ddtminG,djindG,dj1indG,n1);
    }
    n1=nBlocks;
    reduceVoxelCapture<<<1,nThreads,2*nThreads*sizeof(int)+nThreads*sizeof(float)>>>
      (ddtminG,djindG,dj1indG,n1);
    // index 0 is winning voxel. 
    if (xrand1(g1)<rX){
      // nucleation occurs in voxel
      updateVoxelsPart4<<<1,1>>>(dg,dvox,dgid,dvstate,dctheta,troids,dexts,dvs2c,
				 ddtminG,dj1indG);
      HandleError( cudaPeekAtLastError() );
    } else {
      nBlocks=ntot/nThreads+1;
      updateVoxelsPart5<<<nBlocks,nThreads>>>(dg,dvstate,dtempval,dexts,ddtminG,ntot);
      updateVoxelsPart6<<<1,1>>>(dg,dvox,dgid,dvstate,dctheta,troids,dexts,dvs2c,dtinc,
				 ddtminG,djindG,dj1indG);
      HandleError( cudaPeekAtLastError() );
      }
  } // while (isf)
  HandleError(cudaFree(ddtminG));
  HandleError(cudaFree(djindG));
  HandleError(cudaFree(dj1indG));
}
void VoxelsCA::WriteToHDF1(const std::string &filename, const Grid &g, const double *tempcurr)
{
  // writes gID, vState, cTheta per voxel
  int Ntot = g.nX[0]*g.nX[1]*g.nX[2],icase;
  std::string hdf5Filename = filename + ".h5";
  std::vector< float> TempOut(Ntot,0),IPFmapBD(3*Ntot,0), IPFmapx(3*Ntot,0), IPFmapy(3*Ntot,0),cth(4*nGrain,0);
  double vBD[3]={0.0,0.0,1.0},omega,ax[3],vCD[3],rRot[3][3],
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
    TempOut[j] = tempcurr[j];
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
  nVoxT = g.nX[0]*g.nX[1]*g.nX[2];
  ncth=4*nGrain;
  js = 0;
  jc = nVoxT;
  std::vector<float> dX(g.dX,g.dX+3);
  adios2::ADIOS adios;
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
  hdf5Writer.Put<int>(dimsa, g.nX);
  hdf5Writer.Put<float>(dxa, dX.data());
  hdf5Writer.Put<int>(gida, gID);
  hdf5Writer.Put<int>(vStatea, vState);
  hdf5Writer.Put<float>(TempOuta, TempOut.data());
  hdf5Writer.Put<float>(IPFmapBDa, IPFmapBD.data());
  hdf5Writer.Put<float>(IPFmapxa, IPFmapx.data());
  hdf5Writer.Put<float>(IPFmapya, IPFmapy.data());
  hdf5Writer.Put<float>(angAx, cth.data());
  hdf5Writer.Close();
} // end WriteToHDF1
