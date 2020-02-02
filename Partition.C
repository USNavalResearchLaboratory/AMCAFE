#include "Grid.h"
#include "Partition.h"
#include "iostream"
#include "metis.h"
#include <numeric>
#include <algorithm>
#include <vector>
#include <math.h>
#include "mpi.h"
#include "VoxelsCA.h"

#include "fstream"

//constructor
Partition::Partition(Grid &g,int & myidIn,int & nprocsIn)
{
  _xyz = &g;
  myid = myidIn;
  nprocs=nprocsIn;
} // end constructor

void Partition::PartitionGraph()
{

  if (myid==0){
    std::vector<int> cellNeigh,ghostId,cellGhost;
    icellidLocAll.assign(nprocs,std::vector<int>(0));
    idx_t myoptions[METIS_NOPTIONS];
    myoptions[METIS_OPTION_NUMBERING] = 0; // 0 for 0 based, 1 for 1 based #ing system
    std::vector<int> n=_xyz->nX;
    int nelem,nDim,j,cc,cc1,cc2,nneigh,nelemLoc,jn,j1,j2,j3,nedge;
    std::vector<int> :: iterator tmp;
    std::vector<int> neigh,y1;
    nelem = std::accumulate(n.begin(),n.end(),1,std::multiplies<int>());
    std::vector<int> myeptr(nelem+1);
    cc=0;
    for (int j=0;j<nelem;++j){
      _xyz->ComputeNeighborhood(j,_xyz->neighOrder,neigh);
      myeptr[j] = cc;
      cc+=neigh.size();
    }
    myeptr[nelem]=cc;
    nedge = cc;
    //idx_t myeind[nedge];
    std::vector<int> myeind(nedge);
    cc=0;
    for (int j=0;j<nelem;++j){
      _xyz->ComputeNeighborhood(j,_xyz->neighOrder,neigh);
      for (int jn=0;jn<neigh.size();++jn){
	myeind[cc] = neigh[jn];
	cc+=1;
      }
    }
    idx_t ne=nelem,ncon=1,nn = 1,objval;
    std::vector<int> epart(nelem);
    for (int j=0;j<_xyz->nDim;++j){nn*=(n[j]+1);}
    // if (nprocs < 8){
    if (nprocs > 0){
      METIS_PartGraphRecursive(&ne,&ncon,&myeptr[0],&myeind[0],NULL,NULL,NULL,&nprocs,NULL,NULL,NULL,&objval,&epart[0]);
    } else{
      METIS_PartGraphKway(&ne,&ncon,&myeptr[0],&myeind[0],NULL,NULL,NULL,&nprocs,NULL,NULL,NULL,&objval,&epart[0]);
    }
    for (int jid=nprocs-1;jid > -1;--jid){
      cc=0;
      cc1=0;
      icellidLoc.assign(nelem,0);
      cellNeigh.assign(nedge,0);
      for (int j=0;j<nelem;++j){
        if (epart[j]==jid){
          icellidLoc[cc] = j;
          cc+=1;
          for (int j1=myeptr[j];j1<myeptr[j+1];++j1){
            cellNeigh[cc1] = myeind[j1];
            cc1+=1;
          }
        }
      }
      ncellLoc = cc;
      cellNeigh.resize(cc1);
      // make cellNeigh unique and check which proc owns it
      std::sort(cellNeigh.begin(),cellNeigh.end());
      tmp = std::unique(cellNeigh.begin(),cellNeigh.end());
      cellNeigh.resize(std::distance(cellNeigh.begin(),tmp));
      nneigh = cellNeigh.size();
      cellGhost.assign(nneigh,0);
      ghostId.assign(nneigh,0);
      cc=0;
      for (int j=0;j<nneigh;++j){
        if (epart[cellNeigh[j]] !=jid){
          cellGhost[cc] = cellNeigh[j];
          ghostId[cc] = epart[cellNeigh[j]];
          cc+=1;
        }
      }
      ghostId.resize(cc);
      cellGhost.resize(cc);
      nGhost = cc;
      nelemLoc = nGhost + ncellLoc;
      icellidLoc.resize(nelemLoc);
      y1.assign(nGhost,0);
      std::iota(y1.begin(),y1.end(),0);
      std::sort(y1.begin(),y1.end(),[&cellGhost](int i1,int i2){return cellGhost[i1]<cellGhost[i2];});
      for (int j=0;j<nGhost;++j){
        icellidLoc[j+ncellLoc]=cellGhost[y1[j]];
      };
      ineighProcId = ghostId;
      std::sort(ineighProcId.begin(),ineighProcId.end());
      tmp = std::unique(ineighProcId.begin(),ineighProcId.end());
      ineighProcId.resize(std::distance(ineighProcId.begin(),tmp));
      nneighProc = ineighProcId.size();
      ineigh2Locptr.assign(nneighProc+1,0);
      ineigh2LocVals.assign(nGhost,0);
      ineigh2Locptr[0]=0;
      cc1=0;
      cc=0;
      for (int j=0;j<nneighProc;++j){
	for (int j1 =0;j1<nGhost;++j1){
	  if (ghostId[y1[j1]]==ineighProcId[j]){
	    ineigh2LocVals[cc1] =  std::distance(icellidLoc.begin(),std::find(icellidLoc.begin(),icellidLoc.end(),cellGhost[y1[j1]]));
	    cc1+=1;
	  }
	}
	cc+=1;
	ineigh2Locptr[cc] = cc1;
      }
      iloc2Neighptr.assign(nneighProc+1,0);
      iloc2NeighVals.assign(nedge,0);
      cc2=0;
      for (int jn=0;jn<nneighProc;++jn){
        cellNeigh.assign(nedge,0);
        cc1=0;
        for (int j=0;j<nelem;++j){
          if (epart[j]==ineighProcId[jn]){
            for (int j1=myeptr[j];j1<myeptr[j+1];++j1){
              cellNeigh[cc1] = myeind[j1];
              cc1+=1;
            }
          }
        }
        cellNeigh.resize(cc1);
	std::sort(cellNeigh.begin(),cellNeigh.end());
        tmp = std::unique(cellNeigh.begin(),cellNeigh.end());
        cellNeigh.resize(std::distance(cellNeigh.begin(),tmp));
        nneigh = cellNeigh.size();
        iloc2Neighptr[jn]=cc2;
        for (int j=0;j<nneigh;++j){
          if (epart[cellNeigh[j]] ==jid){
            iloc2NeighVals[cc2] =std::distance(icellidLoc.begin(),std::find(icellidLoc.begin(),icellidLoc.end(),cellNeigh[j]));
            cc2+=1;
          }
        }
      } // for jn
      iloc2Neighptr[nneighProc]=cc2;
      iloc2NeighVals.resize(cc2);
      // compute pointIDs for cells
      cc2 = _xyz->nnodePerCell;
      ipointidLoc.assign(cc2*ncellLoc,0);
      if (_xyz->nDim ==2){
	for (int j =0;j<ncellLoc;++j){
	  j2 = floor(icellidLoc[j]/n[0]); 
	  j1 = icellidLoc[j] - n[0]*j2;
	  jn = (n[0]+1)*j2 + j1;
	  ipointidLoc[cc2*j] = jn;
	  ipointidLoc[cc2*j+1] = jn+1;
	  ipointidLoc[cc2*j+2] = jn+(n[0]+1)+1;
	  ipointidLoc[cc2*j+3] = jn+(n[0]+1);
	}
      } else {
	for (int j =0;j<ncellLoc;++j){
	  j3 = floor(icellidLoc[j]/(n[0]*n[1]));
	  j2 = floor( (icellidLoc[j] - n[0]*n[1]*j3)/n[0]);
	  j1 = icellidLoc[j] - n[0]*n[1]*j3 - n[0]*j2;
	  jn = (n[0]+1)*(n[1]+1)*j3 + (n[0]+1)*j2 + j1;
	  ipointidLoc[cc2*j] = jn;
	  ipointidLoc[cc2*j+1] = jn+1;
	  ipointidLoc[cc2*j+2] = jn+(n[0]+1)+1;
	  ipointidLoc[cc2*j+3] = jn+(n[0]+1);
	  ipointidLoc[cc2*j+4] = jn + (n[0]+1)*(n[1]+1);
	  ipointidLoc[cc2*j+5] = jn+1  + (n[0]+1)*(n[1]+1);
	  ipointidLoc[cc2*j+6] = jn+(n[0]+1)+1 + (n[0]+1)*(n[1]+1);
	  ipointidLoc[cc2*j+7] = jn+(n[0]+1) + (n[0]+1)*(n[1]+1);
	}
      } // if (_xyz->nDim ==2
      std::sort(ipointidLoc.begin(),ipointidLoc.end());
      tmp =std::unique(ipointidLoc.begin(),ipointidLoc.end());
      ipointidLoc.resize(std::distance(ipointidLoc.begin(),tmp));
      npointLoc = ipointidLoc.size();
      iconnectivityLoc.assign(_xyz->nnodePerCell*ncellLoc,0);

      if (_xyz->nDim==2){
	for (int j=0;j<ncellLoc;++j){
	  j2 = floor(icellidLoc[j]/n[0]); 
	  j1 = icellidLoc[j] - n[0]*j2;
	  jn = (n[0]+1)*j2 + j1;
	  iconnectivityLoc[cc2*j]=std::distance(ipointidLoc.begin(),std::find(ipointidLoc.begin(),ipointidLoc.end(),jn));
	  iconnectivityLoc[cc2*j+1]=std::distance(ipointidLoc.begin(),std::find(ipointidLoc.begin(),ipointidLoc.end(),jn+1));
	  iconnectivityLoc[cc2*j+2]=std::distance(ipointidLoc.begin(),std::find(ipointidLoc.begin(),ipointidLoc.end(),jn+1+n[0]+1));
	  iconnectivityLoc[cc2*j+3]=std::distance(ipointidLoc.begin(),std::find(ipointidLoc.begin(),ipointidLoc.end(),jn+n[0]+1));
	} // for (int j...
      } else {
	for (int j=0;j<ncellLoc;++j){
	  j3 = floor(icellidLoc[j]/(n[0]*n[1]));
	  j2 = floor( (icellidLoc[j] - n[0]*n[1]*j3)/n[0]);
	  j1 = icellidLoc[j] - n[0]*n[1]*j3 - n[0]*j2;
	  jn = (n[0]+1)*(n[1]+1)*j3 + (n[0]+1)*j2 + j1;
	  iconnectivityLoc[cc2*j]=std::distance(ipointidLoc.begin(),std::find(ipointidLoc.begin(),ipointidLoc.end(),jn));
	  iconnectivityLoc[cc2*j+1]=std::distance(ipointidLoc.begin(),std::find(ipointidLoc.begin(),ipointidLoc.end(),jn+1));
	  iconnectivityLoc[cc2*j+2]=std::distance(ipointidLoc.begin(),std::find(ipointidLoc.begin(),ipointidLoc.end(),jn+1+n[0]+1));
	  iconnectivityLoc[cc2*j+3]=std::distance(ipointidLoc.begin(),std::find(ipointidLoc.begin(),ipointidLoc.end(),jn+n[0]+1));
	  iconnectivityLoc[cc2*j+4]=std::distance(ipointidLoc.begin(),std::find(ipointidLoc.begin(),ipointidLoc.end(),jn+(n[0]+1)*(n[1]+1)));
	  iconnectivityLoc[cc2*j+5]=std::distance(ipointidLoc.begin(),std::find(ipointidLoc.begin(),ipointidLoc.end(),jn+1+(n[0]+1)*(n[1]+1)));
	  iconnectivityLoc[cc2*j+6]=std::distance(ipointidLoc.begin(),std::find(ipointidLoc.begin(),ipointidLoc.end(),jn+1+n[0]+1+(n[0]+1)*(n[1]+1)));
	  iconnectivityLoc[cc2*j+7]=std::distance(ipointidLoc.begin(),std::find(ipointidLoc.begin(),ipointidLoc.end(),jn+n[0]+1+(n[0]+1)*(n[1]+1)));
	}
      } // if (_xyz->nDim==2...
      icellidLocAll[jid] = icellidLoc;
      // send ineigh2Locptr,ineigh2LocVals, iloc2Neighptr, iloc2NeighVals, icellidLoc,ineighProcId,nneighProc,nGhost,ncellLoc
      if (jid >0){
        // SEND Data
	MPI_Send(&nGhost,1,MPI_INT,jid,0,MPI_COMM_WORLD);
        MPI_Send(&ncellLoc,1,MPI_INT,jid,0,MPI_COMM_WORLD);
        MPI_Send(&nneighProc,1,MPI_INT,jid,0,MPI_COMM_WORLD);
        MPI_Send(&npointLoc,1,MPI_INT,jid,0,MPI_COMM_WORLD);
        MPI_Send(&icellidLoc[0],(nGhost+ncellLoc),MPI_INT,jid,0,MPI_COMM_WORLD);
        MPI_Send(&ipointidLoc[0],npointLoc,MPI_INT,jid,0,MPI_COMM_WORLD);
        MPI_Send(&ineighProcId[0],nneighProc,MPI_INT,jid,0,MPI_COMM_WORLD);
        MPI_Send(&iloc2Neighptr[0],(nneighProc+1),MPI_INT,jid,0,MPI_COMM_WORLD);
        MPI_Send(&ineigh2Locptr[0],(nneighProc+1),MPI_INT,jid,0,MPI_COMM_WORLD);
        MPI_Send(&iloc2NeighVals[0],iloc2Neighptr[nneighProc],MPI_INT,jid,0,MPI_COMM_WORLD);
        MPI_Send(&ineigh2LocVals[0],ineigh2Locptr[nneighProc],MPI_INT,jid,0,MPI_COMM_WORLD);
        MPI_Send(&iconnectivityLoc[0],_xyz->nnodePerCell*ncellLoc,MPI_INT,jid,0,MPI_COMM_WORLD);
      } // if (jid >0)
    } // for jid
  } // if (myid==0)
  if (myid > 0){
    // receive data
    MPI_Recv(&nGhost,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(&ncellLoc,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(&nneighProc,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(&npointLoc,1,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    icellidLoc.assign((nGhost+ncellLoc),0);
    ineighProcId.assign(nneighProc,0);
    iloc2Neighptr.assign(nneighProc+1,0);
    ineigh2Locptr.assign(nneighProc+1,0);
    ipointidLoc.assign(npointLoc,0);
    MPI_Recv(&icellidLoc[0],(nGhost+ncellLoc),MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(&ipointidLoc[0],npointLoc,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(&ineighProcId[0],nneighProc,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(&iloc2Neighptr[0],(nneighProc+1),MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(&ineigh2Locptr[0],(nneighProc+1),MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    iloc2NeighVals.assign(iloc2Neighptr[nneighProc],0);
    ineigh2LocVals.assign(ineigh2Locptr[nneighProc],0);
    iconnectivityLoc.assign(_xyz->nnodePerCell*ncellLoc,0);
    MPI_Recv(&iloc2NeighVals[0],iloc2Neighptr[nneighProc],MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(&ineigh2LocVals[0],ineigh2Locptr[nneighProc],MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    MPI_Recv(&iconnectivityLoc[0],_xyz->nnodePerCell*ncellLoc,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
  } // if (myid > 0)
} // end PartitionGraph()

void Partition::PartitionGraph2(){
  /*
    - Creates the partition of the domain by evenly dividing the domain among all the processors,
    Each domain is a contiguous region where indices go x,y,z. Thus, imagine a block where smallest
    dimension is the z direction. 
    - This approach does not require METIS since creating the partition "manually"- this enables not having
    to do all on 1 processor and MPI_sending to all others

   */
  int nelemT = _xyz->nX[0]*_xyz->nX[1]*_xyz->nX[2],nelemLoc;
  int i1,i2;
  i1 = ceil( (double)nelemT / (double)nprocs);
  i2 = floor( (double)nelemT / (double)i1);
  std::vector<int> cellNeigh,ghostId,cellGhost,neigh, 
    n = _xyz->nX, iv(nprocs+1,0);
  for (int j=0;j<(nprocs+1);++j){
    if (j < (i2+1)){iv[j] = i1*j;}
    if (j>i2 && j<nprocs){iv[j] = i2*i1 +
        floor( (double)(nelemT-i2*i1)/(double)(nprocs-i2));}
    if (j==nprocs){iv[j]=nelemT;}
  } // for (int j...
  ncellLoc = iv[myid+1] - iv[myid];
  icellidLoc.assign(ncellLoc,0);
  cellNeigh.assign(26*ncellLoc,0);
  int j1,j2,j3,jn,cc=0,cc1=0,nneigh,cc2;
  std::vector<int> :: iterator tmp;
  for (int j=iv[myid];j<iv[myid+1];++j){
    icellidLoc[cc] = j;
    cc+=1;
    _xyz->ComputeNeighborhoodMooreFirst(j,neigh);
    i1=neigh.size();
    for (int j1=0;j1<i1;++j1){cellNeigh[cc1+j1]=neigh[j1];}
    cc1+=i1;
  } // for (int j ...
  cellNeigh.resize(cc1);
  // make cellNeigh unique and check which proc owns it  
  std::sort(cellNeigh.begin(),cellNeigh.end());
  tmp = std::unique(cellNeigh.begin(),cellNeigh.end());
  cellNeigh.resize(std::distance(cellNeigh.begin(),tmp));
  nneigh = cellNeigh.size();
  cellGhost.assign(nneigh,0);
  ghostId.assign(nneigh,0);
  cc=0;
  for (int j=0;j<nneigh;++j){
    i1 = std::distance(iv.begin(),std::upper_bound(iv.begin(),
	   iv.end(),cellNeigh[j])) - 1; // i1= proc owns cellNeigh[j]
    if (i1 !=myid){
      cellGhost[cc] = cellNeigh[j];
      ghostId[cc] = i1;
      cc+=1;
    } // if (i1 !=...
  } // for (int j ...
  ghostId.resize(cc);
  cellGhost.resize(cc);
  nGhost = cc;
  nelemLoc = nGhost + ncellLoc;
  icellidLoc.resize(nelemLoc);
  for (int j=0;j<nGhost;++j){
     icellidLoc[j+ncellLoc]=cellGhost[j];
  }
  ineighProcId = ghostId;
  std::sort(ineighProcId.begin(),ineighProcId.end());
  tmp = std::unique(ineighProcId.begin(),ineighProcId.end());
  ineighProcId.resize(std::distance(ineighProcId.begin(),tmp));
  nneighProc = ineighProcId.size();
  ineigh2Locptr.assign(nneighProc+1,0);
  ineigh2LocVals.assign(nGhost,0);
  ineigh2Locptr[0]=0;
  cc1=0;
  for (int j=0;j<nneighProc;++j){
    for (int j1 =0;j1<nGhost;++j1){
      if (ghostId[j1]==ineighProcId[j]){
	ineigh2LocVals[cc1] =  std::distance(icellidLoc.begin(),
		   std::find(icellidLoc.begin(),icellidLoc.end(),cellGhost[j1]));
	cc1+=1;
      }
    }
    ineigh2Locptr[j+1] = cc1;
  }
  iloc2Neighptr.assign(nneighProc+1,0);
  iloc2NeighVals.assign(ncellLoc*nneighProc,0);
  cc2=0;
  for (int jn=0;jn<nneighProc;++jn){
    i1 = iv[ineighProcId[jn]+1] - iv[ineighProcId[jn]];
    cellNeigh.assign(26*i1,0);
    cc1=0;
    for (int j=iv[ineighProcId[jn]];j<iv[ineighProcId[jn]+1];++j){
      _xyz->ComputeNeighborhoodMooreFirst(j,neigh);
      i1=neigh.size();
      for (int j1=0;j1<i1;++j1){cellNeigh[cc1+j1]=neigh[j1];}
      cc1+=i1;
    } // for (int j ..
    cellNeigh.resize(cc1);
    std::sort(cellNeigh.begin(),cellNeigh.end());
    tmp = std::unique(cellNeigh.begin(),cellNeigh.end());
    cellNeigh.resize(std::distance(cellNeigh.begin(),tmp));
    nneigh = cellNeigh.size();
    for (int j=0;j<nneigh;++j){
      i1 = std::distance(iv.begin(),std::upper_bound(iv.begin(),
	  iv.end(),cellNeigh[j])) - 1; // i1= proc owns cellNeigh[j]
      if (i1 ==myid){
	iloc2NeighVals[cc2] =std::distance(icellidLoc.begin(),
		 std::find(icellidLoc.begin(),icellidLoc.begin()+ncellLoc,cellNeigh[j]));
	cc2+=1;
      }      
    }
    iloc2Neighptr[jn+1]=cc2;
  } // for (int jn ...
  iloc2NeighVals.resize(cc2);
  cc2 = _xyz->nnodePerCell;
  ipointidLoc.assign(cc2*ncellLoc,0);
  for (int j =0;j<ncellLoc;++j){
    j3 = floor(icellidLoc[j]/(n[0]*n[1]));
    j2 = floor( (icellidLoc[j] - n[0]*n[1]*j3)/n[0]);
    j1 = icellidLoc[j] - n[0]*n[1]*j3 - n[0]*j2;
    jn = (n[0]+1)*(n[1]+1)*j3 + (n[0]+1)*j2 + j1;
    ipointidLoc[cc2*j] = jn;
    ipointidLoc[cc2*j+1] = jn+1;
    ipointidLoc[cc2*j+2] = jn+(n[0]+1)+1;
    ipointidLoc[cc2*j+3] = jn+(n[0]+1);
    ipointidLoc[cc2*j+4] = jn + (n[0]+1)*(n[1]+1);
    ipointidLoc[cc2*j+5] = jn+1  + (n[0]+1)*(n[1]+1);
    ipointidLoc[cc2*j+6] = jn+(n[0]+1)+1 + (n[0]+1)*(n[1]+1);
    ipointidLoc[cc2*j+7] = jn+(n[0]+1) + (n[0]+1)*(n[1]+1);
  }
  std::sort(ipointidLoc.begin(),ipointidLoc.end());
  tmp =std::unique(ipointidLoc.begin(),ipointidLoc.end());
  ipointidLoc.resize(std::distance(ipointidLoc.begin(),tmp));
  npointLoc = ipointidLoc.size();
  iconnectivityLoc.assign(cc2*ncellLoc,0);
  for (int j=0;j<ncellLoc;++j){
    j3 = floor(icellidLoc[j]/(n[0]*n[1]));
    j2 = floor( (icellidLoc[j] - n[0]*n[1]*j3)/n[0]);
    j1 = icellidLoc[j] - n[0]*n[1]*j3 - n[0]*j2;
    jn = (n[0]+1)*(n[1]+1)*j3 + (n[0]+1)*j2 + j1;
    iconnectivityLoc[cc2*j]=std::distance(ipointidLoc.begin(),
      std::find(ipointidLoc.begin(),ipointidLoc.end(),jn));
    iconnectivityLoc[cc2*j+1]=std::distance(ipointidLoc.begin(),
      std::find(ipointidLoc.begin(),ipointidLoc.end(),jn+1));
    iconnectivityLoc[cc2*j+2]=std::distance(ipointidLoc.begin(),
      std::find(ipointidLoc.begin(),ipointidLoc.end(),jn+1+n[0]+1));
    iconnectivityLoc[cc2*j+3]=std::distance(ipointidLoc.begin(),
      std::find(ipointidLoc.begin(),ipointidLoc.end(),jn+n[0]+1));
    iconnectivityLoc[cc2*j+4]=std::distance(ipointidLoc.begin(),
      std::find(ipointidLoc.begin(),ipointidLoc.end(),jn+(n[0]+1)*(n[1]+1)));
    iconnectivityLoc[cc2*j+5]=std::distance(ipointidLoc.begin(),
      std::find(ipointidLoc.begin(),ipointidLoc.end(),jn+1+(n[0]+1)*(n[1]+1)));
    iconnectivityLoc[cc2*j+6]=std::distance(ipointidLoc.begin(),
      std::find(ipointidLoc.begin(),ipointidLoc.end(),jn+1+n[0]+1+(n[0]+1)*(n[1]+1)));
    iconnectivityLoc[cc2*j+7]=std::distance(ipointidLoc.begin(),
      std::find(ipointidLoc.begin(),ipointidLoc.end(),jn+n[0]+1+(n[0]+1)*(n[1]+1)));
  } // for (int j...
} // end PartitionGraph2()

void Partition::PassInformation(std::vector<int> &ivec)
{
  // passes and receives ghost values of ivec among all processors
  std::vector<int> ivecsend,ivecrecv;
  int nsend,nrecv;
  MPI_Request reqs[2*nneighProc];
  MPI_Status stats[2*nneighProc];
  for (int j2=0;j2<nneighProc;++j2){
    nsend = iloc2Neighptr[j2+1]-iloc2Neighptr[j2];
    ivecsend.assign(nsend,0);
    for (int j=0;j<nsend;++j){ivecsend[j]=ivec[iloc2NeighVals[iloc2Neighptr[j2]+j]];}
    nrecv = ineigh2Locptr[j2+1]-ineigh2Locptr[j2];
    ivecrecv.assign(nrecv,0);
    MPI_Sendrecv(&ivecsend[0],nsend,MPI_INT,ineighProcId[j2],0,
     &ivecrecv[0],nrecv,MPI_INT,ineighProcId[j2],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    for (int j=0;j<nrecv;++j){ivec[ineigh2LocVals[ineigh2Locptr[j2]+j]]=ivecrecv[j];}
  } // for j2
} // end PassInformation(std::vector<int> &ivec)
void Partition::PassInformation(std::vector<double> &vec)
{
  // passes and receives ghost values of ivec among all processors
  std::vector<double> vecsend,vecrecv;
  int nsend,nrecv;
  MPI_Request reqs[2*nneighProc];
  MPI_Status stats[2*nneighProc];
  for (int j2=0;j2<nneighProc;++j2){
    nsend = iloc2Neighptr[j2+1]-iloc2Neighptr[j2];
    vecsend.assign(nsend,0);
    for (int j=0;j<nsend;++j){vecsend[j]=vec[iloc2NeighVals[iloc2Neighptr[j2]+j]];}
    nrecv = ineigh2Locptr[j2+1]-ineigh2Locptr[j2];
    vecrecv.assign(nrecv,0);
    MPI_Sendrecv(&vecsend[0],nsend,MPI_DOUBLE,ineighProcId[j2],0,
		 &vecrecv[0],nrecv,MPI_DOUBLE,ineighProcId[j2],MPI_ANY_TAG,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    for (int j=0;j<nrecv;++j){vec[ineigh2LocVals[ineigh2Locptr[j2]+j]]=vecrecv[j];}
  } // for j2
} // end PassInformation(std::vector<double> &vec)
