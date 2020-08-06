// member function definitions for Grid.C

#include "Param.h"
#include "mpi.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>

Param::Param(std::string &filInName,int & myid,int & numprocs)
{

readFile(filInName,myid,numprocs);

} // end constructor
void Param::readFile(std::string &filInName, int &myid,int & numprocs)
{

  std::ifstream filIn;
  std::string inputData,keyword;
  int k=0,n=0;
  if (myid==0){
    filIn.open(filInName.c_str());
    while (!filIn.eof()){
      char c;
      filIn.get(c);
      if (c == '#') {
	// skip comment                                                                                                                                              
	while (filIn && c != '\n') {
	  filIn.get(c);
	}
      }

      if (k >= n) {
	n = 2*n;
	const int m = 10000;
	if (n < m) n = m;
	inputData.resize(n);
      }

      inputData[k++] = c;
    }
    filIn.close();
    inputData.resize(k);
  }
  
  MPI_Bcast(&k,1,MPI_INT,0,MPI_COMM_WORLD);
  inputData.resize(k);
  MPI_Bcast(&inputData[0],k,MPI_CHAR,0,MPI_COMM_WORLD);
  std::istringstream simInput(inputData);
  simInput >> keyword;
  while(simInput){
    std::transform(keyword.begin(),keyword.end(),keyword.begin(),
                   [](unsigned char s){return std::tolower(s);});
    if (keyword=="dx") {
      dX.resize(3);
      simInput >> keyword;
      dX[0]=std::stod(keyword);
      dX[1]=dX[0];
      dX[2]=dX[0];
    }
    if (keyword=="bmv") {
      simInput >> keyword;
      bmV=std::stod(keyword);
    }
    if (keyword=="bmp") {
      simInput >> keyword;
      bmP=std::stod(keyword);
    }
    if (keyword=="nx") {
      nX.resize(3);
      simInput >> keyword;
      nX[0]=std::stoi(keyword);
      simInput >> keyword;
      nX[1]=std::stoi(keyword);
      simInput >> keyword;
      nX[2]=std::stoi(keyword);
    }
    if (keyword=="r") {
      simInput >> keyword;
      rnuc=std::stod(keyword);
    }
    if (keyword=="lt") {
      simInput >> keyword;
      layerT=std::stod(keyword);
    }
    if (keyword=="tl") {
      simInput >> keyword;
      tL=std::stod(keyword);
    }
    if (keyword=="ts") {
      simInput >> keyword;
      tS=std::stod(keyword);
    }
    if (keyword=="bhatch") {
      simInput >> keyword;
      bhatch=std::stod(keyword);
    }
    if (keyword=="meltparam") {
      meltparam.resize(4);
      simInput >> keyword;
      meltparam[0]=std::stod(keyword);
      simInput >> keyword;
      meltparam[1]=std::stod(keyword);
      simInput >> keyword;
      meltparam[2]=std::stod(keyword);
      simInput >> keyword;
      meltparam[3]=std::stod(keyword);
    }
    if (keyword=="patternid") {
      simInput >> keyword;
      patternID=std::stoi(keyword);
    }
    if (keyword=="outint") {
      simInput >> keyword;
      outint=std::stoi(keyword);
    }

    simInput >> keyword;
  } // while(simInput)
} // end readFile
