1. Description
===============

This is an implementation of the Cellular Automata Finite Element (CAFE) algorithm to simulate solidification of additively manufactured metals. Description of the code can be found in

Teferra, Kirubel, and David J. Rowenhorst. "Optimizing the cellular automata finite element model for additive manufacturing to simulate large microstructures." Acta Materialia 213 (2021): 116930

Simulation results in the paper can reproduced using the input files in the examples directory. Below contains compilation instructions for various systems. These instructions mostly serve as guidelines to help compile on your own system. In addition to compilation instructions, the src/ directoriescontain the makefiles used to compile the code.


2. Compiling and running the AMCAFE code
==========================================

The personal notes below show compilation notes for 3 systems titled neocortex, gaffney, and onyx (gaffney and onyx are DoD HPC machines, neocortex is ubuntu OS). The code is a C++ code that requires an MPI wrapper as well as external packages ADIOS2 and metis. The notes below give some guidance on how ADIOS2 can be compiled. METIS is a straightforward compilation. A static metis library can easily be compiled following instructions on the website. The variables in makefile in this folder need to specify the paths to the required libraries, executables, and include files. These must be adjusted to your system. After which, the executable can be built by typing "make cafe" in the terminal. The code is run by

mpirun -n <number_of_processors> ./cafe <input_file_name>

The executation command may need minor adjustments depending on your compiler. The *sh files in the folder give examples of how to run. The examples in paper mentioned above can be run by input file SDX1_.in (example 1) and SDXY1_.in (example 2)


For questions on compiling and running, please email: kirubel.teferra@nrl.navy.mil




2.1 compiling ADIOS2 to be able to output HDF5
===============================================

dependencies: HDF5 parallel and cmake

neocortex:
HDF5
1) untarred hdf5 and cd'd to directory /usr/local
2) chmod 755 for the directory
3) did the following command
CC=mpicc ./configure --enable-parallel
make -j 12
make install

* note that mpicc refers to a petsc build already on neocortex:
/opt/petsc/arch-linux2-c-opt/bin/mpicc

CMAKE

option 1

1) downloaded cmake source and untarred, cd'd to directory
2) needed the openssl library (had an error first)
sudo apt get install libssl-dev
3) then built cmake as
./bootstrap
sudo make -j 12
make install

the did another chmod -R 755 . for the entire directory

option 2

download the cmake-version-.sh online then do
sudo sh cmake-$version.$build-Linux-x86_64.sh --prefix=/opt/cmake
then you can put the path with the cmake executable in PATH or just give full path when use cmake. Then, not sure if necessary but I make cmake RWX for everyone

ADIOS2

I am using this cmake: /opt/cmake/cmake-3.19.0-Linux-x86_64/bin/cmake
also, i make sure i'm using this mpicc: /opt/petsc/arch-linux2-c-opt/bin/mpicc, by specifying that directory in my PATH variable in ~/.profile

I created a folder /usr/local/ADIOS2 with sudo then gave 777 permissions

1) git cloned it then created subdirectory adios2-build and cd'ed to it
cmake -DCMAKE_INSTALL_PREFIX=/usr/local/ADIOS2 -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_MPI=ON -DADIOS2_USE_HDF5=ON -DHDF5_ROOT=/usr/local/hdf5-1.12.0/ ../../ADIOS2
2) make -j 12
3) make install

Then had to again give 777 permissions in all created folders

The adios-config file is in the ../ADIOS/bin directory; to determine flags to compile your application with
4) ./adios-config --cxx-flags
5) ./adios-config --cxx-libs

gaffney:
=========

already had a module for cmake and HDF5 parallel so needed to load those modules,
1) module load Cmake/3.18.1
2) module load hdf5-parallel/intel-18.1.163/1.10.5
3) module load compiler/intelmpi/2019.5.281 
4) module load gcc/9.2.0 (necessary to provide paths for cmake)
5) compile ADIOS2
5a) make directory for ADIOS2, git cloned to directory, mkdir adios2-build and cd'ed to it
5b) cmake -DCMAKE_INSTALL_PREFIX=/p/home/kteferra/Documents/software/ADIOS2/ -DADIOS2_USE_MPI=ON -DADIOS2_USE_HDF5=ON -DHDF5_ROOT=/app/hdf5-parallel/1.10.5-intel-2018.1.163-intelmpi -DADIOS2_USE_Fortran=NO ../../ADIOS
(notice that it only works if I turn off fortran)
5c) Make -j 16
5d) Make install
The adios-config file is in the ../ADIOS/bin directory; to determine flags to compile your application with
5e) ./adios-config --cxx-flags
5f) ./adios-config --cxx-libs


lastly when compile (AND EXECUTE) an application that uses adios2, you have to make sure you have the following modules loaded
module swap compiler/intel/2019.4.243 compiler/intel/2019.5.281
module swap mpt/2.20 compiler/intelmpi/2019.5.281
module load hdf5-parallel/intel-18.1.163/1.10.5
module load gcc/9.2.0

ONYX:
===========

As onyx is a CRAY system ADIOS2 compiled as a static library. Here are the steps:
1) load appropriate system modules
1a) module swap PrgEnv-cray PrgEnv-intel/6.0.5
1b) module load intel/19.0.1.144
1c) module load cray-hdf5-parallel/1.10.5.0
1d) module load gcc/8.3.0

2)cmake:
2a) downloaded and untared cmake-3.18 then did
2b) ./bootstrap
2c) make -j 12
2d) make install

3) ADIOS2
3a) created ADIOS2 directroy, git cloned ADIOS2 to directory, created subdir adios2-build and cc'ed to it
3b) ../../cmake-3.18.2/bin/cmake -DCMAKE_INSTALL_PREFIX=/p/home/kteferra/Documents/software/ADIOS2 -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_MPI=ON -DADIOS2_USE_HDF5=ON ../../ADIOS2
3c) make -j 12
3d) make install

you can find compiler flags by running executable adios2-config in ADIOS2/bin

20200526:
if system does not have MPI wrapper, you can download petsc tar file, untar and install with:
./configure --with-cc=gcc --with-cxx=g++ --with-debugging=0 COPTFLAGS='-O3 -march=native -mtune=native' 
CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native'
 --download-mpich --download-metis


On DoD HPC system: better to use intel compiler than gcc and use as many existing modules as possible as they are optimized for system


