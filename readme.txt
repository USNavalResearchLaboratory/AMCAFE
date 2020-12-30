Compiling and running the AMCAFE code:

The personal notes below show compilation notes for 3 systems titled neocortex, gaffney, and onyx (gaffney and onyx are DoD HPC machines), and can serveas guidelines for your system. The code is a C++ code that requires an MPI wrapper as well as external packages ADIOS2 and metis. The notes below give some guidance on how ADIOS2 can be compiled. METIS is a straightforward compilation. A static metis library can easily be compiled following instructions on the website. The variables in makefile in this folder need to specify the paths to the required libraries, executables, and include files. These must be adjusted to your system. After which, the executable can be built by typing "make cafe" in the terminal. The code is run by

mpirun -n <number_of_processors> ./cafe <input_file_name>

The executation command may need minor adjustments depending on your compiler. The *sh files in the folder give examples of how to run. The examples in paper XXXX can be run by input file SDX1_.in (example 1) and SDXY1_.in (example 2)


For questions on compiling and running, please email: kirubel.teferra@nrl.navy.mil
For usage of the code, please cite: XXXX



20200914:

creating a fork in /opt/cafe that will create branches for others to work on

20200821: compiling ADIOS2 to be able to output HDF5

adios2 needs, HDF5 parallel and cmake

neocortex:
HDF5
1) untarred hdf5 and cd'd to directory /usr/local
2) chmod 777 for the directory
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

the did another chmod -R 777 . for the entire directory

option 2

download the cmake-version-.sh online then do
sudo sh cmake-$version.$build-Linux-x86_64.sh --prefix=/opt/cmake
then you can put the path with the cmake executable in PATH or jsut give full path when use cmake. Then, not sure if necessary but I make cmake RWX for everyone

ADIOS2

I am using this cmake: /opt/cmake/cmake-3.19.0-Linux-x86_64/bin/cmake
also, i make sure i'm using this mpicc: /opt/petsc/arch-linux2-c-opt/bin/mpicc, by specifying that directory in my PATH variable in ~/.profile

I created a folder /usr/local/ADIOS2 with sudo then gave 777 permissions

1) git cloned it then created subdirectory adios2-build and cd'ed to it
cmake -DCMAKE_INSTALL_PREFIX=/usr/local/ADIOS2 -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_MPI=ON -DADIOS2_USE_HDF5=ON -DHDF5_ROOT=/usr/local/hdf5-1.12.0/ ../../ADIOS2
2) make -j 12
3) make install

Then had to again give 777 permissions in all folders

The adios-config file is in the ../ADIOS/bin directory; to determine flags to compile your application with
4) ./adios-config --cxx-flags
5) ./adios-config --cxx-flags

gaffney:

already had a module for cmake and HDF5 parallel so needed to load those modules,
1) module load Cmake/3.18.1
2) module load hdf5-parallel/intel-18.1.163/1.10.5
 then I loaded the
3) module load compiler/intelmpi/2019.5.281 
but also needed to load the gcc module to provide paths for cmake as
4) module load gcc/9.2.0
then i git cloned the ADIOS2, mkdir adios2-build and cd'ed to it
5) cmake -DCMAKE_INSTALL_PREFIX=/p/home/kteferra/Documents/software/ADIOS2/ -DADIOS2_USE_MPI=ON -DADIOS2_USE_HDF5=ON -DHDF5_ROOT=/app/hdf5-parallel/1.10.5-intel-2018.1.163-intelmpi -DADIOS2_USE_Fortran=NO ../../ADIOS
(notice that it only works if I turn off fortran)
6) Make -j 16
7) Make install
The adios-config file is in the ../ADIOS/bin directory; to determine flags to compile your application with
8) ./adios-config --cxx-flags
9) ./adios-config --cxx-flags


lastly when compile (AND EXECUTE) an application that uses adios2, you have to make sure you have the following modules loaded
module swap compiler/intel/2019.4.243 compiler/intel/2019.5.281
module swap mpt/2.20 compiler/intelmpi/2019.5.281
module load hdf5-parallel/intel-18.1.163/1.10.5
module load gcc/9.2.0

ONYX:
As onyx is a CRAZY system ADIOS2 compiled as a static library. Here are the steps:
1) module swap PrgEnv-cray PrgEnv-intel/6.0.5
2) module load intel/19.0.1.144
3) module load cray-hdf5-parallel/1.10.5.0
4) module load gcc/8.3.0

cmake:
untared cmake-3.18 then did
1) ./bootstrap
2) make -j 12
3) make install

then git cloned adios2, created subdir adios2-build and cc'ed to it
1) ../../cmake-3.18.2/bin/cmake -DCMAKE_INSTALL_PREFIX=/p/home/kteferra/Documents/software/ADIOS2 -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_MPI=ON -DADIOS2_USE_HDF5=ON ../../ADIOS2
2) make -j 12
3) make install

you can find compiler flags by running executable adios2-config in ADIOS2/bin

20200526:
if system does not have MPI wrapper, you can download petsc tar file, untar and install with:
./configure --with-cc=gcc --with-cxx=g++ --with-debugging=0 COPTFLAGS='-O3 -march=native -mtune=native' 
CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native'
 --download-mpich --download-metis

That being said, on HPC it probably is better to use an intel compiler rather than GNU and to
use any MPI library that is compiled already on the system by system administrators. That is likely
to be more optimized.


