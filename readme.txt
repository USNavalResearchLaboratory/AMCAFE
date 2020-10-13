notes for compilation and executing the NRL cafe code AMCAFE.

In order to compile and execute:
1) need metis
   - go to metis library and install static library of metis
   - http://glaros.dtc.umn.edu/gkhome/metis/metis/download
2) need adios2
   - there are various steps depending on existing environment
   - adios2 needs cmake, HDF5-parallel, gcc
   - follow instructions: https://adios2.readthedocs.io/en/latest/
3) customize makefile 
   - change variables to appropriate paths
4) type "make cafe" in terminal
5) run cafe executable: specific command depends on compiler used


below are some notes on DoD HPC systems:

neocortex:
compile HDF5
1) untarred hdf5 and cd'd to directory
2) did the following command
CC=mpicc ./configure --enable-parallel
make -j 12
make install

* note that mpicc refers to a mpich build already on neocortex.
  another build should have been ok.

compile CMAKE
1) downloaded cmake source and untarred, cd'd to directory
2) needed the openssl library (had an error first)
sudo apt get install libssl-dev
3) then built cmake as
./bootstrap
sudo make -j 12
make install

compile ADIOS2
1) git cloned it then created subdirectory adios2-build and cd'ed to it
cmake -DCMAKE_INSTALL_PREFIX=/home/kteferra/Documents/research/software/ADIOS2 -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_MPI=ON -DADIOS2_USE_HDF5=ON -DHDF5_ROOT=/home/kteferra/Documents/research/software/hdf5-1.12.0/ ../../ADIOS2
2)sudo make -j 12
3) make install

The adios-config file is in the ../ADIOS/bin directory; to determine flags to compile your application with
8) ./adios-config --cxx-flags
9) ./adios-config --cxx-flags

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
9) ./adios-config --cxx-libs

However, for some reason, the cxx-libs does not give all the necessary libraries. In order to get it to work i needed to use this command in the makefile
for the adioslibrary
adioslib = -Wl,-rpath,/p/home/kteferra/Documents/software/ADIOS2/lib64 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_cxx11_mpi.so.2.6.0 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_cxx11.so.2.6.0 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_core_mpi.so.2.6.0 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_core.so.2.6.0 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_evpath.so /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_ffs.so.1.6.0 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_atl.so.2.2.1 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_dill.so.2.4.1 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_enet.so.1.3.14 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_taustubs.so -Wl,-rpath-link,/p/home/kteferra/Documents/software/ADIOS2/lib64

lastly when compile (AND EXECUTE) an application that uses adios2, you have to make sure you have the following modules loaded
module swap compiler/intel/2019.4.243 compiler/intel/2019.5.281
module swap mpt/2.20 compiler/intelmpi/2019.5.281
module load hdf5-parallel/intel-18.1.163/1.10.5
module load gcc/9.2.0

see the Makefile


ONYX:
1) module swap PrgEnv-cray PrgEnv-intel/6.0.5
2) module load intel/19.0.1.144
3) module load cray-hdf5-parallel/1.10.5.0
4) module load gcc/8.3.0

then git cloned adios2, created subdir adios2-build and cc'ed to it
1) ../../cmake-3.18.2/bin/cmake -DCMAKE_INSTALL_PREFIX=/p/home/kteferra/Documents/software/ADIOS2 -DADIOS2_USE_Fortran=OFF -DADIOS2_USE_MPI=ON -DADIOS2_USE_HDF5=ON ../../ADIOS2
2) make -j 12
3) make install

this created static libraries. Onyx is a cray system and this system by default creates static librarys. compiler flags can be found by running executable adios2-config in ADIOS2/bin
