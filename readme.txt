20200914:

creating a fork in /opt/cafe that will create branches for others to work on

20200821: compiling ADIOS2 to be able to output HDF5

adios2 needs, HDF5 parallel and cmake

neocortex:
HDF5
1) untarred hdf5 and cd'd to directory
2) did the following command
CC=mpicc ./configure --enable-parallel
make -j 12
make install

* note that mpicc refers to a mpich build already on neocortex.
  another build should have been ok.

CMAKE
1) downloaded cmake source and untarred, cd'd to directory
2) needed the openssl library (had an error first)
sudo apt get install libssl-dev
3) then built cmake as
./bootstrap
sudo make -j 12
make install

ADIOS2
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
I'm not really sure how this happened, but I compiled a static library in onyx whereas the others are dynamic, but this is what i did
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

this created static libraries. again, you can find compiler flags by running executable adios2-config in ADIOS2/bin






20200526:
note that in the makefile in this repository there are links to 
compiler wrapping MPI libraries in a petsc directory and a link to 
metis library. I've typically compiled these libraries by downloading 
petsc zip file, going to directory, and doing the following command (or at least some minor
variation- you can go to the petsc website to see):

./configure --with-cc=gcc --with-cxx=g++ --with-debugging=0 COPTFLAGS='-O3 -march=native -mtune=native' 
CXXOPTFLAGS='-O3 -march=native -mtune=native' FOPTFLAGS='-O3 -march=native -mtune=native'
 --download-mpich --download-metis

That being said, on HPC it probably is better to use an intel compiler rather than GNU and to
use any MPI library that is compiled already on the system by system administrators. That is likely
to be more optimized. In that case, you must find the compiler that is wrapped around an mpi library.
An example of me doing this is with the HPC system ONYX which is a cray system.




This is a git repository for the 3D cafeMPI code. Please refer to the 
comments in the commits to know about the details of each version of the
code. This repository is likely to be cloned by my various machines, such as
ghidora, lochness, gaffney, neocortex, etc. The initial commit of this code
is from neocortex 8/22/2019


Since this is the working repository, I cloned it in ghidorah for keeping a backup by doing

git clone ssh://kteferra@neocortex.nrl.navy.mil:/home/kteferra/Documents/research/projects/AMICME/codes/CA/cafeMPI3D/

The above is only done once, and then everything I commit (in neocortex), i pull it in ghidorah

git pull

this is done in ghidorah at:
/Users/kteferra/Documents/research/projects/AMICME/codes/CA/cafeMPI3D/neocortex



I'm doing most/all of the editing in neocortex so i created a branch for neocortex that 
has the evolution of the code (commits), certainly as of 20191007

-- below is old because i dont push fromo neocortex to anything any more. instead i use
neocortex as the main repository and i pull from other places (i.e., ghidorah) that i set 
up as clones.---

When i push updates from neocortex to godzilla2:
1) make sure godzilla2 is mounted (see line "sudo mount ..." in ~/.bashrc,
   note you need sudo privileges to mount- i havent looked into how to do
   this without sudo yet)

2) type 'git push origin noecortex' (note that first time type: 'git push -u origin neocortex'

note that you probably need to do 'sudo git push ...'
   


--
Note that I added a new main file called main1Scale.C as well as a temperature field file
TempFieldScale.C and TempFieldScale.h. Also, there is a new line in the Makefile to compile
this and a new run file runCafeScale.sh to run this. This is to enable an automated way to 
do a scale test of the simulation in terms of parallelization scalability. This will be done
in HPC. The testing will be the time it takes to run a typical time step. A very large base 
plate is used and then time between time step 1 and time step 2 will be used for comparison
as a function of number of processors.
