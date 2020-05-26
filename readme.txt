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
