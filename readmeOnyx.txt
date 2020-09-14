Onyx is a cray system which means it doesn't work with mpirun or mpiexec. you cannot bring in your own mpi wrappers and must use theirs
and compile with cc

In order to get this to work i had to change from cray compiler to gnu, i did this as

module swap PrgEnv-cray PrgEnv-intel

(note that I initially did the following to use the gnu c++ compiler: "module swap PrgEnv-cray PrgEnv-gnu". However it was shockingly slow on 
a test run in debug mode. On the other hand, the intel compiler is much faster - even faster than when i was using the mpi wrapper based on 
the petsc version i compiled)

I had to bring in metis and compile that as a static library and link that as a static link when compiled code. The
only noteworthy thing about compiling metis is that i needed to keep the metis.h file bits to 32 rather than switch to 64 (that becomes
self evident when compiling metis)

Also, i needed to dynamic link the stdc++ library so that when i build the executable for cafe i needed to have -lstdc++

Also, when run code need to submit with "aprun"

Please note that the compiler "cc" points to (based on the PrgEnv) is in fact a wrapper with MPI. You can look at the version of petsc used
as "module avail cray-petsc". Also if you do "module avail" you'll see what mpich is connected. So look at the makefile to see how the code
is compiled. You do not need to specifiy any linking to MPI b/c its already in the environment based on the PrgEnv. They do this so that it
is much more optimized, and I notice this when using the intel compiler (not the gnu though).

