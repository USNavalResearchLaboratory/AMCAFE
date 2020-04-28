#!/bin/bash

export numProcs=24
/opt/petsc/arch-linux2-c-opt/bin/mpiexec -n $numProcs /home/kteferra/Documents/research/projects/AMICME/codes/CA/cafeMPI3D/cafeScale $numProcs

exit
