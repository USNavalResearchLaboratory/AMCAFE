#!/bin/bash

export numProcs=24
/opt/petsc/arch-linux2-c-opt/bin/mpiexec -n $numProcs ./cafeScale $numProcs

exit
