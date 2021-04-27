#!/bin/bash

#valgrind --leak-check=yes /opt/petsc/arch-linux2-c-opt/bin/mpiexec -n 24 ./cafe singleScanNW.in
valgrind --tool=callgrind /opt/petsc/arch-linux2-c-opt/bin/mpiexec -n 24 ./cafe singleScanNW.in

exit
