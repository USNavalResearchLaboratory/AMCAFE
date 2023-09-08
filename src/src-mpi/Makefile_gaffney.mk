##!/bin/sh
CPP = mpiicpc #/p/app/intel/parallel_studio_xe_2019_update5/impi/2019.5.281/intel64/bin/mpiicpc
CPPFLAGS = -std=c++11 -DADIOS2_USE_MPI -isystem /p/home/kteferra/Documents/software/ADIOS2/include
#CPPOPTTFLAGS = -Wall -Werror
CPPOPTFLAGS = -O3
CPPINCLUDE = -I /p/home/kteferra/Documents/software/metis-5.1.0/include -I /p/home/kteferra/Documents/software/gsl-2.7.1/include
#CPPINCLUDE = -I /p/app/intel/parallel_studio_xe_2019_update5/compilers_and_libraries/linux/mpi/intel64/include -I /p/home/kteferra/Documents/software/metis-5.1.0/include

METISLIB = /p/home/kteferra/Documents/software/metis-5.1.0/build/Linux-x86_64/libmetis/libmetis.a

adioslib = -Wl,-rpath,/p/home/kteferra/Documents/software/ADIOS2/lib64 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_cxx11_mpi.so.2.6.0 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_cxx11.so.2.6.0 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_core_mpi.so.2.6.0 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_core.so.2.6.0 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_evpath.so /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_ffs.so.1.6.0 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_atl.so.2.2.1 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_dill.so.2.4.1 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_enet.so.1.3.14 /p/home/kteferra/Documents/software/ADIOS2/lib64/libadios2_taustubs.so -Wl,-rpath-link,/p/home/kteferra/Documents/software/ADIOS2/lib64 -lz -ldl -lm

GSLLIB = -Wl,-rpath,/p/home/kteferra/Documents/software/gsl-2.7.1/lib -L/p/home/kteferra/Documents/software/gsl-2.7.1/lib/ -lgsl -lgslcblas -lm

modulesload =module swap compiler/intel/2019.4.243 compiler/intel/2019.5.281; \
        module swap mpt/2.20 compiler/intelmpi/2019.5.281; \
        module load hdf5-parallel/intel-18.1.163/1.10.5; \
        module load gcc/9.2.0;

#SOURCES = $(wildcard *.C)
SOURCES = main1.C Grid.C BasePlate.C TempField.C Partition.C SampleOrientation.C VoxelsCA.C Utilities.C
OBJECTS = $(SOURCES:.C=.o)


cafe: $(OBJECTS)
	@rm -f $@
	$(modulesload) $(CPP) -o $@ $^ $(METISLIB) $(adioslib) $(GSLLIB)
clean:
	rm -f *.o

#--------------------------------------------------------
# pattern rules for creating objects
.SUFFIXES: .C  # define the suffixes

%.o : %.C
	$(modulesload) $(CPP) -c $(CPPFLAGS) $(CPPOPTFLAGS) $(CPPINCLUDE) $< -o $@

# pattern rules for creating objects
#--------------------------------------------------------
