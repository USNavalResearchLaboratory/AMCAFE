# #!/bin/sh

CPP = nvcc
#CPP2 = /opt/petsc/arch-linux2-c-opt/bin/mpic++
CPP2=g++

VIRT_ARCH=compute_60
REAL_ARCH=sm_60


# if want double precision, comment first one and uncomment second 
# (and vice versa if want single precision
#DOUBLEPRECISION=-DUSE_DOUBLES
DOUBLEPRECISION=

CUDAFLAGS=  --gpu-architecture=$(VIRT_ARCH) --gpu-code=$(REAL_ARCH) 

#CUDACOMPILEFLGS=-dc -G 
CUDACOMPILEFLGS=-dc -opt-info inline -Xptxas -O100 

#CPPFLAGS = -std=c++11 -DADIOS2_USE_MPI -isystem /home/kteferra/Documents/research/software/ADIOS2/include
CPPFLAGS = -std=c++11 -DADIOS2_USE_MPI -isystem /usr/local/ADIOS2/include

CPPOPTFLAGS = -g -O5 
#CPPOPTFLAGS = -Wall
#CPPOPTFLAGS = -O3

CPPINCLUDE = -I /opt/petsc/arch-linux2-c-opt/include/ 
METISLIB = /opt/petsc/arch-linux2-c-opt/lib/libmetis.so
#adioslib=-Wl,-rpath,/home/kteferra/Documents/research/software/ADIOS2/lib /home/kteferra/Documents/research/software/ADIOS2/lib/libadios2_cxx11_mpi.so.2.6.0 /home/kteferra/Documents/research/software/ADIOS2/lib/libadios2_cxx11.so.2.6.0 -Wl,-rpath-link,/home/kteferra/Documents/research/software/ADIOS2/lib


#adioslib=-Wl,-rpath,/usr/local/ADIOS2/lib /usr/local/ADIOS2/lib/libadios2_cxx11_mpi.so.2.6.0 /usr/local/ADIOS2/lib/libadios2_cxx11.so.2.6.0 /usr/local/ADIOS2/lib/libadios2_core_mpi.so.2.6.0 /usr/local/ADIOS2/lib/libadios2_core.so.2.6.0 /usr/local/ADIOS2/lib/libadios2_evpath.so /usr/local/ADIOS2/lib/libadios2_ffs.so.1.6.0 /usr/local/ADIOS2/lib/libadios2_atl.so.2.2.1 /usr/local/ADIOS2/lib/libadios2_dill.so.2.4.1 /usr/local/ADIOS2/lib/libadios2_enet.so.1.3.14 /usr/local/ADIOS2/lib/libadios2_taustubs.so /usr/local/ADIOS2/lib/libadios2_c_mpi.so.2.6.0 /usr/local/ADIOS2/lib/libadios2_c.so.2.6.0 -Wl,-rpath-link,/usr/local/ADIOS2/lib 


adioslib=-Wl,-rpath,/usr/local/ADIOS2/lib /usr/local/ADIOS2/lib/libadios2_cxx11_mpi.so.2.6.0 /usr/local/ADIOS2/lib/libadios2_cxx11.so.2.6.0 -Wl,-rpath-link,/usr/local/ADIOS2/lib
#adioslib=-l -rpath,/usr/local/ADIOS2/lib /usr/local/ADIOS2/lib/libadios2_cxx11_mpi.so.2.6.0 /usr/local/ADIOS2/lib/libadios2_cxx11.so.2.6.0 -Wl,-rpath-link,/usr/local/ADIOS2/lib

#SOURCES = $(wildcard *.C)
SOURCES = main1.cu Grid.cu BasePlate.cu TempField.cu SampleOrientation.cu VoxelsCA.cu
OBJECTS = $(SOURCES:.cu=.o)


cafe: $(OBJECTS)
	@rm -f $@
	$(CPP) $(CUDAFLAGS) -dlink $(OBJECTS) -o gpucode.o
	$(CPP2) $(CPPFLAGS) $(CPPOPTFLAGS) -o $@  gpucode.o $^ $(METISLIB) $(adioslib) /usr/local/cuda/lib64/*so  /usr/local/cuda/lib64/*.a
#	$(CPP2) $(CPPFLAGS) $(CPPOPTFLAGS) -o $@  gpucode.o $^ $(METISLIB) $(adioslib) -L/usr/local/cuda/lib64 -lcudart

clean:
	rm -f *.o

#--------------------------------------------------------
# pattern rules for creating objects
.SUFFIXES: .cu  # define the suffixes

%.o : %.cu
	$(CPP) -c $(CPPFLAGS) $(DOUBLEPRECISION) $(CUDAFLAGS) $(CPPOPTFLAGS) $(CUDACOMPILEFLGS) $(CPPINCLUDE) $< -o $@

# pattern rules for creating objects
#--------------------------------------------------------
