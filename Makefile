# #!/bin/sh

CPP = /opt/petsc/arch-linux2-c-opt/bin/mpic++
CPPFLAGS = -std=c++11

#CPPOOPTTFLAGS = -O2 -Wall
CPPOOPTTFLAGS = -Wall -Werror
#CPPOOPTTFLAGS = -O3

CPPINCLUDE = -I /opt/petsc/arch-linux2-c-opt/include/ 
METISLIB = /opt/petsc/arch-linux2-c-opt/lib/libmetis.so

#SOURCES = $(wildcard *.C)
SOURCES = main1.C Grid.C BasePlate.C TempField.C Param.C Partition.C SampleOrientation.C VoxelsCA.C 
SOURCES2 = main1Scale.C Grid.C BasePlate.C TempFieldScale.C Partition.C SampleOrientation.C VoxelsCA.C 
SOURCES3 = main1val1.C Grid.C BasePlateval1.C TempFieldval1.C Partition.C SampleOrientation.C VoxelsCAval1.C
OBJECTS = $(SOURCES:.C=.o)
OBJECTS2 = $(SOURCES2:.C=.o)
OBJECTS3 = $(SOURCES3:.C=.o)


cafe: $(OBJECTS)
	@rm -f $@
	$(CPP) -o $@ $^ $(METISLIB)

cafeScale: $(OBJECTS2)
	@rm -f $@
	$(CPP) -o $@ $^ $(METISLIB)

cafeval1: $(OBJECTS3)
	@rm -f $@
	$(CPP) -o $@ $^ $(METISLIB)

clean:
	rm -f *.o

#--------------------------------------------------------
# pattern rules for creating objects
.SUFFIXES: .C  # define the suffixes

%.o : %.C
	$(CPP) -c $(CPPFLAGS) $(CPPOPTFLAGS) $(CPPINCLUDE) $< -o $@

# pattern rules for creating objects
#--------------------------------------------------------
