# #!/bin/sh

CPP = /opt/petsc/arch-linux2-c-opt/bin/mpic++
CPPFLAGS = -std=c++11
#CPPOOPTTFLAGS = -O2 -Wall
#CPPOOPTTFLAGS = -O2
CPPOOPTTFLAGS = -O3

CPPINCLUDE = -I /opt/petsc/arch-linux2-c-opt/include/ 
METISLIB = /opt/petsc/arch-linux2-c-opt/lib/libmetis.so

#SOURCES = $(wildcard *.C)
SOURCES = main1.C Grid.C BasePlate.C TempField.C Partition.C SampleOrientation.C VoxelsCA.C 
OBJECTS = $(SOURCES:.C=.o)

cafe: $(OBJECTS)
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
