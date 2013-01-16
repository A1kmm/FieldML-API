#######################################################################
## Makefile for a couple of FieldML test applications: 
##
## simple -> simplistic test of a couple of basic FieldML functions
##
## write  -> tests the HDF5 and non-HDF5 output FieldML functions 
##
## parallel -> tests the parallel HDF5 option in FieldML
##
##   Mark Cheeseman, NIWA
##   January 7, 2012
#######################################################################

include ../Makefile.inc

FIELDML_LIBS = -L../LIBS -lfieldml -lfieldml_io 

INCS = -I../core/src -I../io/src  $(XML2_INC) $(PHDF5_INC) 
LIBS = $(FIELDML_LIBS) $(XML2_LIB) $(PHDF5_LIB) $(ZLIB_LIB)  
FLAGS = $(INCS) $(CXXFLAGS)

all: phdf5

simple: ./src/SimpleTest.cpp 
	$(CXX) $(CXXFLAGS) -o ../bin/simple.x ./src/SimpleTest.cpp	

write: ./src/fieldml_test.cpp
	$(CXX) $(FLAGS) -o ../bin/write.x ./src/fieldml_test.cpp $(LIBS)

test: ./src/simple_hdf5_test.c
	$(CXX) $(FLAGS) -o ../bin/hdf5_test.x ./src/simple_hdf5_test.c $(LIBS)

parallel: ./src/fieldml_phdf5_test.cpp
	$(MPICXX) $(FLAGS) -o ../bin/parallel.x ./src/fieldml_phdf5_test.cpp $(LIBS)

phdf5: ./src/parallel_hdf5_test_suite.cpp
	$(MPICXX) $(FLAGS) -o ../bin/phdf5.x ./src/parallel_hdf5_test_suite.cpp $(LIBS)

clean:
	rm ../bin/simple.x ../bin/write.x ../bin/hdf5_test.x ../bin/phdf5.x 
