#######################################################################
## Makefile for FieldML core library
##
##   Mark Cheeseman, NIWA
##   December 10, 2012
#######################################################################

include ../Makefile.inc

FIELDML_LIB = ../LIBS/libfieldml_io.a

OBJS = 	./src/ArrayDataReader.o \
	./src/ArrayDataWriter.o \
	./src/FieldmlIoApi.o \
	./src/FieldmlIoSession.o \
	./src/Hdf5ArrayDataReader.o \
	./src/Hdf5ArrayDataWriter.o \
	./src/InputStream.o \
	./src/OutputStream.o \
	./src/StringUtil.o \
	./src/TextArrayDataReader.o \
	./src/TextArrayDataWriter.o

INCS = -I./src -I../core/src  $(XML2_INC) $(HDF5_INC)
LIBS = $(XML2_LIB) $(HDF5_LIB) $(ZLIB_LIB)
FLAGS = $(INCS) $(CXXFLAGS) -DFIELDML_HDF5_ARRAY

all: $(FIELDML_LIB)

%.o: %.cpp
	$(CXX) $(FLAGS) -c -o $@ $<

$(FIELDML_LIB): $(OBJS)
	ar cr $(FIELDML_LIB) $(OBJS)

clean:
	rm $(OBJS) $(FIELDML_LIB)
