#######################################################################
## Makefile for FieldML core library
##
##   Mark Cheeseman, NIWA
##   December 10, 2012
#######################################################################

include ../Makefile.inc

FIELDML_LIB = ../LIBS/libfieldml.a

OBJS = 	./src/ErrorContextAutostack.o \
	./src/Evaluators.o \
	./src/fieldml_api.o \
	./src/FieldmlDOM.o \
	./src/FieldmlRegion.o \
	./src/FieldmlSession.o \
	./src/fieldml_structs.o \
	./src/fieldml_write.o \
	./src/ImportInfo.o \
	./src/ObjectStore.o \
	./src/SimpleBitset.o \
	./src/string_const.o \
	./src/String_InternalLibrary.o \
	./src/String_InternalXSD.o

INCS = -I./src $(XML2_INC)
LIBS = $(XML2_LIB)
FLAGS = $(INCS) $(CXXFLAGS)

all: $(FIELDML_LIB)

%.o: %.cpp
	$(CXX) $(FLAGS) -c -o $@ $<

$(FIELDML_LIB): $(OBJS)
	ar cr $(FIELDML_LIB) $(OBJS)

clean:
	rm $(OBJS) $(FIELDML_LIB)
