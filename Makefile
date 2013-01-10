#######################################################################
## Central Makefile for FieldML 
##
##   Mark Cheeseman, NIWA
##   December 10, 2012
#######################################################################

include Makefile.inc

all: fieldml fieldml_io tests

fieldml: 
	cd core; $(MAKE) all; cd ..

fieldml_io: 
	cd io; $(MAKE) all; cd ..

tests: 
	cd test; $(MAKE) all; cd ..

clean:
	cd core; $(MAKE) clean; cd ..
	cd io; $(MAKE) clean; cd ..
	cd test; $(MAKE) clean; cd ..
