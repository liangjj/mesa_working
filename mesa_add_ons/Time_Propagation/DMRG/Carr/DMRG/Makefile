#  @(#)Makefile	DMRG
#
#                Makefile for making the DMRG code
#
# You should set your environment variables F90, F90FLAGS and MD_LIB to the 
# fortran compiler, the flags you set and the location of the machine dependent 
# fortran library.  You can define them here as well but I find that less convenient.
#
# There is a grand makefile in the DMRG directory that with the command make will 
# compile and load  everything.  Similarly a make clean will remove all but the sources.
# You can also go into the individual subdirectoies and do individual makes.
#
# The DMRG directory contains Modules, Sources and Main.  The last directory, Main, 
# is where the executable is created and where all the externals get satisfied.
#
# Note that it is important to have defined the MD_LIB or its equivalent as that is
# the library where all of the BLAS and LAPACK routines are assumed to be contained.
#

#
# Define Directories
#
DMRG = ../DMRG
	MAIN = $(DMRG)/Main 
	SOURCES = $(DMRG)/Sources 
	MODULES = $(DMRG)/Modules
#
#
DMRG_ALL = \
	$(MODULES) \
	$(SOURCES) \
	$(MAIN)


all: $(DMRG_ALL) 

.RECURSIVE: $(DMRG_ALL)

$(DMRG_ALL): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

directory: FORCE
	cd $(MODULES) ; $(MAKE) $(MFLAGS) directory
	cd $(SOURCES) ; $(MAKE) $(MFLAGS) directory
	cd $(MAIN) ; $(MAKE) $(MFLAGS) directory

clean: FORCE
	cd $(MODULES) ; $(MAKE) $(MFLAGS) clean
	cd $(SOURCES) ; $(MAKE) $(MFLAGS) clean
	cd $(MAIN) ; $(MAKE) $(MFLAGS) clean

FORCE:
