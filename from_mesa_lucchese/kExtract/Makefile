#
# Makefile for programs in kExtract
#

BINDIR = ../bin

PROGS = ExtractData ExtractDipMat ExtractConfig ExtractIP ExtractCIVec

.PHONY : all
all : $(PROGS)

$(PROGS) : %: %.f90 ExtractMod.o
	 $(FC90) $(FFLAGS90) -o $(BINDIR)/$(@) ../kExtract/$(*).f90 ../kExtract/ExtractMod.o

ExtractMod.o : ExtractMod.f90
	 $(FC90) $(FFLAGS90) -c ExtractMod.f90
 
