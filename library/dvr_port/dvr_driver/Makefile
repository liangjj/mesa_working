#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M1100
#
#               Suffices
.SUFFIXES: .f90 .mod $(SUFFIXES) 
#
#MESALIB = $(MESA_LIB)/mesalib.a
#POTLIB = $(MESA_LIB)/potlib/potlib_f90.a
#MDLIB = /usr/lib/libcxml.a
DVRLIB = $(MESA_LIB)/dvr_port/dvrlib_f90.a
FDLIB = $(MESA_LIB)/dvr_port/fdlib_f90.a
POTLIB = $(MESA_LIB)/dvr_port/potlib_f90.a
UTILLIB = $(MESA_LIB)/dvr_port/utillib.a
MODDIR = $(MESA_LIB)/dvr_port/Modules
MODLIB = $(MESA_LIB)/dvr_port/modlib_f90.a
BLASLIB = $(MESA_LIB)/blaslib.a
LOCAL_BIN = ../dvr_driver
FC = fort
FFLAGS = -c -O4 -I$(MODDIR)
LD = $(FC)
LDFLAGS =

SRCS = \
	driver.f90  

.f90.o:
	$(FC) $(FFLAGS) $<

all: $(LOCAL_BIN)/dvr
$(LOCAL_BIN)/dvr: $(SRCS:.f90=.o)   
	$(LD) $(LDFLAGS) $(SRCS:.f90=.o)  $(DVRLIB) $(FDLIB) \
					  $(MODLIB) $(POTLIB) $(UTILLIB)\
					  $(BLASLIB)  \
                                          -o $(LOCAL_BIN)/dvr

sources: $(SRCS)
	co -l RCS/*.f90