#  @(#)Makefile	5.1 11/6/94
#
# Makefile for making Mesalib
# This assumes the following directory system as well as makefiles
#
# You should set your environment variables for the fortran compiler
# and the various flags and libraries needed.
#
TOP_DIR = ./
include $(TOP_DIR)/Makefile.inc
.SUFFIXES: .f90 .mod $(SUFFIXES)

MAKE = make
MFLAGS = 
RCS = RCS

LIB_DIRECTORY = \
	Module 

all: $(LIB_DIRECTORY)

.RECURSIVE: $(LIB_DIRECTORY) 
$(LIB_DIRECTORY): FORCE
	cd $@ ; $(MAKE) $(MFLAGS)

directory: FORCE
	cd Module ; $(MAKE) $(MFLAGS) directory

clean: FORCE
	cd Module ; $(MAKE) $(MFLAGS) clean

FORCE:

