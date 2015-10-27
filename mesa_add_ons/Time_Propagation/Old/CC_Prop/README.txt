# You should set your environment variables F90, F90FLAGS and MD_LIB to the fortran compiler, the flags
# you set and the location of the machine dependent fortran.  You can define them here as well but I
# find that less convenient.
#
# There is a grand makefile in the CC_Prop directory that with the command make will compile and load
# everything.  Similarly a make clean will remove all but the sources.  You can also go into the individual
# subdirectoies and do individual makes.
#
# The CC_Prop directory contains the sources, modules and main code.  The DVR_Library directory builds the
# FEDVR library which computes everthing needed to compute the one-body operators in the FEDVR basis as well
# as some other things like the Poisson solver for the two-electron integrals.  The Utilities directory is a hack
# of the old Mesalib directories.  Its mostly all of it except  IOsys.  I like to use these routines over and over again
# as they do lots of good things and rewriting code is not my favorite pastime.  This directory would basically allow
# a user to function independently of BLAS, LAPACK etc. if they needed to since it has fortran versions of all of
# that and a heck of a lot more such as the LANL CLAMS library and the SLATEC library.
#
# This Makefile defines the major directories and then recursively goes into each of them in the proper order
# and compiles and loads what is needed.  The last directory, Main, is where the executable is created and where all the
# externals get satisfied.
#
