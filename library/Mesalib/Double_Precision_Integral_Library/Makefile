#   @(#)Makefile	5.1 11/6/94 
#
#    Makefile for the double precision integral library of MesaLib
#    
#    There are four routines partitioned off into SRCDP.
#    These need to be compiled as double-precision routines on both
#    32 and 64 bit machines.
#    
.SUFFIXES: .f90 .mod $(SUFFIXES)
TOP_DIR = $(MESALIB)
include $(TOP_DIR)/Makefile.inc
GET = sccs get
REL = 
SRCS= \
	dfunc.f \
	dnode.f \
	droot.f \
	dsmit.f

OBJECTS = \
	$(SRCS:.f=.o) 

all: $(DPINTLIB)
$(DPINTLIB): $(OBJECTS) 
	$(AR) $(ARFLAGS) $(DPINTLIB) $(OBJECTS)
	$(RANLIB) $(DPINTLIB)
	$(AR) $(ARFLAGS) $(MESA_LIB) $(OBJECTS)
	$(RANLIB) $(MESA_LIB)

sources: $(SRCS) 
$(SRCS):
	$(GET) $(RELEASE) $@

clean:
	rm -f *.o
	rm -f core
