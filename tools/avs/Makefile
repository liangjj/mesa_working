#
# Makefile for mesareader2
#

INC_FILE=$(ROOT)/usr/avs/include/Makeinclude
include $(INC_FILE)

CCHEMLIB= -lchem_c
FCHEMLIB= -lchem_f
AVS_LIBS = $(ROOT)/usr/avs/lib
INC = $(ROOT)/usr/include
BASELIBS=-lgeom -lutil -lm $(LASTFLIBS)
CMSI_inc=-I/usr/avs/MSI/include
CMSI_lib= -L/usr/avs/MSI/library -lmsichem
F77_FLOWLIBS=-L$(AVS_LIBS) $(CMSI_lib) $(FCHEMLIB) $(CCHEMLIB) -lflow_f $(BASELIBS) 
F77_SIMLIBS=-L$(AVS_LIBS) -lsim_f $(BASELIBS)
F77_BIND = $(ROOT)/usr/avs/bin/f77_binding
F77_INC = $(ROOT)/usr/avs/include
F77FLAGS=$(FFLAGS) +E1 +ppu -I/usr/include -I$(F77_INC)
LIBCHEM_INC=-I$(ROOT)/usr/avs/chemistry/include -I$(ROOT)/usr/avs/include
CHEM_INC=$(ROOT)/usr/avs/include/chemistry
CFLOWLIBS=  -L$(AVS_LIBS) $(CMSI_lib) $(CCHEMLIB) -lflow_c $(BASELIBS) 
CSIMLIBS=  -L$(AVS_LIBS) -lsim_c $(CCHEMLIB) $(BASELIBS) 
CFLAGS= $(AVS_INC) $(AOPTCFLAGS) $(LOCAL_CFLAGS) -I$(CHEM_INC) $(CMSI_inc)

all: mesareader2

mesareader2:	mesareader2.o myftocfoo.o mesalib.a
	$(FC) -o mesareader2  mesareader2.o myftocfoo.o mesalib.a $(F77_FLOWLIBS) 

mesareader2.o:	mesareader2.f
	$(FC) $(F77FLAGS) -c mesareader2.f

molvoyeur:	molvoyeur.o
	$(CC) -o molvoyeur molvoyeur.o $(CFLOWLIBS)

.c.o:	
	$(CC) $(CFLAGS) -c $<

