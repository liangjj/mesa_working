#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
.SUFFIXES: .f90 .mod $(SUFFIXES)

LD = $(FC)
LDFLAGS = 
RCS = RCS
CO = co -l
TEST = /home/mesa/mesa/Time_Propagation/Coupled_Channel_TD_Propagator/Test
SRCS =	factl.f90 \
	dfactl.f90 \
	legend.f90
.f90.o:
	$(FC) $(FFLAGS) -I$(MODULES) $<

all: $(TEST)
$(TEST): $(SRCS:.f90=.o)

files: $(SRCS) $(MODSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o 
