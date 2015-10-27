#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = $(MESA_HOME)/Time_Propagation/CC_Prop
include $(TOP_DIR)/Makefile.inc
SRCS =	CC_Prop_Main.f90
OBJS = $(SRCS:.f90=.o) 

all: $(PROP_DRIVER)

$(PROP_DRIVER): $(OBJS) $(LIB) 
	$(LD90) $(LDFLAGS) $(OBJS) $(LIB) $(MD_LIB) -o $(PROP_DRIVER)

files: $(SRCS) $(MODSRCS)
	$(CO) $(RCS)/*.f90

clean:
	rm -f *.o *.mod *~
