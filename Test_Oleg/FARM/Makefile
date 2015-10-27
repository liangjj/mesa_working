#@(#)Makefile   1.4   8/2/91
#                                                                                                                                 
#  Makefile for the propagator code                                                                                                    
#                                                                                                                    
# This Makefile will compile, and load the libraries necessary to produce the executable farm_x
#
TOP_DIR = $(MESA_HOME)/Test_Oleg
include $(TOP_DIR)/Makefile.inc
SRCS = farm1.f
OBJS = $(SRCS:.f=.o)
all: $(FARM_X)

$(FARM_X): $(OBJS) 
	$(LD90) $(LDFLAGS) $(OBJS) $(MD_LIB) $(GOTOBLASLIB) -o $(FARM_X)

files: $(SRCS)
	$(CO) $(RCS)/*.f

clean:
	rm -f *.o *.mod *~
