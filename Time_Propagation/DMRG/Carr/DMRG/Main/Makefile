#   @(#)Makefile	1.4   8/2/91
#
#  Makefile for the propagator code
#
TOP_DIR = ..
include $(TOP_DIR)/Makefile.inc

CPPSRCS = \
	MPDFortran.F90 

all: $(DMRG_DRIVER)

$(DMRG_DRIVER): $(CPPSRCS:.F90=.o) $(DMRG_MOD_LIB) $(DMRG_SOURCES_LIB)
	$(LD90) $(LDFLAGS) -o $(DMRG_DRIVER) $(CPPSRCS:.F90=.o) \
		$(DMRG_MOD_LIB) $(DMRG_SOURCES_LIB) $(MD_LIB)

files: $(SRCS)
	$(CO) $(RCS)/*.F90

clean:
	rm -f *.o 
