#   @(#)Makefile	1.4   8/1/91
#
# Makefile for M942
#
FC = $(FC)
FFLAGS = $(FFLAGS)

LD = $(FC)
LDFLAGS =

BINDIR = $(MESA_BIN)
MESALIB = $(MESA_LIB)/mesalib.a
BLIULIB=$(MESA_LIB)/itlib/bliulib.a
MDLIB = $(MDLIB)

SRCS = \
	../m942/eshft.f \
	../m942/frmopt.f \
	../m942/ham.f \
	../m942/hampp.f \
	../m942/hamqp.f \
	../m942/hamqq.f \
	../m942/hfill.f \
	../m942/hqpfill.f \
	../m942/hmul.f \
	../m942/hmulpp.f \
	../m942/hmulqp.f \
	../m942/incmem.f \
	../m942/lstwrt.f \
	../m942/m942.f \
	../m942/optdir.f \
	../m942/optit.f \
	../m942/outmem.f \
	../m942/pvectr.f \
	../m942/rdham.f \
	../m942/rdhqq.f \
	../m942/rdhqp.f \
	../m942/svectr.f \
	../m942/vcomp.f \
	../m942/wrtbufq.f \
	../m942/wrtham.f \
	../m942/wrthqq.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m942

$(BINDIR)/m942: $(SRCS:.f=.o) $(BLIULIB) $(MESALIB) 
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(BLIULIB) $(MESALIB) $(MDLIB) \
	                               -o $(BINDIR)/m942

