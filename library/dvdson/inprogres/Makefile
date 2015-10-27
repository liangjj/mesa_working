#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the davidson eigenvalue and linear system routines
#
DVDLIB = /usr1/mesa/library/dvdson/dvdlib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FC = f77
FFLAGS = -c 
GET = sccs get
REL = 

SRCS = \
	/usr1/mesa/library/dvdson/cdham.f \
	/usr1/mesa/library/dvdson/cdvd.f \
	/usr1/mesa/library/dvdson/cdvr2h0.f \
	/usr1/mesa/library/dvdson/ch02dvr.f \
	/usr1/mesa/library/dvdson/cindvd.f \
	/usr1/mesa/library/dvdson/clavec.f \
	/usr1/mesa/library/dvdson/clares.f \
	/usr1/mesa/library/dvdson/cnvec.f \
	/usr1/mesa/library/dvdson/creinit.f \
	/usr1/mesa/library/dvdson/crmres.f \
	/usr1/mesa/library/dvdson/ctrial.f \
	/usr1/mesa/library/dvdson/cua2d.f \
	/usr1/mesa/library/dvdson/cua3d.f \
	/usr1/mesa/library/dvdson/cuta2d.f \
	/usr1/mesa/library/dvdson/cuta3d.f \
	/usr1/mesa/library/dvdson/dvr2h0.f \
	exvec.f \
	/usr1/mesa/library/dvdson/eigord.f \
	/usr1/mesa/library/dvdson/frmres.f \
	/usr1/mesa/library/dvdson/h02dvr.f \
	/usr1/mesa/library/dvdson/hamil.f \
	/usr1/mesa/library/dvdson/honv.f \
	/usr1/mesa/library/dvdson/lares.f \
	/usr1/mesa/library/dvdson/lavec.f \
	/usr1/mesa/library/dvdson/lindvd.f \
	/usr1/mesa/library/dvdson/lsolve.f \
	/usr1/mesa/library/dvdson/newvec.f \
	/usr1/mesa/library/dvdson/rdham.f \
	/usr1/mesa/library/dvdson/rdiag.f \
	/usr1/mesa/library/dvdson/reinit.f \
	/usr1/mesa/library/dvdson/renorm.f \
	rpa.f \
	/usr1/mesa/library/dvdson/rpadiag.f \
	/usr1/mesa/library/dvdson/rpainit.f \
	/usr1/mesa/library/dvdson/rpaonv.f \
	/usr1/mesa/library/dvdson/rpavec.f \
	rsdvd.f \
	/usr1/mesa/library/dvdson/setind.f \
	/usr1/mesa/library/dvdson/tstovl.f \
	/usr1/mesa/library/dvdson/ua2d.f \
	/usr1/mesa/library/dvdson/ua3d.f \
	/usr1/mesa/library/dvdson/uta2d.f \
	/usr1/mesa/library/dvdson/uta3d.f \
	/usr1/mesa/library/dvdson/vtrial.f 

.f.o:
	$(FC) $(FFLAGS) $<

all: $(DVDLIB)
$(DVDLIB): $(SRCS:.f=.o)
	$(AR) $(ARFLAGS) $(DVDLIB) $(?:.f=.o)



