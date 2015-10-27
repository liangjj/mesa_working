#   @(#)Makefile	1.9   8/7/91
#
#  Makefile for the CLAMS (common Los Alamos mathematical subroutines)
#  portion of MesaLib
#
MYLIB = ../mylib.a
MESALIB = ../mesalib.a
AR = -ar
ARFLAGS = rvu
RANLIB = ranlib

FFLAGS = -c 
GET = sccs get
REL = 

SRCS = \
	alnrel.f\
	bakvec.f \
	balanc.f \
	balbak.f \
	bisect.f \
	bsaxpy.f \
	bscattr.f \
	bscopy.f \
	bsdot.f \
	c9lgmc.f \
	carg.f \
	catan.f \
	cbabk2.f \
	cbabk22.f \
	cbabk3.f \
	cbal.f \
	ccot.f \
	cdiv.f \
	cgamma.f \
	cgeev.f \
	cgevlrv.f \
	clngam.f \
	clnrel.f \
	comqr.f \
	comqr2.f \
	corth.f \
	cpsi.f \
	csevl.f \
	csroot.f \
	dcabs1.f \
	fdump.f \
	figi.f \
	figi2.f \
	gamlim.f \
	gamma.f \
	hqr.f \
	hqr2.f \
	i1mach.f \
	inits.f \
	interv.f \
	imtqlv.f \
	j4save.f \
	lzit.f \
	lzhes.f \
	orthes.f \
	ortran.f \
	psi.f \
	rabs.f \
	r1mach.f \
	r9lgmc.f \
	s88fmt.f \
	scopym.f \
	sgbfa.f \
	sgbsl.f \
	sgeev.f \
	sgmm.f \
	sgmv.f \
	sgtcon.f \
        sgtrmv.f \
	sgttrf.f \
	sgttrs.f \
	sgtsl.f \
	slacon.f \
	slangt.f \
	slassq.f \
	snrm2.f \
	sqrdc.f \
	sqrsl.f \
	sqrst.f \
	srot.f  \
	srotg.f \
	sspev.f \
	ssvdc.f \
	strdi.f \
	strsl.f \
	thtw.f \
	tinvit.f \
	tqlrat.f \
	tridib.f \
	xerclr.f \
	xermsg.f \
	xerprn.f \
	xersve.f \
	xerhlt.f \
	xercnt.f \
	xgetua.f

.f.o:
	$(FC) $(FFLAGS) $<

all: $(MYLIB)
$(MYLIB): $(SRCS:.f=.o) 
	$(AR) $(ARFLAGS) $(MYLIB) $(?:.f=.o)
	$(RANLIB) $(MYLIB)

all: $(MESALIB)
$(MESALIB): $(SRCS:.f=.o) 
	$(AR) $(ARFLAGS) $(MESALIB) $(?:.f=.o) 
	$(RANLIB) $(MESALIB)


