#   %W%   %G%
#
# Makefile for cknew
#
# FC = g77
FFLAGS = $(FFLAGSOZ)

# LD = g77
# LDFLAGS =

BINDIR = ../bin
BLASLIB= $(BLASUSE)
GET = sccs get
REL = 

SRCS = \
	../kcknew/bbtrans.f \
	../kcknew/bforthog.f \
	../kcknew/cgeco.f \
	../kcknew/cgefa.f \
	../kcknew/cgefs.f \
	../kcknew/cgesl.f \
	../kcknew/denom.f \
	../kcknew/dlabad.f \
	../kcknew/dladiv.f \
	../kcknew/dlamc1.f \
	../kcknew/dlamc2.f \
	../kcknew/dlamc3.f \
	../kcknew/dlamc4.f \
	../kcknew/dlamc5.f \
	../kcknew/dlamch.f \
	../kcknew/dlapy2.f \
	../kcknew/dlapy3.f \
	../kcknew/fdump.f \
	../kcknew/fforth1.f \
	../kcknew/fforth2.f \
	../kcknew/hbftrans.f \
	../kcknew/i1mach.f \
	../kcknew/ilaenv.f \
	../kcknew/j4save.f \
	../kcknew/lsame.f \
	../kcknew/obftrans.f \
	../kcknew/ovtrans.f \
	../kcknew/r1mach.f \
	../kcknew/rmachr4.f \
	../kcknew/rdbinsqr.f \
	../kcknew/rdbintri.f \
	../kcknew/top.f \
	../kcknew/wrbinsqr.f \
	../kcknew/xerbla.f \
	../kcknew/xercnt.f \
	../kcknew/xerhlt.f \
	../kcknew/xermsg.f \
	../kcknew/xerprn.f \
	../kcknew/xersve.f \
	../kcknew/xgetua.f \
	../kcknew/zgebak.f \
	../kcknew/zgebal.f \
	../kcknew/zgeev.f \
	../kcknew/zgehd2.f \
	../kcknew/zgehrd.f \
	../kcknew/zhseqr.f \
	../kcknew/zlacgv.f \
	../kcknew/zlacpy.f \
	../kcknew/zladiv.f \
	../kcknew/zlahqr.f \
	../kcknew/zlahrd.f \
	../kcknew/zlange.f \
	../kcknew/zlanhs.f \
	../kcknew/zlarf.f \
	../kcknew/zlarfb.f \
	../kcknew/zlarfg.f \
	../kcknew/zlarft.f \
	../kcknew/zlarfx.f \
	../kcknew/zlascl.f \
	../kcknew/zlaset.f \
	../kcknew/zlassq.f \
	../kcknew/zlatrs.f \
	../kcknew/ztrevc.f \
	../kcknew/zung2r.f \
	../kcknew/zunghr.f \
	../kcknew/zungqr.f 

SRCSV = \
	../kcknew/cknew.for

all: $(BINDIR)/xcknew

$(SRCS:%.f=%.o): %.o: %.f
	$(FC) $(FFLAGS) $<

$(SRCSV:%.for=%.o): %.o: %.for
	sed -e "s/#maxnbfkohn/$(MAXNBFKOHN)/g" <$(*).for | \
	sed -e "s/#maxsmall/$(MAXSMALL)/g" | \
	sed -e "s/#maxbig/$(MAXBIG)/g" | \
	sed -e "s/#maxltop/$(MAXLTOP)/g" | \
	sed -e "s/#maxlmtop/$(MAXLMTOP)/g" | \
	sed -e "s/#maxchan/$(MAXCHAN)/g" >$(*).f
	$(FC) $(FFLAGS) $(*).f
	rm -f $(*).f

$(BINDIR)/xcknew: $(SRCSV:%.for=%.o) $(SRCS:.f=.o)
	$(LD) $(LDFLAGS) $(SRCSV:%.for=%.o) $(SRCS:.f=.o) $(BLASLIB) -o $(BINDIR)/xcknew


sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/cknew.f
	cat $(SRCS) > ../source/cknew.f

print: link
	lpr ../source/cknew.f

clean:
	rm -f *.o core
