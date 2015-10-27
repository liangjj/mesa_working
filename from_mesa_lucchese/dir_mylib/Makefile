#    %W%   %G%
#
#    Makefile for the character library of MesaLib
#    
MYLIB = ./mylib.a
AR = -ar
RANLIB = ranlib
# FC = g77
# FFLAGS = -c
FFLAGS = $(FFLAGSC)
GET = sccs get
REL = 

SRCS = \
	./COUL90.f \
	./alnrel.f \
	./d1mach.f \
	./dgamln.f \
	./JWKB.f \
	./RICBES.f \
	./SBESJY.f \
	./c9lgmc.f \
	./carg.f \
	./cgamma.f \
	./cgeco.f \
	./cgefa.f \
	./cgefs.f \
	./cgesl.f \
	./class.f \
	./clngam.f \
	./clnrel.f \
	./csevl.f \
	./gamfun.f \
	./gaussq.f \
	./gbslve.f \
	./gbtql2.f \
	./givens.f \
	./fdump.f \
	./i1mach.f \
	./inits.f \
	./j4save.f \
	./r1mach.f \
	./xerbla.f \
	./xercnt.f \
	./xerhlt.f \
	./xermsg.f \
	./xerprn.f \
	./xersve.f \
	./xgetua.f \
	./zabs.f \
	./zacai.f \
	./zacon.f \
	./zairy.f \
	./zasyi.f \
	./zbesh.f \
	./zbesj.f \
	./zbesy.f \
	./zbinu.f \
	./zbknu.f \
	./zbuni.f \
	./zbunk.f \
	./zdiv.f \
	./zexp.f \
	./zgebak.f \
	./zgebal.f \
	./zgeev.f \
	./zgehd2.f \
	./zgehrd.f \
	./zhseqr.f \
	./zkscl.f \
	./zlacgv.f \
	./zlacpy.f \
	./zladiv.f \
	./zlahqr.f \
	./zlahrd.f \
	./zlange.f \
	./zlanhs.f \
	./zlarf.f \
	./zlarfb.f \
	./zlarfg.f \
	./zlarft.f \
	./zlarfx.f \
	./zlascl.f \
	./zlaset.f \
	./zlassq.f \
	./zlatrs.f \
	./zlog.f \
	./zmlri.f \
	./zmlt.f \
	./zrati.f \
	./zs1s2.f \
	./zseri.f \
	./zshch.f \
	./zsqrt.f \
	./zuchk.f \
	./zunhj.f \
	./zuni1.f \
	./zuni2.f \
	./zunik.f \
	./zunk1.f \
	./zunk2.f \
	./zuoik.f \
	./zwrsk.f \
	./ztrevc.f \
	./zung2r.f \
	./zunghr.f \
	./zungqr.f 
 
.f.o:
	$(FC) $(FFLAGS) $<

$(MYLIB): $(SRCS:.f=.o) 
	$(AR) rvu $(MYLIB) $(?:.f=.o)
	$(RANLIB) $(MYLIB)

sources: $(SRCS) 
$(SRCS):
	$(GET) $(RELEASE) $@


