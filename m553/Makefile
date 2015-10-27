#    @(#)Makefile	5.1  11/6/94
#
# Makefile for M553
#
# Modified 4/12/90; load otherlibs before mesalib.
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
MDLIB =
OTHERLIBS = \
	../library/m330lib.a \
	../library/m820lib.a \
	../library/m830lib.a \
	../library/m902lib.a \
	../library/m903lib.a
GET = sccs get
REL = 

SRCS = \
	../m553/ainex.f \
	../m553/aint.f \
	../m553/bigdav.f \
	../m553/bigdmp.f \
	../m553/bigdrv.f \
	../m553/bighvc.f \
	../m553/bigues.f \
	../m553/blfold.f \
	../m553/chknrm.f \
	../m553/cint.f \
	../m553/corden.f \
	../m553/damps.f \
	../m553/dampt.f \
	../m553/davsdr.f \
	../m553/davson.f \
	../m553/dcan.f \
	../m553/expand.f \
	../m553/fixbuf.f \
	../m553/fixval.f \
	../m553/fock.f \
	../m553/formah.f \
	../m553/gethes.f \
	../m553/getrab.f \
	../m553/guess.f \
	../m553/hesvec.f \
	../m553/hfden.f \
	../m553/hvec.f \
	../m553/intrn.f \
	../m553/jtran.f \
	../m553/ktran.f \
	../m553/linear.f \
	../m553/m553.f \
	../m553/makdia.f \
	../m553/maksml.f \
	../m553/makvec.f \
	../m553/mcabcx.f \
	../m553/mcabix.f \
	../m553/mcaibx.f \
	../m553/mcaugh.f \
	../m553/mcbcku.f \
	../m553/mcbigj.f \
	../m553/mcbigk.f \
	../m553/mccg1.f \
	../m553/mccg1e.f \
	../m553/mccga.f \
	../m553/mccgb.f \
	../m553/mccgc.f \
	../m553/mccgpk.f \
	../m553/mcconv.f \
	../m553/mccor.f \
	../m553/mcdmsq.f \
	../m553/mcdriv.f \
	../m553/mcdxaa.f \
	../m553/mcfden.f \
	../m553/mcfgrd.f \
	../m553/mcflag.f \
	../m553/mcfmtr.f \
	../m553/mcfock.f \
	../m553/mcfrmu.f \
	../m553/mcg1.f \
	../m553/mcg1e.f \
	../m553/mcgcor.f \
	../m553/mcgrad.f \
	../m553/mcgrd.f \
	../m553/mcgupk.f \
	../m553/mcgvca.f \
	../m553/mcgvcb.f \
	../m553/mch2h.f \
	../m553/mchess.f \
	../m553/mchsdr.f \
	../m553/mciiii.f \
	../m553/mciter.f \
	../m553/mcjac.f \
	../m553/mcjcbi.f \
	../m553/mckac.f \
	../m553/mcktac.f \
	../m553/mclden.f \
	../m553/mcledr.f \
	../m553/mcmix.f \
	../m553/mcmxdr.f \
	../m553/mcmxvc.f \
	../m553/mcprti.f \
	../m553/mcrdj.f \
	../m553/mcrdk.f \
	../m553/mcsqdm.f \
	../m553/mcsqf.f \
	../m553/mcstv.f \
	../m553/mcstvt.f \
	../m553/mctdrt.f \
	../m553/mcthii.f \
	../m553/mcthij.f \
	../m553/mctrsp.f \
	../m553/mcvcc.f \
	../m553/mcvec.f \
	../m553/mczero.f \
	../m553/movham.f \
	../m553/mxthij.f \
	../m553/mysort.f \
	../m553/natodr.f \
	../m553/natral.f \
	../m553/newvec.f \
	../m553/outtrn.f \
	../m553/pm553.f \
	../m553/printm.f \
	../m553/printr.f \
	../m553/prthss.f \
	../m553/prtint.f \
	../m553/projct.f \
	../m553/pschmd.f \
	../m553/rdiagh.f \
	../m553/rdsgrd.f \
	../m553/rdsqah.f \
	../m553/rdtgrd.f \
	../m553/rdtrah.f \
	../m553/reordc.f \
	../m553/reordr.f \
	../m553/scdm1.f \
	../m553/scfden.f \
	../m553/sch2h.f \
	../m553/schmdt.f \
	../m553/schvec.f \
	../m553/scthii.f \
	../m553/scthij.f \
	../m553/sorti.f \
	../m553/srtdrt.f \
	../m553/sxthij.f \
	../m553/symord.f \
	../m553/trabcx.f \
	../m553/trlag.f \
	../m553/trn1e.f \
	../m553/trn2e.f \
	../m553/vecout.f \
	../m553/vsamul.f \
	../m553/wtdiag.f \
	../m553/zerodm.f
 
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m553

$(BINDIR)/m553: $(SRCS:.f=.o) $(MESALIB) $(OTHERLIBS)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(OTHERLIBS) $(MESALIB)\
             $(MDLIB) -o $(BINDIR)/m553

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m553.f
	cat $(SRCS) > ../source/m553.f

print: link
	lpr ../source/m553.f

clean:
	rm -f *.o core
