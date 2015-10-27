#   @(#)Makefile	1.3   8/1/91
#
# Makefile for M551
#
# Modified 4/12/90; load otherlibs before mesalib.
#
FC = f77
FFLAGS = -c

LD = f77
LDFLAGS =

BINDIR = ../bin
MESALIB = ../library/mesalib.a
OTHERLIBS = \
	../library/m330lib.a \
	../library/m820lib.a \
	../library/m830lib.a \
	../library/m902lib.a \
	../library/m903lib.a
GET = sccs get
REL = 

SRCS = \
	../m551/bigdav.f \
	../m551/bigdmp.f \
	../m551/bigdrv.f \
	../m551/bighvc.f \
	../m551/bigues.f \
	../m551/blfold.f \
	../m551/chknrm.f \
	../m551/corden.f \
	../m551/damps.f \
	../m551/dampt.f \
	../m551/davsdr.f \
	../m551/davson.f \
	../m551/dcan.f \
	../m551/expand.f \
	../m551/fixbuf.f \
	../m551/fixval.f \
	../m551/fock.f \
	../m551/formah.f \
	../m551/gethes.f \
	../m551/getrab.f \
	../m551/guess.f \
	../m551/hesvec.f \
	../m551/hfden.f \
	../m551/hvec.f \
	../m551/intrn.f \
	../m551/jtran.f \
	../m551/ktran.f \
	../m551/linear.f \
	../m551/m551.f \
	../m551/makdia.f \
	../m551/maksml.f \
	../m551/makvec.f \
	../m551/mcabcx.f \
	../m551/mcabix.f \
	../m551/mcaibx.f \
	../m551/mcaugh.f \
	../m551/mcbcku.f \
	../m551/mcbigj.f \
	../m551/mcbigk.f \
	../m551/mccg1.f \
	../m551/mccg1e.f \
	../m551/mccga.f \
	../m551/mccgb.f \
	../m551/mccgc.f \
	../m551/mccgpk.f \
	../m551/mcconv.f \
	../m551/mccor.f \
	../m551/mcdmsq.f \
	../m551/mcdriv.f \
	../m551/mcdxaa.f \
	../m551/mcfden.f \
	../m551/mcfgrd.f \
	../m551/mcflag.f \
	../m551/mcfmtr.f \
	../m551/mcfock.f \
	../m551/mcfrmu.f \
	../m551/mcg1.f \
	../m551/mcg1e.f \
	../m551/mcgcor.f \
	../m551/mcgrad.f \
	../m551/mcgrd.f \
	../m551/mcgupk.f \
	../m551/mcgvca.f \
	../m551/mcgvcb.f \
	../m551/mch2h.f \
	../m551/mchess.f \
	../m551/mchsdr.f \
	../m551/mciiii.f \
	../m551/mciter.f \
	../m551/mcjac.f \
	../m551/mcjcbi.f \
	../m551/mckac.f \
	../m551/mcktac.f \
	../m551/mclden.f \
	../m551/mcledr.f \
	../m551/mcmix.f \
	../m551/mcmxdr.f \
	../m551/mcmxvc.f \
	../m551/mcprti.f \
	../m551/mcrdj.f \
	../m551/mcrdk.f \
	../m551/mcsqdm.f \
	../m551/mcsqf.f \
	../m551/mcstv.f \
	../m551/mcstvt.f \
	../m551/mctdrt.f \
	../m551/mcthii.f \
	../m551/mcthij.f \
	../m551/mctrsp.f \
	../m551/mcvcc.f \
	../m551/mcvec.f \
	../m551/mczero.f \
	../m551/movham.f \
	../m551/mxthij.f \
	../m551/mysort.f \
	../m551/natodr.f \
	../m551/natral.f \
	../m551/newvec.f \
	../m551/outtrn.f \
	../m551/pm551.f \
	../m551/printm.f \
	../m551/printr.f \
	../m551/prthss.f \
	../m551/prtint.f \
	../m551/pschmd.f \
	../m551/rdiagh.f \
	../m551/rdsgrd.f \
	../m551/rdsqah.f \
	../m551/rdtgrd.f \
	../m551/rdtrah.f \
	../m551/reordc.f \
	../m551/reordr.f \
	../m551/scdm1.f \
	../m551/scfden.f \
	../m551/sch2h.f \
	../m551/schmdt.f \
	../m551/schvec.f \
	../m551/scthii.f \
	../m551/scthij.f \
	../m551/sorti.f \
	../m551/srtdrt.f \
	../m551/sxthij.f \
	../m551/symord.f \
	../m551/trabcx.f \
	../m551/trlag.f \
	../m551/trn1e.f \
	../m551/trn2e.f \
	../m551/vecout.f \
	../m551/zerodm.f
 
.f.o:
	$(FC) $(FFLAGS) $<

all: $(BINDIR)/m551

$(BINDIR)/m551: $(SRCS:.f=.o) $(MESALIB) $(OTHERLIBS)
	$(LD) $(LDFLAGS) $(SRCS:.f=.o) $(OTHERLIBS) $(MESALIB)\
             -o $(BINDIR)/m551

sources: $(SRCS)
$(SRCS):
	$(GET) $(REL) $@

link: sources	
	rm -f ../source/m551.f
	cat $(SRCS) > ../source/m551.f

print: link
	lpr ../source/m551.f

clean:
	rm -f *.o core
