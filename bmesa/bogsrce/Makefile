#   @(#)Makefile	1.4   8/2/91
#
# Makefile for M8020
#
CC = gcc
CFLAGS = -ansi

LD = gcc
LDFLAGS = -lm
BINDIR = /newdsk/mesa/lam/bogsrce
MESCTLIB = /newdsk/mesa/library/mesctlib.a
VECLIB = /usr/lib/libvec.a
F2CLIB = /newdsk/usr/local/lib/gcc-lib/hppa1.1-hp-hpux9.01/2.6.3/libf2c.a
FLIB = /usr/lib/libf.a
# F2CINC = /newdsk/usr/local/lib/gcc-lib/hppa1.1-hp-hpux9.01/2.6.3/include
# F2CINC = /usr/include
# GET = sccs get

SRCS = gauher2.o \
       cdvd_new.o \
       cdvr2h0.o \
       ch02dvr.o \
       cnvec.o \
       cua2d.o \
       cua3d.o \
       cuta2d.o \
       cuta3d.o \
       eigord.o \
       lanc.o \
       rdiag.o \
       tstovl.o \
       gp3d_bog_barry.o
.f.c.o:
	$(CC) $(CFLAGS) $<

all: $(BINDIR)/a.out

$(BINDIR)/a.out: $(SRCS:.f.c=.o) $(MESCTLIB) $(VECLIB) $(FLIB) $(CLIB) $(F2CLIB)
	$(LD) $(SRCS:.f.c=.o) $(VECLIB) $(MESCTLIB) $(FLIB) $(CLIB) $(F2CLIB) $(LDFLAGS)
