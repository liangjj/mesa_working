#
# $Header: /msrc/proj/mss/tcgmsg/ipcv4.0/Makefile.proto,v 1.1 1994/02/23 16:54:02 d3g681 Exp $
#
    TOPDIR = /n/srv2/MESA/library/tcgmsg
    LIBDIR = $(TOPDIR)/ipcv4.0
     SHELL = /bin/sh
########################################
# Should not need to modify below here
########################################

# for RCS
COFLAGS = -l

LIBOBJ =  brdcst.o    nodeid.o    setdbg.o                          \
          checkbyte.o llog.o      nxtval.o    shmem.o     synch.o   \
          cluster.o   memalign.o  signals.o                         \
          usleep.o    srmover.o   pbegin.o    snd.o       waitall.o \
          drand48.o   mtime.o     pbeginf.o   sockets.o             \
          error.o     nnodes.o    sema.o      stats.o               \
          xdrstuff.o  waitcom.o   random.o    niceftn.o \
	  mdtob.o     mitob.o     mitod.o   \
	  mdtoi.o     globalop.o  evlog.o   \
          evon.o      pfilecopy.o $(EXTRA)

 INCLUDES = cluster.h  sema.h     signals.h  sndrcvP.h  srftoc.h    \
            farg.h     shmem.h    sndrcv.h   sockets.h  xdrstuff.h  \
            msgtypesc.h           evlog.h    defglobals.h

.SUFFIXES:	
.SUFFIXES: .o .f .c .s

.f.o:
	$(FC) -c $(FFLAGS) $*.f


all:	parallel hello.x test.x testf.x parse hello.p test.p testf.p toplot


libtcgmsg.a:	$(LIBOBJ)
	/bin/rm -f libtcgmsg.a
	ar rc libtcgmsg.a $(LIBOBJ)
	ranlib libtcgmsg.a

hello.x:	hello.o libtcgmsg.a
	$(CC) $(CFLAGS) -o $@ hello.o libtcgmsg.a $(LIBS)	

parallel:	parallel.o libtcgmsg.a
	$(CC) $(CFLAGS) -o $@ parallel.o libtcgmsg.a $(LIBS)

test.x:	test.o checkbyte.o libtcgmsg.a
	$(CC) $(CFLAGS) -o $@ test.o libtcgmsg.a $(LIBS)

testf.x:	testf.o libtcgmsg.a
	$(FC) $(FFLAGS) -o $@ testf.o libtcgmsg.a $(LIBS)

parse:	parse.o strdup.o
	$(CC) -O -o $@ parse.o strdup.o
parse.o:	parse.c
	$(CC) -O -c parse.c

toplot:	toplot.o strdup.o
	$(CC) -g -o $@ toplot.o strdup.o -lplot
toplot.o:	toplot.c
	$(CC) -O -c toplot.c

testpf.x:	testpf.o libtcgmsg.a
	$(FC) $(FFLAGS) -o $@ testpf.o libtcgmsg.a $(LIBS)

$(LIBOBJ) parallel.o test.o testf.o:	$(INCLUDES)

testf.o:	msgtypesf.h
globalop.o:	msgtypesc.h

msgtypesf.h:	msgtf.h
	cp msgtf.h msgtypesf.h


lint:	
	lint $(LINTFLAGS) \
	brdcst.c    error.c    mtime.c    pbegin.c signals.c synch.c \
	checkbyte.c nnodes.c   pbeginf.c  snd.c    test.c \
	cluster.c   llog.c     nodeid.c   sema.c   sockets.c waitall.c \
	usleep.c    memalign.c nxtval.c   setdbg.c stats.c \
	drand48.c   srmover.c  parallel.c shmem.c  strdup.c  strtok.c \
	mdtob.c     mitob.c    mitod.c    mdtoi.c  random.c   \
	xdrstuff.c  evon.c     evlog.c    globalop.c niceftn.c

clean:
	/bin/rm -f parallel hello.x test.x testf.x *.o *.a core a.out
	/bin/rm -f events events.??? log.??? parse toplot TAGS *.log

realclean:	clean
	/bin/rm -f Makefile *~ *.p *.log msgtf.h msgtypesf.h #*#

hello.p test.p testf.p:
	echo "`whoami` `hostname` 4 `pwd`/hello.x /tmp" > hello.p
	echo "`whoami` `hostname` 4 `pwd`/test.x /tmp" > test.p
	echo "`whoami` `hostname` 4 `pwd`/testf.x /tmp" > testf.p



