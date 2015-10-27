*deck @(#)tdm02.f	5.1  11/6/94
      subroutine tdm02(arc,wt,nlwks,ijadd,kadd,ladd,
     #     orbsym,refwt,refarc,b,refb,ints,c,s,
     #     nrows,norbs,nlevs,orbfrm,nsym,nmax,nwks,nnp,
     #     irowsv,jrowsv,segsv,pagesv,iwtsv,jwtsv,traksv,acoef,bcoef,
     $     lnbuf,ibuf,rbuf,cutoff,ngroup,imngrp,imxgrp,jmngrp,jmxgrp,
     $     sirow,sjrow,ops,dunit)
c
c***begin prologue     tdm02
c***date written       870730   (yymmdd)
c***revision date      870806   (yymmdd)
c     modified 6 august 1987  by pws at lanl
c     modifying to handle more than one block of integrals, and to
c     write out a random diagonalization tape.
c
c***keywords           density matrix
c***author             saxe, paul (lanl)
c***source             @(#)tdm02.f	5.1   11/6/94
c
c***purpose            to form the one- and two-particle density
c                      matrices.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       tdm02
c
      implicit integer (a-z)
c
      character*(*) ops
      character*16 dunit
      integer arc(4,nrows)
      integer wt(4,nrows)
      integer nlwks(nrows)
      integer ijadd(nnp)
      integer kadd(norbs,nsym)
      integer ladd(norbs,nsym)
      integer orbsym(norbs)
      integer refwt(nlevs)
      integer refarc(nlevs)
      integer b(nrows)
      integer refb(nlevs)
      integer irowsv(nlevs)
      integer jrowsv(nlevs)
      integer segsv(nlevs)
      integer pagesv(nlevs)
      integer iwtsv(nlevs)
      integer jwtsv(nlevs)
      integer traksv(nlevs)
      integer imngrp(ngroup)
      integer imxgrp(ngroup)
      integer jmngrp(ngroup)
      integer jmxgrp(ngroup)
      real*8 ints(nmax)
      real*8 c(nwks)
      real*8 s(nwks)
      real*8 acoef(nlevs)
      real*8 bcoef(nlevs)
c
      common /io/ inp,iout
c
c
      call dmlps(arc,wt,nlwks,ijadd,kadd,ladd,orbsym,b,ints,sirow,
     $     sjrow,nrows,norbs,nlevs,orbfrm,nsym,nmax,nwks,nnp,irowsv,
     $     jrowsv,segsv,pagesv,iwtsv,jwtsv,traksv,acoef,bcoef,
     $     c,s,imngrp,imxgrp,jmngrp,jmxgrp,
     $     ngroup,dunit,'"guga transition density matrix"')
c
c
      return
      end
