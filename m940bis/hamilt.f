*deck  @(#)hamilt.f	1.3 7/30/91
      subroutine hamilt(arc,wt,nlwks,ijadd,kadd,ladd,
     #     orbsym,refwt,refarc,b,refb,ints,c,
     #     nrows,norbs,nlevs,orbfrm,nsym,nmax,nwks,nnp,
     #     irowsv,jrowsv,segsv,pagesv,iwtsv,jwtsv,traksv,acoef,bcoef,
     $     lnbuf,ibuf,rbuf,cutoff,ngroup,imngrp,imxgrp,
     $     jmngrp,jmxgrp,sirow,sjrow,ops,unit,dsk,prtflg,title,prnt)
c
c***begin prologue     hamilt
c***date written       870730   (yymmdd)
c***revision date      910606   (yymmdd)
c
c     modified by bhl to handle p and q spaces correctly and
c     the changes incorporated into latest rlm version 6 june 1991.
c
c     modified 6 august 1987  by pws at lanl
c     modifying to handle more than one block of integrals, and to
c     write out a random diagonalization tape.
c
c***keywords           hamiltonian
c***author             saxe, paul (lanl)
c***source             @(#)hamilt.f	1.1   11/30/90
c
c***purpose            to form, and write out, the hamiltonian matrix.
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       hamilt
c
      implicit integer (a-z)
c
      character*16 unit
      character*(*) prtflg, dsk
      logical prnt
      character*(*) ops, title
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
      integer ibuf(2,lnbuf)
      integer imngrp(ngroup)
      integer imxgrp(ngroup)
      integer jmngrp(ngroup)
      integer jmxgrp(ngroup)
      real*8 ints(nmax)
      real*8 c(nwks)
      real*8 acoef(nlevs)
      real*8 bcoef(nlevs)
      real*8 rbuf(lnbuf)
      real*8 cutoff
      dimension title(3)
c
      common /io/ inp,iout
c
c
c
c     ----- if there is an 'exchange level' set , we only wish the
c           exchange part of the hamiltonina.
c
      levp=intkey(ops,'ci=exchange',0,' ')
c
      call iosys('create integer '//title(1)//' on hamiltonian',
     1           -1,0,0,' ')
c
c     ----- this loop is backwards because we are accumulating
c           the diagonals in c
c
      nonz=0
      n=0
c
      ioff=0
      joff=0
      if (levp.eq.0) then
          if(sirow.eq.0) then
             call iosys('read integer nwksq from '//dsk,1,nwksq,0,' ')
             sirow=2
             sjrow=2
             ioff=nwksq
             joff=nwksq
             call loops(arc,wt,nlwks,ijadd,kadd,ladd,orbsym,b,
     1                  ints,sirow, sjrow,nrows,norbs,nlevs,
     2                  orbfrm,nsym,nmax,nwks,nnp,irowsv,jrowsv,
     3                  segsv,pagesv,iwtsv,jwtsv,traksv,acoef,
     4                  bcoef,rbuf,ibuf,lnbuf,cutoff,n,nonz,c,
     5                  imngrp,imxgrp,jmngrp,jmxgrp,ngroup,unit,
     6                  ioff,joff,ops,title)
             sirow=1  
             sjrow=2
             joff=0
             ioff=nwksq
             call loops(arc,wt,nlwks,ijadd,kadd,ladd,orbsym,b,
     1                  ints,sirow,sjrow,nrows,norbs,nlevs,
     2                  orbfrm,nsym,nmax,nwks,nnp,irowsv,jrowsv,
     3                  segsv,pagesv,iwtsv,jwtsv,traksv,acoef,
     4                  bcoef,rbuf,ibuf,lnbuf,cutoff,n,nonz,c,
     5                  imngrp,imxgrp,jmngrp,jmxgrp,ngroup,unit,
     6                  ioff,joff,ops,title)
             sirow=2
             sjrow=1
             joff=0
             ioff=nwksq
             call loops(arc,wt,nlwks,ijadd,kadd,ladd,orbsym,b,
     1                  ints,sirow,sjrow,nrows,norbs,nlevs,
     2                  orbfrm,nsym,nmax,nwks,nnp,irowsv,jrowsv,
     3                  segsv,pagesv,iwtsv,jwtsv,traksv,acoef,
     4                  bcoef,rbuf,ibuf,lnbuf,cutoff,n,nonz,c,
     5                  imngrp,imxgrp,jmngrp,jmxgrp,ngroup,unit,
     6                  ioff,joff,ops,title)
             sirow=1
             sjrow=1
             ioff=0
             joff=0
             call loops(arc,wt,nlwks,ijadd,kadd,ladd,orbsym,b,
     1                  ints,sirow,sjrow,nrows,norbs,nlevs,
     2                  orbfrm,nsym,nmax,nwks,nnp,irowsv,jrowsv,
     3                  segsv,pagesv,iwtsv,jwtsv,traksv,acoef,
     4                  bcoef,rbuf,ibuf,lnbuf,cutoff,n,nonz,c,
     5                  imngrp,imxgrp,jmngrp,jmxgrp,ngroup,unit,
     6                  ioff,joff,ops,title)
          else
             call loops(arc,wt,nlwks,ijadd,kadd,ladd,orbsym,b,
     1                  ints,sirow,sjrow,nrows,norbs,nlevs,
     2                  orbfrm,nsym,nmax,nwks,nnp,irowsv,jrowsv,
     3                  segsv,pagesv,iwtsv,jwtsv,traksv,acoef,
     4                  bcoef,rbuf,ibuf,lnbuf,cutoff,n,nonz,c,
     5                  imngrp,imxgrp,jmngrp,jmxgrp,ngroup,unit,
     6                  ioff,joff,ops,title)
          endif
      else
          write (iout,10) levp
 10       format (5x,'exchange only, level=',i3)
c
          call exchng(arc,wt,nlwks,ijadd,kadd,ladd,orbsym,b,
     1                ints,sirow,sjrow,nrows,norbs,nlevs,
     2                orbfrm,nsym,nmax,nwks,nnp,irowsv,jrowsv,
     3                segsv,pagesv,iwtsv,jwtsv,traksv,acoef,
     4                bcoef,rbuf,ibuf,lnbuf,cutoff,n,nonz,c,
     5                imngrp,imxgrp,jmngrp,jmxgrp,ngroup,levp,
     6                unit,title)
      endif
c
c      if (sirow.eq.1.and.sjrow.eq.2.and.levp.eq.0) then
c
c        ----- hpq part. now form j>i walks -----
c
c          sirow=2
c          sjrow=1
c          call loops(arc,wt,nlwks,ijadd,kadd,ladd,orbsym,b,
c     1               ints,sirow,sjrow,nrows,norbs,nlevs,
c     2               orbfrm,nsym,nmax,nwks,nnp,irowsv,jrowsv,
c     3               segsv,pagesv,iwtsv,jwtsv,traksv,acoef,
c     4               bcoef,rbuf,ibuf,lnbuf,cutoff,n,nonz,c,
c     5               imngrp,imxgrp,jmngrp,jmxgrp,ngroup,unit,
c     6               ioff,joff,ops)
c      endif
c
c     ----- flush the final buffers -----
c
      if (n.gt.0) then
          nonz=nonz+n
          call iosys('write integer '//title(1)//' to hamiltonian'//
     $               ' without rewinding',2*n,ibuf,0,' ')
          call iosys('write integer '//title(1)//' to hamiltonian'//
     $               ' without rewinding',wptoin(n),rbuf,0,' ')
      endif
c
      call iosys('endfile '//title(1)//' on hamiltonian',0,0,0,' ')

c
      call iosys('write integer '//title(2)//' to hamiltonian',
     1            1,nonz,0,' ')
      call iosys('write real '//title(3)//' to hamiltonian',
     1            nwks,c,0,' ')
      call iosys('write real "h '//title(3)//'" to rwf',nwks,c,0,' ')
      if (prtflg.ne.'minimum') then
         nnpwks=nwks*(nwks+1)/2
         write (iout,20) nnpwks,cutoff,nonz+nwks,
     $                   100.0*(nonz+nwks)/nnpwks
 20      format(5x,'total number of matrix elements:',i12,
     $        /,5x,'effective zero cutoff          :',e12.1,
     $        /,5x,'number exceeding cutoff        :',i12,
     $             '(',f5.1,'%)')
      end if
c
c
      if(prnt) then
         call rdham(rbuf,ibuf,c,lnbuf,nwks,title)
      endif
      call iosys('rewind all on hamiltonian read-and-write',0,0,0,' ')
      return
      end







