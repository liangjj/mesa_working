*deck @(#)srt2e.f	5.1  11/6/94
      subroutine srt2e(values,nnp,ntriang,c,nbf,norbs,t1,t2,numij,
     #                 asort,lnsort,ngroup,
     #                 nmax,nsym,ijgrp,ijadd,kadd,ladd,bftorb,h,
     #                 levfrm,val,lab,
     #                 bin,lenbin,toguga,moval,tunit,ops)
c
c***begin prologue     srt2e
c***date written       yymmdd   (yymmdd)
c***revision date      910725   (yymmdd)
c   25 july     1991   rlm at lanl
c      adding unit to the calling list.  it may be 'zints' or 'tints' depending
c      upon the circumstances of the run.
c
c***keywords
c***author             lengsfield, byron (llnl)
c***source             @(#)srt2e.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       srt2e
c
c
c***purpose: pre-sort integrals to guga-vector order
c
c
c paul saxe                 21 august 1984                  lanl
c
      implicit integer (a-z)
c
      character*(*) ops
      character*(*) tunit
      character*8 unit
      logical toguga
      real*8 moval(numij,ntriang)
      real*8 values(nnp,ntriang),c(nbf,norbs),t1(norbs,norbs)
      real*8 t2(norbs,nbf),asort(*),h(numij)
      real*8 val(lenbin)
      integer lab(lenbin),bin(lenbin)
      integer ijgrp(numij),ijadd(numij),kadd(norbs,nsym)
      integer ladd(norbs,nsym),bftorb(norbs)
c
      common /io/inp,iout
c
c
      unit=tunit
c
      len=numij*numij+numij+5000
      len=min(len,5000000)
      call iosys('open scr as scratch',len,0,0,' ')
c
c
      call iosys('read real "mo one-electron integrals" from '//unit,
     $            numij,h,0,' ')
c
      call trtosq(t1,h,norbs,numij)
c
      do 910 ii=1,norbs
         ip=bftorb(ii)
         call scatter(norbs,t2(1,ip),bftorb,t1(1,ii))
  910 continue
c
      call sqtotr(h,t2,norbs,numij)
c
      call iosys('write real "mo one-electron integrals" to scr',
     $            numij,h,0,' ')
c
      len=numij*numij
      call iosys('create real "mo two-electron integrals" on scr',
     $            len,0,0,' ')
c
c     ----- loop through kl triangles of integrals, transforming -----
c
      call iosys('rewind "mo two-electron integrals" on '//unit,
     $            0,0,0,' ')
c
c     ----- initialize the sort to guga order -----
c
      pt=0
c
c     ----- loop through ij triangles of integrals, transforming -----
c
      maxij=0
      ij=0
      do 31 i=1,norbs
         io=bftorb(i)
         ia=io*(io-1)/2
         do 30 j=1,i
            jo=bftorb(j)
            if(jo.gt.io) then
             ja=jo*(jo-1)/2
             indx=ja+io
            else
             indx=ia+jo
            end if
            ij=ij+1
c
c           ----- check that this triangle of integrals is in core ----
c
            if (ij.gt.maxij) then
               minij=maxij+1
               maxij=min(numij,maxij+ntriang)
               lnread=(maxij-minij+1)*nnp
               call iosys('read real "mo two-electron integrals"'//
     #                    ' from '//unit//' without rewinding',lnread,
     $                    values,0,' ')
            end if
c
            call trtosq(t1,values(1,ij-minij+1),norbs,numij)
c
            do 911 ii=1,norbs
               ip=bftorb(ii)
               call scatter(norbs,t2(1,ip),bftorb,t1(1,ii))
  911       continue
c
            call sqtotr(t1,t2,norbs,numij)
c
            ioff=(indx-1)*numij
            call iosys('write real "mo two-electron integrals" to scr',
     $                  numij,t1,ioff,' ')
c
c
 30      continue
 31   continue
c
c
      return
      end
