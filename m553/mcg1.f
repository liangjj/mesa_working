*deck @(#)mcg1.f	5.1  11/6/94
      subroutine mcg1(nk,nob,lokk,lenk,mix,cm,
     $     nf35,r,g,buf,lbufso,rabcx,incor,sg)
c
c***begin prologue     mcg1
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mcg1.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c-----------------------------------------------------------------------
c
c   vectorized version 5/15/86  bhl
c
c
c --- description     this routine makes contributions from a block
c                     of 1-symmetry to the updated integrals g(ijkl).
c                     (ij kl) are stored in canonical order.
c
c                         gt(i,j,k,l)=cm(a,i)*rabcx(a,j,k,l)
c
c                         gt is then symmetrized to form g
c
c
c --- input
c
c     nk              no of k active orbitals.
c     nob             total no. of k orbitals
c     lock(nk)        lock(k) starting postion of the kth vector in
c                      the arrays mix and cm.
c     lenk(nk)        length of the kth vector.
c     mix(--)         indices of vector components.
c     cm(--)          vector components.
c     nf35            fortran no for dataset containing transformed
c                     integrals.
c
c --- working storage
c
c     r(nob,nk)
c
c --- output
c
c     g(2)
c
c-----------------------------------------------------------------------
c
c***references
c
c***routines called    (none)
c
c***end prologue       mcg1
c
      implicit real*8(a-h,o-z)
c
      dimension lokk(2),lenk(2),mix(2),cm(nob,nk)
      dimension g(*),sg(*)
      real*8 buf(lbufso),rabcx(*)
      dimension r(nob,nk)
      common /pcpack/ ipkt, nhext, ipkes
c
      common /io/ inp,iout
c
      nsk=nk*nk
      nnk=nk*(nk+1)/2
c
      if(incor.eq.0) then
cccc
c   buffer in rabcx array
cccc
         call iosys('rewind abcx on rwf',0,0,0,' ')
c
         nr=nk*nob
         lpass=lbufso/nr
         npass=(nnk-1)/lpass+1
         ix=1
         nni=nnk
c
         do 5 i=1,npass
c
            call iosys('read real abcx from rwf without rewinding',
     $           lbufso,buf,0,' ')
            nm=min(lpass,nni)
            nni=nni-nm
            mk=nm*nk
c
            call ebtc(g(ix),cm,buf,nk,nob,mk)
            ix=ix+nk*mk
c
 5       continue
c
c
      else
cccc
c    rabcx array resides in core
cccc
         mk=nnk*nk
c
c
         call ebtc(g,cm,rabcx,nk,nob,mk)
c
c
      end if
c
c
      ix=1
      jx=1
      do 10 i=1,nnk
         call blfold(sg(ix),g(jx),nk)
         jx=jx+nsk
         ix=ix+nnk
 10   continue
c
      call blfold(g,sg,nnk)
c
c
      return
 6000 write (iout,9001)
 9001 format(//'0****** mcg1',6x,' error reading ordered integrals')
      call lnkerr(' ')
      end
