*deck @(#)tag1.f	5.1  11/6/94
      subroutine tag1(nao,nob,ta,
     1      r,g,buf,lbufso,rabcx,incor,sg,tg,ndf)
cc
cc
      implicit real*8(a-h,o-z)
      dimension ta(*)
      dimension g(*),sg(*),tg(*)
      real*8 buf(lbufso),rabcx(*)
      dimension r(*)
      common /pcpack/ ipkt, nhext, ipkes
      common/io/inp,iout
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
c     nao              no of k active orbitals.
c     nob             total no. of k orbitals
c     lock(nk)        lock(k) starting postion of the kth vector in
c                      the array cm.
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
c
      naonao=nao*nao
      nnao=nao*(nao+1)/2
      naonob=nao*nob
      ntnao=nnao*nnao
      nsnao=nnao*naonao
c
      if(incor.eq.0) then
cccc
c   buffer in rabcx array
cccc
      call iosys('rewind abcx on rwf',0,0,0,' ')
c
      lpass=lbufso/naonob
      npass=(nnao-1)/lpass+1
c
      if(lpass.lt.1) then
        call lnkerr(' m1001: buffer size too small in tag1')
      endif
c
      ix=1
      nni=nnao
c
      do 5 i=1,npass
c
       nm=min(lpass,nni)
       nni=nni-nm
       mk=nm*nao
c
       lread=nm*naonob
c
      call iosys('read real abcx from rwf without rewinding',
     #              lread,buf,0,' ')
c
      mx=1
      jx=ix
      do 1 j=1,ndf
       call ebtc(tg(jx),ta(mx),buf,nao,nob,mk)
       mx=mx+naonob
       jx=jx+nsnao
   1  continue
c
      ix=ix+nao*mk
c
   5  continue
c
      kx=1
      do 11 i=1,ndf
       ix=1
       jx=kx
      do 10 j=1,nnao
       call blfold(sg(ix),tg(jx),nao)
       ix=ix+nnao
       jx=jx+naonao
  10  continue
c
      call beapat(g(kx),sg,nnao)
c
      kx=kx+nnao*nnao
c
  11  continue
c
      else
cccc
c    rabcx array resides in core
cccc
      mk=nnao*nao
      naondf=nao*ndf
      naonob=nao*nob
c
c
      kx=1
      mx=1
      do 100 k=1,ndf
c
        call ebtc(tg,ta(mx),rabcx,nao,nob,mk)
c
        mx=mx+naonob
c
      ix=1
      jx=1
      do 110 j=1,nnao
       call blfold(sg(ix),tg(jx),nao)
       ix=ix+nnao
       jx=jx+naonao
 110  continue
c
      call beapat(g(kx),sg,nnao)
c
      kx=kx+nnao*nnao
c
 100  continue
cc
      end if
c
c
      return
 6000 write (iout,9001)
 9001 format(//'0****** mcg1',6x,' error reading ordered integrals')
      call lnkerr(' ')
      end
