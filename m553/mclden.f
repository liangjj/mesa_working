*deck @(#)mclden.f	5.1  11/6/94
      subroutine mclden(ipqrs,mk,ml,id,nstep,nij,
     $                  mjk,ntot,idtype,nu,isqr,nblkd,mblkd,lblkd)
c
c***begin prologue     mclden
c***date written       871022   (yymmdd)
c***revision date      900417   (yymmdd)
c
c                      april 17, 1990 rlm at lanl
c                      removing necessity to call as an entry point.
c***keywords
c***author             lengsfield, byron (brl)
c***source             @(#)mclden.f	5.1   11/6/94
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue       mclden
c
      implicit real*8 (a-h,o-z)
c
      integer mblkd(51,2),lblkd(2)
c
      common /io/ inp,iout
c
c
c
      isqr=0
      nblkd=lblkd(ipqrs)
      idtype  =  mblkd(1,nblkd)
      nu      =  mblkd(2,nblkd)
      nj      =  mblkd(8,nblkd)
      nk      =  mblkd(9,nblkd)
c
      if(mjk.ne.0)go to 100
c
c   search for the proper coulomb dm block
c
      ix=9
      do 10 i=1,nj
         ix=ix+1
         mkk  =  mblkd(ix,nblkd)
         ix=ix+1
         mll  =  mblkd(ix,nblkd)
         ix=ix+1
         id   =  mblkd(ix,nblkd)
         ix=ix+1
         ntot =  mblkd(ix,nblkd)
         ix=ix+1
         nij  =  mblkd(ix,nblkd)
         ix=ix+2
         isqr =  mblkd(ix,nblkd)
         if(mkk.eq.mk.and.mll.eq.ml)go to 20
 10   continue
c
      write(iout,15) ipqrs,nblkd,mk,ml
 15   format(/,'0  unable to find the density matrices to ',
     $     'match the coulomb  integrals',/,'  ipqrs  nblk  mk  ml',
     $     '         ',6(2x,i4))
      call lnkerr(' ')
c
 20   if(isqr.eq.1) go to 30
      nstep=1
      return
c
 30   continue
      nstep=-nij
      return
c
 100  continue
c
c   search for the proper exchange dm block
c
      ix=23
      do 110 i=1,nk
         ix=ix+1
         mkk   =  mblkd(ix,nblkd)
         ix=ix+1
         mll   =  mblkd(ix,nblkd)
         ix=ix+1
         id    =  mblkd(ix,nblkd)
         ix=ix+1
         ntot  =  mblkd(ix,nblkd)
         ix=ix+1
         nij   =  mblkd(ix,nblkd)
         ix=ix+2
         isqr  =  mblkd(ix,nblkd)
         if(mkk.eq.mk.and.mll.eq.ml)go to 120
 110  continue
c
      write(iout,115) ipqrs,nblkd,mk,ml
 115  format(/,'0  unable to find the density matrices to ',
     $     'match the exchange integrals',/,'  ipqrs  nblk  mk  ml',
     $     '  mr  ms ',6(2x,i4))
      call lnkerr(' ')
c
 120  if(isqr.eq.1) go to 125
      nstep=1
      return
c
 125  continue
      nstep=nij
      return
c
c
      end
