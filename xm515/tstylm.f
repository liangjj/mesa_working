*deck @(#)tstylm.f	1.1  4/25/95
      subroutine tstylm(lmax,n,nlm,ylm,wts,ptlm,scr,mxerr)
c***begin prologue     tstylm.f
c***date written       940315  
c***revision date      11/6/94      
c
c***keywords           spherical harmonics, orthogonality 
c***author             martin, richard(lanl)
c***source             @(#)tstylm.f	1.1   4/25/95
c***purpose            tests numerical representation of spherical
c                      harmonics for orthogonality. 
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       tstylm.f
      implicit none
c     --- input variables -----
      integer n,nlm,lmax
c     --- input arrays (unmodified) ---
      integer ptlm(0:lmax,-lmax:lmax)
      real*8 ylm(n,nlm),wts(n)
c     --- input arrays (scratch) ---
      real*8 scr(n)
c     --- output arrays ---
c     --- output variables ---
      real*8 mxerr
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer l,m,lp,mp
      real*8 ovrlp,sdot
      real*8 zero,one,thresh,diff
c
      parameter (zero=0.0d+00,one=1.0d+00,thresh=1.0d-11)
c
      common/io/inp,iout
c
c     --- test the orthonormality of the spherical harmonics.
      mxerr=zero
      do 40 l=0,lmax
         do 30 m=-l,l
c
            do 20 lp=0,l
               do 10 mp=-lp,lp
                  call vmul(scr,wts,ylm(1,ptlm(lp,mp)),n)
                  ovrlp=sdot(n,scr,1,ylm(1,ptlm(l,m)),1)
                  if(l.eq.lp.and.m.eq.mp) then
                     diff=abs(ovrlp-one)
                     if(diff.gt.thresh) then
c                        write(iout,*) 'l,m,lp,mp,ovrlap',
c     $                                 l,m,lp,mp,ovrlp
                     endif
                  else
                     diff=abs(ovrlp)
                     if(diff.gt.thresh) then
c                        write(iout,*) 'l,m,lp,mp,ovrlap',
c     $                                 l,m,lp,mp,ovrlp
                     endif
                  endif
                  mxerr=max(mxerr,diff)
   10          continue
   20       continue
c
   30    continue
   40 continue
c
c
      return
      end
