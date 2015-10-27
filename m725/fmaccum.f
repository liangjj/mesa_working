*deck @(#)fmaccum.f	5.1  11/28/95
      subroutine fmaccum(d,dk,out,nbf)
c***begin prologue     fmaccum.f
c***date written       940524   (yymmdd)  
c***revision date      2/5/95
c   may 13, 1994       russo at lanl
c
c***keywords           form sum needed for DFT gradient
c***author             RUSSO, thomas (lanl)
c***source             @(#)fmaccum.f	5.1   11/28/95
c***purpose            form sum(nu) P(mu,nu)DK(mu,nu)
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       fmaccum.f
      implicit none
c     --- input variables -----
      integer nbf
c     --- input arrays (unmodified) ---
      real*8 d(nbf,nbf),dk(nbf,nbf)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 out(nbf)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer mu,nu
c
c
      call rzero(out,nbf)
      do 1 mu=1,nbf
         do 2 nu=1,nbf
            out(mu)=out(mu)+d(mu,nu)*dk(mu,nu)
 2       continue 
 1    continue 
c
c
      return
      end
