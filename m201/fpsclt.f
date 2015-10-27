*deck @(#)fpsclt.f	5.2  4/18/95
      subroutine fpsclt(nvar,pool0,pool1,d2var,xi,f,fzero,f1,lambda,
     $                  exit)
c***begin prologue     fpsclt.f
c***date written       yymmdd  
c***revision date      4/18/95      
c
c***keywords           
c***author             binkley, et al. (g82)
c***source             @(#)fpsclt.f	5.2   4/18/95
c***purpose            
c***description
c     
c    
c
c***references
c
c***routines called
c
c***end prologue       fpsclt.f
      implicit none
c     --- input variables -----
      integer nvar,lambda
      logical exit
      real*8 f,fzero
c     --- input arrays (unmodified) ---
      real*8 d2var(nvar)
      real*8 xi(nvar),f1(4)
c     --- input arrays (scratch) ---
c     --- output arrays ---
      real*8 pool0(nvar),pool1(nvar)
c     --- output variables ---
c     --- scratch arrays ---
c     --- local variables ---
      integer inp,iout
      integer i
      real*8 alpha,zero,one,two
c
      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00)
c
      common/io/inp,iout
c
 1000 format(5x,'compute scaling factor')
 1010 format(5x,'after extrapolation, alpha=',e20.10)
c
c
      if(lambda.eq.1) then
         write(iout,1000)
      endif
c
c     --- save the value of the function, and increment the counter.
      f1(lambda)=f
      lambda=lambda+1
c
c     --- determine the scaling.  if we have enough points, determine
c         the scaling from a parabolic fit, otherwise step along this
c         direction one more time.
      if(lambda.le.2) then
         alpha=float(lambda)
      else
c        --- use parabolic fit to determine the scaling.
         alpha=one-(f1(2)-fzero)/(two*(f1(2)+fzero-two*f1(1)))
         write(iout,1010) alpha
      endif
c
c     --- set up pool0 from pool1 in xi-space, and transform to x-space.
      do 40 i=1,nvar
         pool1(i)=pool1(i)*sqrt(d2var(i))
   40 continue
      call vwxs(pool0,pool1,xi,alpha,1,nvar)
      do 50 i=1,nvar
         pool0(i)=pool0(i)/sqrt(d2var(i))
         pool1(i)=pool1(i)/sqrt(d2var(i))
   50 continue
c
c     --- determine if we can exit from the scaling sequence.  if so,
c         set the exit flag.
      if(lambda.le.2) then
         exit=.false.
      else
         exit=.true.
      endif
c
c
      return
      end
