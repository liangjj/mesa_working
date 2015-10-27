*deck %W%   %G%
      subroutine omega1(i,j,k,thetak,phik,omega,acoef1)
c***begin prologue     omega1
c***date written       920601  (yymmdd)
c***revision date      yymmdd
c
c***keywords           
c***author             martin, richard (lanl) 
c***source             %W%   %G%
c***purpose            compute type1 angular integrals for ecp evaluation 
c***description
c***references
c                      l.e.mcmurchie and e.r.davidson,
c                         j.comp.phys., 44, 289(1981).
c
c***routines called
c                      ylm
c
c***end prologue       omega1
c
c
      implicit integer(a-z)
c
c     ----- arguments unchanged -----
      integer i,j,k
      real*8 thetak,phik
      complex*16 acoef1(*)
c     ----- arguments returned -----
      real*8 omega(0:i+j+k)
c     ----- external functions -----
      complex*16 ylm
c     ----- local variables -----
      complex*16 t1,t2,sum
c
      real*8 zero
c
      parameter (zero=0.0d+00)
c
c
      do 10 lambda=0,i+j+k
         omega(lambda)=zero
   10 continue
c
c     ----- parity of i+j+k must match lambda -----
      if(mod(i+j+k,2).eq.0) then
         ls=0
      else
         ls=1
      endif
c
      ind=1
      do 30 lambda=ls,i+j+k,2
c        ----- do the mu=0 term first -----
         sum=ylm(thetak,phik,lambda,0)*acoef1(ind)
         do 20 mu=1,lambda
c           ----- the negative mu contributions are accounted
c                 for by the conjugation -----
            ind=ind+1
            t1=ylm(thetak,phik,lambda,mu)
            t2=acoef1(ind)
            sum=sum +conjg(t1)*t2 +t1*conjg(t2)
   20    continue
         omega(lambda)=real(sum)
   30 continue
c
c
      return
      end
