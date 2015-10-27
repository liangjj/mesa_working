*deck %W%   %G%
      subroutine omega2(a,b,c,l,thetak,phik,omega,ltop,lmax,
     $                  acoef2)
c***begin prologue     omega2
c***date written       920601  (yymmdd)
c***revision date      yymmdd
c
c***keywords           
c***author             martin, richard (lanl) 
c***source             %W%   %G%
c***purpose            compute type2 angular integrals for ecp evaluation 
c***description
c***references
c                      l.e.mcmurchie and e.r.davidson
c
c***routines called
c                      ylm
c
c***end prologue       omega2
c
c
      implicit integer(a-z)
c
c     ----- arguments unchanged -----
      integer a,b,c,l
      real*8 thetak,phik
      complex*16 acoef2(*)
c     ----- arguments returned -----
      complex*16 omega(0:lmax,0:lmax+a+b+c)
c
c     ----- local variables -----
      complex*16 ylm
      complex*16 t1,t2
c
      real*8 zero
c
      parameter (zero=0.0d+00)
c
c
      do 10 m=0,l
         do 10 lambda=0,l+a+b+c
            omega(m,lambda)=dcmplx(zero,zero)
   10 continue
c
      lambeg=max(l-i-j-k,0)
c
c     ----- parity of lambda must match l+i+j+k -----
      if(mod(l+i+j+k,2).eq.0) then
         if(mod(lambeg,2).ne.0) lambeg=lambeg+1
      else
         if(mod(lambeg,2).eq.0) lambeg=lambeg+1
      endif
c     find pointer into the acoef2 array.
      ind=1
      do 170 m=0,l
         do 160 lambda=lambeg,l+i+j+k,2
c           ----- do the mu=0 term first -----
            omega(m,lambda)=ylm(thetak,phik,lambda,0)
     $                          *acoef2(ind)
            do 150 mu=1,lambda
c              ----- the negative mu terms are accounted
c                    for by the conjugation operation -----
               ind=ind+1
               t1=ylm(thetak,phik,lambda,mu)
               t2=acoef2(ind)
               omega(m,lambda)=omega(m,lambda)
     $                         +t1*conjg(t2)+conjg(t1)*t2
  150       continue
  160    continue
  170 continue
c
c
      return
      end
