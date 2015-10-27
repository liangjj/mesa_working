*deck %W%  %G%
      function qcomp (n,l,alpha,rk,t)
      implicit real*8(a-h,o-z)
c
c     this routine controls the computation of q(n,l)......
c     arguments are alpha, xk, and t=xk**2/(4.*alpha).
c     there are no restrictions on the magnitude of t.
c     in order to prevent overflow the result returned in qcomp
c     is actually exp(-t)*q(n,l).
c     qalt and qpow are restricted to (n+l).le.22, (2l).le.20,
c     but this may be increased by simply increasing the dfac array.
c     qasy is valid for arbitrary (n,l).
c
c     ----- function returned -----
      real*8 qcomp
c     ----- arguments unchanged -----
      integer n,l
      real*8 alpha,rk,t
c     ----- arguments returned -----
c     none
c
c     ----- local variables -----
      real*8 tmin(9)
      real*8 fifteen
      real*8 expab
c
      data tmin /31.0d+00,28.0d+00,25.0d+00,23.0d+00,22.0d+00,
     $           20.0d+00,19.0d+00,18.0d+00,15.0d+00/
c
      parameter (fifteen=15.0d+00)
c
c     ----- set bias -----
      expab=exp(argab)
c
c     ----- evaluate appropriate series -----
c     if (((n+l).and.1).eq.0.and.(n.gt.l)) then
c     if ((and(n+l,1).eq.0).and.(n.gt.l)) then
      if ((mod(n+l,2).eq.0).and.(n.gt.l)) then
c        ----- n+l even and n>l -----
         qcomp=qalt(n,l,alpha,rk,t)
      else
c        ----- n+l odd or n>=l -----
         if(n.lt.9) then
            if (t.lt.tmin(n+1)) then
               qcomp=qpow(n,l,alpha,rk,t,expab)
            else
               qcomp=qasy(n,l,alpha,rk,t)
            endif
         else
            if(t.lt.fifteen)then
               qcomp=qpow(n,l,alpha,rk,t,expab)
            else
               qcomp=qasy(n,l,alpha,rk,t)
            endif
         endif
      endif
c
c
      return
      end
