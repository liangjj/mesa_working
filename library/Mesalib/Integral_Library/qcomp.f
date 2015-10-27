*deck @(#)qcomp.f	5.1  11/6/94
      function qcomp (n,l)
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
      integer and
      real*8 tmin(9)
      real*8 fifteen
      common/dfac/dfac(30)
      common/qstore/q(13,11),alpha,rk,t
      common/argab/argab,expab
c
      data tmin /31.0d+00,28.0d+00,25.0d+00,23.0d+00,22.0d+00,
     $           20.0d+00,19.0d+00,18.0d+00,15.0d+00/
      save tmin
c
      parameter (fifteen=15.0d+00)
c
      argab=-t
      expab=exp(argab)
c     if (((n+l).and.1).eq.0.and.(n.gt.l)) then
      if ((and(n+l,1).eq.0).and.(n.gt.l)) then
c        ----- n+l even and n>l -----
         qcomp=qalt(n,l)
      else
c        ----- n+l odd or n>=l -----
         if(n.lt.9) then
            if (t.lt.tmin(n+1)) then
               qcomp=qpow(n,l)
            else
               qcomp=qasy(n,l)
            endif
         else
            if(t.lt.fifteen)then
               qcomp=qpow(n,l)
            else
               qcomp=qasy(n,l)
            endif
         endif
      endif
c
c
      return
      end
