*deck @(#)qasy.f	5.1  11/6/94
      function qasy (n,l)
      implicit real*8(a-h,o-z)
c
c     asymptotic form of q(n,l).......
c     valid for arbitrary (n,l).
c
      real*8 zero,one,two,four
      real*8 tenm12,round
c
      common/dfac/dfac(30)
      common/qstore/q(13,11),alpha,rk,t
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
c
      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00,four=4.0d+00)
      parameter (tenm12=1.0d-12,round=0.4d-14)
c
      xkp=rk**(n-2)
      prefac=xkp*sqpi2/sqrt((2.*alpha)**(2*n-1))
c
      sum=one
      to=one
      fac1=l-n+2
      fac2=1-l-n
      xc=one
   10 tn=to*fac1*fac2/(four*xc*t)
      if (tn.eq.zero) go to 20
         sum=sum+tn
         if (abs(tn/sum).lt.tenm12) go to 20
         fac1=fac1+two
         fac2=fac2+two
         xc=xc+one
         to=tn
      go to 10
   20 qasy=prefac*sum
c
c    following stmt corrects for truncation error.
      qasy=qasy*(one+xc*round)
      return
      end
