*deck %W%  %G%
      function qasy (n,l,alpha,rk,t)
      implicit real*8(a-h,o-z)
c
c     asymptotic form of q(n,l).......
c     valid for arbitrary (n,l).
c
c
c     ----- function returned -----
      real*8 qasy
c     ----- arguments unchanged -----
      integer n,l
      real*8 alpha,rk,t
c
c     ----- local variables -----
      real*8 zero,one,two,four
      real*8 tenm12,round
      real*8 xkp,prefac,sum,to,xc,tn,fac1,fac2
c     ----- common -----
      real*8 pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
c
      parameter (zero=0.0d+00,one=1.0d+00,two=2.0d+00,four=4.0d+00)
      parameter (tenm12=1.0d-12,round=0.4d-14)
c
c
      xkp=rk**(n-2)
      prefac=xkp*sqpi2/sqrt((two*alpha)**(2*n-1))
c
      sum=one
      to=one
      fac1=float(l-n+2)
      fac2=float(1-l-n)
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
