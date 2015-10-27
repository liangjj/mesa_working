*deck @(#)qpow.f	5.1  11/6/94
      function qpow (n,l)
      implicit real*8(a-h,o-z)
c     power series for q(n,l).......
c
      common/argab/argab,expab
      common/dfac/dfac(30)
      common/qstore/q(13,11),alpha,rk,t
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
c
      if (l.eq.0) xkp=1.
      if (l.ne.0) xkp=rk**l
      prefac=expab*xkp/sqrt((2.*alpha)**(n+l+1))
      xnum=l+n-1
      xden=l+l+1
      term=dfac(l+n+1)/dfac(l+l+3)
      sum=term
      xj=0.
   10 xnum=xnum+2.
      xden=xden+2.
      xj=xj+1.
      term=term*t*xnum/(xj*xden)
      sum=sum+term
      if ((term/sum).gt.1.d-14) go to 10
      qpow=prefac*sum
      if (mod((l+n),2).eq.0) qpow=qpow*sqpi2
c    following stmt corrects for truncation error.
      qpow=qpow*(1.+xj*0.4d-14)
c
c
      return
      end
