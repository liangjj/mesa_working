*deck @(#)qalt.f	5.1  11/6/94
      function qalt (n,l)
      implicit real*8(a-h,o-z)
c
c     alternating series for q(n,l).......
c
      common/dfac/dfac(30)
      common/qstore/q(13,11),alpha,rk,t
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
c
      parameter (one=1.0d+00,two=2.0d+00)
c
      if (l.eq.0) then
         xkp=one
      else
         xkp=rk**l
      endif
      prefac=sqpi2*xkp*dfac(n+l+1)
     $      /(sqrt((2.*alpha)**(n+l+1))*dfac(l+l+3))
      num=l-n+2
      xden=l+l+3
      term=one
      sum=term
      xc=-one
   10 if (num.eq.0) goto 20
         term=term*float(num)*t/(xden*xc)
         xc=xc-one
         sum=sum+term
         num=num+2
         xden=xden+two
      go to 10
   20 qalt=prefac*sum
c
c
      return
      end
