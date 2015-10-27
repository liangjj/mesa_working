*deck %W%  %G%
      function qalt (n,l,alpha,rk,t)
      implicit integer(a-z)
c
c     alternating series for q(n,l).......
c
c     ----- function returned -----
      real*8 qalt
c     ----- arguments unchanged -----
      integer n,l
      real*8 alpha,rk,t
c     -----arguments returned -----
c     none
c
c     ----- local variables -----
      real*8 one,two
      real*8 xkp,prefac,xden,term,sum,xc
c
c     ----- common -----
      real*8 dfac,pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      common/dfac/dfac(23)
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
c
      parameter (one=1.0d+00,two=2.0d+00)
c
c
      if (l.eq.0) then
         xkp=one
      else
         xkp=rk**l
      endif
      prefac=sqpi2*xkp*dfac(n+l+1)
     $      /(sqrt((two*alpha)**(n+l+1))*dfac(l+l+3))
      num=l-n+2
      xden=float(l+l+3)
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
