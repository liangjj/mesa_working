*deck @(#)ldata.f	5.1  11/6/94
      subroutine ldata
      implicit real*8(a-h,o-z)
c
c     load useful constants.
c
      common/dfac/dfac(30)
      common/fact/fac(17),fprod(9,9)
      common/pifac/pi,twopi,fpi,pi3haf,pi5hf2,piquart,sqpi,sqpi2
      common/const/zero,one,two,three,four,five,six,ten
c     data zero/0.0d0/,one/1.0d0/,two/2.0d0/,three/3.0d0/,four/4.0d0/
c     data five/5.0d0/,six/6.0d0/,ten/10.0d0/
c     save zero,one,two,three,four,five,six,ten
c
      zero=0.0d+00
      one=1.0d+00
      two=2.0d+00
      three=3.0d+00
      four=4.0d+00
      five=5.0d+00
      six=6.0d+00
      ten=10.0d+00
c
      pi=four*atan(one)
      twopi=two*pi
      fpi=four*pi
      sqpi=sqrt(pi)
      sqpi2=sqrt(pi/two)
      pi3haf=pi*sqrt(pi)
      pi5hf2=twopi*pi3haf
      piquart=twopi/sqrt(pi5hf2)
      dfac(1)=one
      dfac(2)=one
      dfac(3)=one
      do 10 i=4,30
   10    dfac(i)=dfac(i-2)*float(i-2)
      fac(1)=one
      do 20 i=1,16
   20    fac(i+1)=fac(i)*float(i)
      do 30 l1=1,9
         do 30 k1=1,l1
            k=k1-1
            fprod(k1,l1)=fac(l1+k)/(fac(k1)*fac(l1-k1+1))
   30 continue
c
c
      return
      end
