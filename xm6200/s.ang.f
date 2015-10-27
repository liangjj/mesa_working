h31390
s 00042/00000/00000
d D 1.1 94/02/16 20:34:56 mesa 1 0
c date and time created 94/02/16 20:34:56 by mesa
e
u
U
f e 0
t
T
I 1
*deck %W%  %G%
c***begin prologue     ang
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           ang, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            fill angular quantities in lebedev quadrature
c***references         
c
c***routines called
c***end prologue       ang
      subroutine ang (pt,cthet,sthet,cphi,sphi,phpt,nleb)
      implicit integer (a-z)
      real*8 pt, cthet, sthet, cphi, sphi, phpt, twopi
      dimension pt(3,nleb), cthet(nleb), sthet(nleb), cphi(nleb)
      dimension sphi(nleb), phpt(nleb)
      common /io/ inp, iout
      data twopi /  6.283185307179586d0  /
      cthet(1)=pt(3,1)
      sthet(1)=sqrt(1.d0-cthet(1)*cthet(1))
      sphi(1)=0.d0
      cphi(1)=1.d0
      phpt(1)=0.d0
      cthet(2)=-1.d0
      sthet(2)=sqrt(1.d0-cthet(2)*cthet(2))
      sphi(2)=0.d0
      cphi(2)=1.d0
      phpt(2)=0.d0
      do 10 i=3,nleb
         cthet(i)=pt(3,i)
         sthet(i)=sqrt(1.d0-cthet(i)*cthet(i))
         phpt(i)=atan2(pt(2,i),pt(1,i))
         if (phpt(i).lt.0.d0) then
             phpt(i)=twopi+phpt(i)
         endif    
         cphi(i)=cos(phpt(i))
         sphi(i)=sin(phpt(i))
   10 continue         
      return
      end

E 1
