h46268
s 00023/00000/00000
d D 1.1 94/02/16 20:35:13 mesa 1 0
c date and time created 94/02/16 20:35:13 by mesa
e
u
U
f e 0
t
T
I 1
*deck %W%  %G%
c***begin prologue     xyz
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           xyz, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            fill a vector
c***references         
c
c***routines called
c***end prologue       xyz
      subroutine xyz (pt,x,y,z)
      implicit integer (a-z)
      real*8 pt, x, y, z
      dimension pt(3)
      common /io/ inp, iout
      pt(1)=x
      pt(2)=y
      pt(3)=z
      return
      end

E 1
