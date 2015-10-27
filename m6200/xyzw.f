*deck xyzw.f
c***begin prologue     xyzw
c***date written       930802   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           xyzw, link m6200
c***author             schneider, barry (nsf)
c***source             m6200
c***purpose            fill a vector and weight
c***references         
c
c***routines called
c***end prologue       xyz
      subroutine xyzw (pt,wt,x,y,z,wtval)
      implicit integer (a-z)
      real*8 pt, wt, x, y, z, wtval
      dimension pt(3)
      common /io/ inp, iout
      pt(1)=x
      pt(2)=y
      pt(3)=z
      wt=wtval
      return
      end

