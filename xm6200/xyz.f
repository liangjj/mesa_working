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

