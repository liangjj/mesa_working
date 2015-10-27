*deck @(#)intgl.f	1.1 9/8/91
c***begin prologue     intgl
c***date written       890511   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6005, link 6005, print grid
c***author             schneider, barry (lanl)
c***source             m6005
c***purpose            test integral
c***references         none
c
c***routines called    none
c***end prologue       intgl
      subroutine intgl(vint,vlamda,grid,npnts)
      implicit integer (a-z)
      common /io/ inp, iout
      real *8 grid, rval
      complex *16 vint, vlamda
      dimension vint(2), grid(4,npnts), vlamda(npnts)
      do 10 grpt=1,npnts
         rval= sqrt(grid(1,grpt)*grid(1,grpt)+grid(2,grpt)*grid(2,grpt)+
     1              grid(3,grpt)*grid(3,grpt))
         vint(1)=vint(1)+grid(4,grpt)/(rval*rval)
         vint(2)=vint(2)+grid(4,grpt)*vlamda(grpt)*exp(-rval)/rval
   10 continue
      return
      end
