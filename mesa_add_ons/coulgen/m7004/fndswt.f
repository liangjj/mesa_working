*deck @(#)fndswt.f	1.1 9/8/91
c***begin prologue     fndswt
c***date written       910128   (yymmdd)
c***revision date               (yymmdd)
c***keywords           potential
c***author             schneider, barry(lanl)
c***source             @(#)m6020
c***purpose            calculate switching point for spline
c***
c***
c***references
c
c***routines called    util
c***end prologue
      subroutine fndswt(x,rswtch,n,ptbeg)
      implicit integer (a-z)
      real *8 x, rswtch
      dimension x(n)
      do 10 pt=1,n
         if (x(pt).gt.rswtch) go to 20
   10 continue
   20 ptbeg=pt-1     
      return
      end



