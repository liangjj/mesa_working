*deck @(#)potntl.f	1.1 9/8/91
c***begin prologue     potntl
c***date written       910128   (yymmdd)
c***revision date               (yymmdd)
c***keywords           potential
c***author             schneider, barry(lanl)
c***source             @(#)m6020
c***purpose            calculate potential for one dimensional
c***                   schroedinger equation.
c***
c***
c***references
c
c***routines called    util
c***end prologue
      subroutine potntl(v,r0,l,stp,n)
      implicit integer (a-z)
      dimension v(n)
      real *8 v, stp, rval, r0
      common /io/ inp, iout
      if (r0.eq.0.0d+00) then
          r0=1.d-06
      endif
      rval=r0
      do 10 i=1,n
         v(i)=.5d+00*l*(l+1)/(rval*rval)
         rval=rval+stp
   10 continue
      return
      end
