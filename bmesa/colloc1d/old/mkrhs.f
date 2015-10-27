*deck mkrhs.f
c***begin prologue     mkrhs
c***date written       910803
c***revision date               (yymmdd)
c***keywords           
c***author             schneider, barry(lanl)
c***source             @(#)util
c***purpose            set up rhs for numerov integration of
c***                   one-dimensional schroedinger equation.
c***
c***                   the equation here is y''  + f(x) y = g
c***                   and we are computing g.
c***routines called  
c***end prologue
      subroutine mkrhs(rhs,g,stp,lwr,upr)
      implicit integer (a-z)
      real *8  rhs, g, stp
      character *(*) type
      dimension rhs(*), g(*)
      common /io/ inp, iout
      do 10 i=lwr,upr
         rhs(i) = stp*stp*( g(i+1) +10.d0*g(i) +g(i-1) )
   10 continue
      return
      end



