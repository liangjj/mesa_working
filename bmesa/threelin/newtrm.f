c $Header: newtrm.f,v 1.2 92/12/12 09:34:55 bis Exp $
*deck newtrm.f
c***begin prologue     newtrm
c***date written       921211 (yymmdd)
c***revision date             (yymmdd)
c***keywords           
c***author             schneider, barry(lanl)
c***source             @(#)util
c***purpose            set up modified rhs for one-dimensional
c***                   schroedinger equation.
c***
c***                   the equation here is y  + f(x) y = g
c***                   and we are computing g.
c***routines called  
c***end prologue
      subroutine newtrm(x,g,energy,refe,n,last)
      implicit integer (a-z)
      real *8  x, g, energy, refe, k
      dimension x(0:n), g(0:n)
      common /io/ inp, iout
      k=sqrt(refe)
      do 10 i=0,last
         g(i) = g(i) - ( energy - refe )*sin(k*x(i))
   10     continue
      return
      end



