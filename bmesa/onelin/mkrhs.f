c $Header: mkrhs.f,v 1.2 92/12/12 09:34:50 bis Exp $
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
c***                   the equation here is y'' + f(x) y = g
c***                   and we are computing g.
c***routines called  
c***end prologue
      subroutine mkrhs(v,x,g,energy,n,last,type)
      implicit integer (a-z)
      real *8  v, x, g, energy, k
      character *(*) type
      dimension v(0:n), x(0:n), g(0:n)
      common /io/ inp, iout
      call rzero(g(0),n+1)
      k=sqrt(energy)
      if (type.eq.'kohn') then
          do 10 i=0,last
             g(i) = 2.d0*v(i)*sin(k*x(i))
   10     continue
      elseif (type.eq.'one') then
          do 20 i=0,last
             g(i)=-2.d0
   20     continue
      endif
      return
      end



