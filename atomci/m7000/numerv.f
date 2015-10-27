*deck @(#)numerv.f	1.1 9/8/91
c***begin prologue     numerv
c***date written       910803
c***revision date               (yymmdd)
c***keywords           numerov, integration
c***author             schneider, barry(lanl)
c***source             @(#)util
c***purpose            to integrate one dimensional second order
c***                   differential equation for a regular solution
c***                   using three point numerov scheme.
c***
c***references         numerov method is well known and can be found
c***                   in many texts on numerical anaysis.
c***                                         ''
c***                   the equation here is y  + f(x) y = g
c***routines called  
c***end prologue
      subroutine numerv(psi,r0,r1,driver,energy,pot,stp,n)
      implicit integer (a-z)
      real *8  psi, driver, energy, stp, pot
      real *8  f, g, r0,r1
      dimension psi(n), pot(n), driver(n)
      common /io/ inp, iout
      f(i) = 2.d+00*( energy - pot(i) )
      g(i) = -2.d+00*driver(i)
c----------------------------------------------------------------------c
c               start solution at origin as regular                    c
c               with second point chosen as 1.d+00                     c
c----------------------------------------------------------------------c
      psi(1)=r0
      psi(2)=r1
c----------------------------------------------------------------------c
c               continue using numerov algorithim                      c
c----------------------------------------------------------------------c
      do 10 i=3,n
         psi(i)= ( 24.d+00 -10.d+00*stp*stp*f(i-1) )*psi(i-1)
         psi(i) = psi(i) - ( 12.d+00 + stp*stp*f(i-2) )*psi(i-2)
         psi(i) = psi(i) + stp*stp*( g(i) + 10.d+00*g(i-1) + g(i-2) )
         psi(i) = psi(i)/ ( 12.d+00 + stp*stp*f(i) )
   10 continue
      return
      end
