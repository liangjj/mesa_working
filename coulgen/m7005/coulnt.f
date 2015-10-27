*deck @(#)coulnt.f	1.1 9/8/91
c***begin prologue     coulnt
c***date written       910803
c***revision date               (yymmdd)
c***keywords           numerov, integration, coulomb
c***author             schneider, barry(lanl)
c***source             @(#)util
c***purpose            to integrate one dimensional second order
c***                   differential equation for a positive energy
c***                   coulomb functions using three point numerov scheme.
c***
c***references         numerov method is well known and can be found
c***                   in many texts on numerical anaysis. it is well
c***                   described in Kopal's book on numerical analysis.
c***                                         ''
c***                   the equation here is y  + f(x) y = g
c***routines called  
c***end prologue
      subroutine coulnt(psi,f,g,stp,n,dir)
      implicit integer (a-z)
      real *8  psi, stp
      real *8  f, g
      character *(*) dir
      dimension psi(n), f(n), g(n)
      common /io/ inp, iout
c**********************************************************************c
c               start solution using two values of function            c
c               chosen either by series or asymptotic expansion        c
c               these values are assumed to be in the proper           c
c               positions in psi before entering the routine           c
c**********************************************************************c
      if (dir.eq.'forward') then
          do 10 i=3,n
             psi(i)= ( 24.d+00 -10.d+00*stp*stp*f(i-1) )*psi(i-1)
             psi(i) = psi(i) - ( 12.d+00 + stp*stp*f(i-2) )*psi(i-2)
             psi(i) = psi(i) + stp*stp*( g(i) + 10.d+00*g(i-1) + 
     1                                        g(i-2) )
             psi(i) = psi(i)/ ( 12.d+00 + stp*stp*f(i) )
   10     continue
      elseif (dir.eq.'backward') then
          do 20 i=n-2,1,-1
             psi(i)= ( 24.d+00 -10.d+00*stp*stp*f(i+1) )*psi(i+1)
             psi(i) = psi(i) - ( 12.d+00 + stp*stp*f(i+2) )*psi(i+2)
             psi(i) = psi(i) + stp*stp*( g(i+2) + 10.d+00*g(i+1) + 
     1                                        g(i) )
             psi(i) = psi(i)/ ( 12.d+00 + stp*stp*f(i) )
   20     continue
      else
          call lnkerr('error in integration direction')
      endif
      return
      end



