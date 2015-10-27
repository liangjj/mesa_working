c $Header: phase.f,v 1.2 92/12/12 09:35:03 bis Exp $
*deck phase.f
c***begin prologue     phase
c***date written       910128   (yymmdd)
c***revision date               (yymmdd)
c***keywords           
c***author             schneider, barry(lanl)
c***source             @(#)m6020
c***purpose            calculate phase shift and display solution
c***
c***
c***references
c
c***routines called    util
c***end prologue
      subroutine phase(psi,x,energy,n)
      implicit integer (a-z)
      dimension psi(0:n), x(0:n)
      real *8 psi, x, energy, amplt, pshft, k
      common /io/ inp, iout
      write(iout,*) '      the solution vector'
      write(iout,1) (psi(i),i=0,n)
      k=sqrt(energy)
      amplt=psi(n)/cos(k*x(n))
      write(iout,*) '     the ampltude of the cosine = ',amplt
      pshft=atan(amplt)
      write(iout,*) '     the phase shift = ',pshft
      return
    1 format(5x,5e15.8)
      end
