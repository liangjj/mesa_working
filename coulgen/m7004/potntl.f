*deck @(#)potntl.f	1.1 9/8/91
c***begin prologue     potntl
c***date written       910128   (yymmdd)
c***revision date               (yymmdd)
c***keywords           potential
c***author             schneider, barry(lanl)
c***source             @(#)m6020
c***purpose            calculate potential for coulomb
c***                   equation in reduced co-ordinates.
c***
c***
c***references
c
c***routines called    util
c***end prologue
      subroutine potntl(v,eta,l,xinv,n,type)
      implicit integer (a-z)
      dimension v(n), xinv(n)
      real *8 v, xinv, eta
      character *(*) type
      common /io/ inp, iout
      if (type.eq.'oscillatory') then
          do 10 i=1,n
             v(i) = 1.d0 - 2.d0*eta*xinv(i) -l*(l+1)*xinv(i)*xinv(i)
   10     continue
      elseif (type.eq.'exponential') then
          do 20 i=1,n
             v(i) = -1.d0 - 2.d0*eta*xinv(i) -l*(l+1)*xinv(i)*xinv(i)
   20     continue
      else
          call lnkerr('error in call to potential')
      endif
      return
      end
