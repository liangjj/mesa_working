*deck veff.f
c***begin prologue     veff
c***date written       910128   (yymmdd)
c***revision date               (yymmdd)
c***keywords           effective potential
c***author             schneider, barry(lanl)
c***source             @(#)m6020
c***purpose            calculate effective potential on grid
c***
c***
c***references
c
c***routines called    util
c***end prologue
      subroutine veff(v,pot,r,energy,l,n,prnt)
      implicit integer (a-z)
      dimension v(n), pot(n), r(n)
      real *8 v, pot, r, energy
      logical prnt
      character*80 title
      common /io/ inp, iout
      do 10 i=1,n
         pot(i) =  energy - 2.d0*( v(i) + l*(l+1)/(r(i)*r(i)) )
 10   continue
      if (prnt) then
          title='effective potential on grid'
          call prntrm(title,pot,n,1,n,1,iout)
      endif
      return
      end







