*deck @(#)qcentr.f	5.1  11/6/94
      subroutine qcentr(maxap3,natoms,a,atmchg,v)
      implicit real*8 (a-h,o-z)
c
c     calculate the position of the center of charge and return it
c     in v.
c
      real*8 a(maxap3,3), atmchg(natoms), v(3)
      real*8 zero
c
      parameter (zero=0.d0)
c
c     call rtrace(6hcenter,1)
      do 10 i=1,3
 10      v(i) = zero
      totwt = zero
c
      do 20 iat=1,natoms
         wt = atmchg(iat)
         totwt = totwt + wt
         v(1) = v(1) + wt*a(iat,1)
         v(2) = v(2) + wt*a(iat,2)
         v(3) = v(3) + wt*a(iat,3)
 20   continue
c
      v(1) = v(1) / totwt
      v(2) = v(2) / totwt
      v(3) = v(3) / totwt
c
      return
      end
