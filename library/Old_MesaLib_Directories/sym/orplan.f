*deck @(#)orplan.f	5.1  11/6/94
      subroutine orplan(maxap3,a,b,atmchg,numatm,prmom,praxes,ixyz)
      implicit real*8(a-h,o-z)
c
c     rotate the molecule about cartesian axis ixyz such that the
c     principal axis corresponding to the higher in-plane moment
c     is coincident with the higher priority in-plane cartesian
c     axis (the priority of cartesian axes is z > y > x).
c
      common/tol/ toler,tol2
      dimension a(maxap3,3), b(maxap3,3), atmchg(1), prmom(3),
     $     praxes(3,3), t(3,3)
      data one,two/1.0d+00,2.0d+00/
      save one,two
c
c     call rtrace(6horplan,1)
      halfpi = two * atan(one)
      natoms = numatm - 3
      call secmom(maxap3,natoms,a,atmchg,prmom,praxes)
c
c     this orientation scheme applys only to asymmetric top molecules.
c
      itop = 0
      tst1 = prmom(2) - prmom(3)
      tst2 = prmom(1) - prmom(3)
      tst3 = prmom(1) - prmom(2)
      if (abs(tst1) .lt. tol2) itop = itop + 1
      if (abs(tst2) .lt. tol2) itop = itop + 1
      if (abs(tst3) .lt. tol2) itop = itop + 1
      if (itop .ne. 3) itop = itop + 1
      if (itop .ne. 1) return
c
      i2 = 1 + mod(ixyz,3)
      i3 = 1 + mod(i2  ,3)
c
      goto(20,40,60),ixyz
 20   v2 = praxes(2,3)
      v3 = praxes(3,3)
      goto 80
 40   v2 = praxes(1,3)
      v3 = praxes(3,3)
      goto 80
 60   v2 = praxes(1,2)
      v3 = praxes(2,2)
c
 80   continue
      if (abs(v2) .lt. toler) goto 100
      theta = halfpi
      if (abs(v3) .gt. toler) theta = atan(v2/v3)
      call rotate(maxap3,a,b,numatm,t,ixyz,theta)
      call move(maxap3,b,a,numatm)
c
c     orient the axes properly.  it is assumed that the ixyz axis is
c     properly oriented and its orientation will not be changed.
c
 100  continue
      call orptst(maxap3,a,natoms,izyx)
      if (izyx .ne. 0) goto 120
c
c     non-planar molecules.
c
      iax = 3
      if (ixyz .eq. 3) iax = 2
      call oraxis(maxap3,a,b,natoms,atmchg,iax)
      if (ixyz .eq. 2) call oraxis(maxap3,a,b,natoms,atmchg,2)
      return
c
c     planar molecules ... plane perpindicular to ixyz.
c
 120  continue
      if (izyx .ne. ixyz) goto 140
      call oraxis(maxap3,a,b,natoms,atmchg,i3)
      call oraxis(maxap3,a,b,natoms,atmchg,i2)
      return
c
c     planar molecules ... ixyz included in the plane.
c
 140  continue
      iax = i2
      if (i2 .eq. izyx) iax = i3
      call oraxis(maxap3,a,b,natoms,atmchg,iax)
      return
      end
