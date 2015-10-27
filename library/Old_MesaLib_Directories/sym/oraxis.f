*deck @(#)oraxis.f	5.1  11/6/94
      subroutine oraxis(maxap3,a,b,natoms,atmchg,ixyz)
      implicit real*8(a-h,o-z)
c
c     an axis of rotation or a principal axis may be aligned with a
c     cartesian axis in one of two ways.  this routine decides which
c     way by successivly applying the following three tests:
c     1-- the third moment of charge should be positive,
c     2-- the sum of the projections of the atomic coordinates on
c     the axis should be positive, and
c     3-- the first atom with a non-zero projection on the axis
c     should have a positive projection on the axis.
c     if rotation is neccessary in order to meet one of these criteria
c     it shall be a 180 degree rotation about the axis defined below:
c     reference      axis of
c     axis         rotation
c     x              y
c     y              z
c     z              x
c
      common/tol/ toler,tol2
      dimension a(maxap3,3), b(maxap3,3), atmchg(1), t(3,3)
      data zero,one,four/0.0d+00,1.0d+00,4.0d+00/
      save zero,one,four
c
c     call rtrace(6horaxis,1)
      pi = four * atan(one)
c
c     test the third moment.
c
      test = or3mom(maxap3,a,atmchg,natoms,ixyz)
      if (abs(test) .lt. tol2) goto 25
      if (test .gt. zero) return
      goto 75
c
c     test the sum of the projections of the atomic coordinates on
c     the ixyz axis.
c
 25   continue
      test = zero
      do 30 iat=1,natoms
 30      test = test + a(iat,ixyz)
         if (abs(test) .lt. toler) goto 50
         if (test .gt. zero) return
         goto 75
c
c     find the first atom with a non-zero projection on the axis.
c
 50      continue
         do 55 iat=1,natoms
            test = a(iat,ixyz)
            if (abs(test) .gt. toler) goto 60
 55      continue
         return
 60      if (test .gt. zero) return
c
c     carry out the neccessary rotation.
c
 75      continue
         i2 = 1 + mod(ixyz,3)
         numatm = natoms + 3
         call rotate(maxap3,a,b,numatm,t,i2,pi)
         call move(maxap3,b,a,numatm)
         return
         end
