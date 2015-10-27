*deck @(#)orc2v.f	5.1  11/6/94
      subroutine orc2v(maxap3,a,b,natoms,atmchg)
      implicit real*8(a-h,o-z)
c
c     following mulliken's recommendation --jcp, 23, 1997 (1955)--,
c     planar c2v molecules are oriented with the molecular plane
c     coincident with the yz cartesian plane.
c     non-planar c2v molecules are oriented so as to place the plane
c     containing the largest number of atoms not on the z-axis in the
c     yz plane.
c
      common/tol/ toler,tol2
      dimension a(maxap3,3), b(maxap3,3), atmchg(1), t(3,3), e(3)
      data one,two,heavy/1.0d+00,2.0d+00,2.0d+00/
      save one,two,heavy
c
c     call rtrace(6horc2v ,1)
      halfpi = two * atan(one)
      numatm = natoms + 3
c
c     test for planarity.
c
      call orptst(maxap3,a,natoms,ixyz)
      if (ixyz .eq. 2) goto 160
      if (ixyz .eq. 1) goto 180
c
c     the molecule is non-planar.  compare the number of atoms in
c     the two mirror planes.
c
      numyz = 0
      numxz = 0
      do 20 iat=1,natoms
         x = a(iat,1)
         y = a(iat,2)
         onz = sqrt(x*x + y*y)
         if (onz .lt. toler) goto 20
         if (abs(x) .lt. toler) numyz = numyz + 1
         if (abs(y) .lt. toler) numxz = numxz + 1
 20   continue
      if (numyz-numxz) 120, 40,200
c
c     compare the number of heavy atoms in the two mirror planes.
c
 40   continue
      numyz = 0
      numxz = 0
      do 60 iat=1,natoms
         if (atmchg(iat) .le. heavy) goto 60
         x = a(iat,1)
         y = a(iat,2)
         onz = sqrt(x*x + y*y)
         if (onz .lt. toler) goto 60
         if (abs(x) .lt. toler) numyz = numyz + 1
         if (abs(y) .lt. toler) numxz = numxz + 1
 60   continue
      if (numyz-numxz) 120, 80,200
c
c     put the lowest numbered atom in the two planes on the yz plane.
c
 80   continue
      do 100 iat=1,natoms
         x = a(iat,1)
         y = a(iat,2)
         onz = sqrt(x*x + y*y)
         if (onz .lt. toler) goto 100
         if (abs(x) .lt. toler) goto 200
         if (abs(y) .lt. toler) goto 120
 100  continue
      call orplan(maxap3,a,b,atmchg,numatm,e,t,3)
      return
 120  call rotate(maxap3,a,b,numatm,t,3,halfpi)
      call move(maxap3,b,a,numatm)
      goto 200
c
c     the molecule is planar and in the yz plane.
c
 180  call orplan(maxap3,a,b,atmchg,natoms+3,e,t,1)
      return
c
c     the molecule is planar and in the xz plane.
c
 160  call oryz(maxap3,a,b,natoms,atmchg,ixyz)
      return
c
 200  continue
      call oraxis(maxap3,b,a,natoms,atmchg,2)
      return
      end
