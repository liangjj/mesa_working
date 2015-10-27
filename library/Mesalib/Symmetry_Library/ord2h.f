*deck @(#)ord2h.f	5.1  11/6/94
      subroutine ord2h(maxap3,a,b,natoms,atmchg,ian)
      implicit real*8(a-h,o-z)
c
c     following mullikens's recommendation --jcp, 23, 1997 (1955)--,
c     planar d2h molecules are oriented with:
c     1-- the molecular plane coincident with the yz cartesian plane.
c     2-- the z axis such that it passes through the greatest number
c     of atoms or the greatest number of bonds if the atom
c     criterion is not decisive.
c
      integer ornax
      dimension a(maxap3,3), b(maxap3,3), atmchg(1), ian(1), t(3,3)
      data one,two/1.0d+00,2.0d+00/
      save one,two
c
c     call rtrace(6hord2h ,1)
c     test for planarity.
c
      call orptst(maxap3,a,natoms,ixyz)
      if (ixyz .eq. 0) return
      if (ixyz .eq. 1) goto 20
c
c     put the molecule in the yz plane.
c
      call oryz(maxap3,a,b,natoms,atmchg,ixyz)
c
c     find the axis which should be z and reorient the molecule
c     if necessary.
c
 20   itst = ornax(maxap3,a,natoms,ian)
      if (itst .ne. 2) return
      numatm = natoms + 3
      halfpi = two * atan(one)
      pi = two * halfpi
      call rotate(maxap3,a,b,numatm,t,1,-halfpi)
      call rotate(maxap3,b,a,numatm,t,3,-pi)
      return
      end
