*deck @(#)cirset.f	5.1  11/6/94
      subroutine cirset(maxap3,natoms,a,atmchg,ixyz,nset,npop,
     $     aset,numset)
      implicit real*8(a-h,o-z)
c
c     a "circular-set" of atoms is hereby defined as those atoms
c     lying in a plane which have the same atomic number and
c     which are equidistant from some reference axis perpindicular
c     to the plane.  a proper rotation axis generates a circular-set
c     of atoms.
c
c     this routine searches for circular-sets of atoms.
c     ixyz is the cartesian reference axis.
c     nset(i) gives the number of the set which atom i belongs to.
c     set 0 is defined as the set of on-axis atoms.
c     npop(j) is the number of atoms in set j.
c
      common/tol/ toler,tol2
      dimension a(maxap3,3), nset(1), npop(1), aset(maxap3,3),
     $     atmchg(1)
      data zero/0.0d+00/
      save zero
c
c     call rtrace(6hcirset,1)
      i2 = 1 + mod(ixyz,3)
      i3 = 1 + mod(i2,3)
c
c     aset(i,1):  atmchg(i) and a flag.
c     aset(i,2):  the projection of the i'th atom on the reference
c     axis.
c     aset(i,3):  its distance from the reference axis.
c
      do 20 iat=1,natoms
         aset(iat,1) = atmchg(iat)
         aset(iat,2) = a(iat,ixyz)
         q2 = a(iat,i2)
         q3 = a(iat,i3)
         aset(iat,3) = sqrt(q2*q2+q3*q3)
 20   continue
c
c     define set 0.
c
      do 40 iat=1,natoms
         if (abs(aset(iat,3)) .gt. toler) goto 40
         nset(iat) = 0
         aset(iat,1) = zero
 40   continue
c
c     find the remaining sets.
c
      iset = 0
      do 80 iat=1,natoms
         if (aset(iat,1) .eq. zero) goto 80
         iset = iset + 1
         nset(iat) = iset
         npop(iset) = 1
         an = aset(iat,1)
         ap = aset(iat,2)
         ad = aset(iat,3)
         aset(iat,1) = zero
         j1 = iat + 1
         if (iat .eq. natoms) goto 80
         do 60 jat=j1,natoms
            if (abs(aset(jat,1)-an) .gt. tol2  .or.
     $           abs(aset(jat,2)-ap) .gt. toler .or.
     $           abs(aset(jat,3)-ad) .gt. toler) goto 60
            nset(jat) = iset
            npop(iset) = npop(iset) + 1
            aset(jat,1) = zero
 60      continue
 80   continue
      numset = iset
c
c     restore aset(i,1).
c
      do 100 iat=1,natoms
 100     aset(iat,1) = atmchg(iat)
         return
         end
