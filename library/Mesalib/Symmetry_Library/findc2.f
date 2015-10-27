*deck @(#)findc2.f	5.1  11/6/94
      subroutine findc2(maxap3,a,b,aset,npop,nset,atmchg,natoms,itst)
      implicit real*8(a-h,o-z)
c
c     this routine tests for a set of norder c2 axes perpindicular
c     to the principal (cartesian z assumed) symmetry axis.  if
c     found, one of the c2 axes is left coincident with the y
c     cartesian axis.
c
      common/tol/ toler,tol2
      dimension a(maxap3,3), b(maxap3,3), atmchg(1), npop(1), nset(1),
     $     aset(maxap3,3), t(3,3)
      data half,one,two/0.5d0,1.d0,2.d0/
      save half,one,two
c
c     call rtrace(6hfindc2,1)
      numatm = natoms + 3
      halfpi = two * atan(one)
      pi = two * halfpi
c
      call cirset(maxap3,natoms,a,atmchg,3,nset,npop,aset,numset)
c
c     look for circular-sets in the xy plane.  a c2 axis must pass
c     through the point midway between any atom in the set and one
c     other atom in the set.
c
      iattop = natoms - 1
      do 100 iset=1,numset
         do 15 iat=1,natoms
            if (nset(iat) .eq. iset) goto 20
 15      continue
         goto 100
 20      continue
         if (abs(aset(iat,2)) .gt. toler) goto 100
         j1 = iat + 1
         do 80 jat=j1,natoms
            if (nset(jat) .ne. iset) goto 80
            x = (a(iat,1)+a(jat,1)) * half
            y = (a(iat,2)+a(jat,2)) * half
            theta = halfpi
            if (abs(y) .gt. toler) theta = - atan(x/y)
            call rotate(maxap3,a,b,numatm,t,3,theta)
            call rotate(maxap3,b,aset,natoms,t,2,pi)
            call equiv(maxap3,b,aset,atmchg,natoms,itst)
            if (itst .eq. 0) goto 80
            call move(maxap3,b,a,numatm)
            return
 80      continue
 100  continue
c
c     pick an atom in one of the circular sets not in the xy plane.
c     a c2 axis must bisect the angle formed by the projection of the
c     reference atom in the xy plane and the projection of an atom
c     in the set opposed to the reference set.
c
      call cirset(maxap3,natoms,a,atmchg,3,nset,npop,aset,numset)
      do 160 iset=1,numset
         do 120 iat=1,natoms
            if (nset(iat) .eq. iset) goto 140
 120     continue
         itst = 0
         return
 140     proi = aset(iat,2)
         ani  = aset(iat,1)
         disi = aset(iat,3)
         goto 180
 160  continue
      itst = 0
      return
 180  xi = a(iat,1)
      yi = a(iat,2)
      j1 = iset + 1
      do 200 jset=j1,numset
         do 185 jat=1,natoms
            if (nset(jat) .eq. jset) goto 190
 185     continue
         itst = 0
         return
 190     continue
         if (abs(proi+aset(jat,2)) .gt. toler   .or.
     $        abs(ani -aset(jat,1)) .gt. tol2    .or.
     $        abs(disi-aset(jat,3)) .gt. toler)  goto 200
         goto 220
 200  continue
      itst = 0
      return
 220  continue
      do 240 jat=1,natoms
         if (nset(jat) .ne. jset) goto 240
         x = (xi+a(jat,1)) * half
         y = (yi+a(jat,2)) * half
         theta = halfpi
         if (abs(y) .gt. toler) theta = - atan(x/y)
         call rotate(maxap3,a,b,numatm,t,3,theta)
         call rotate(maxap3,b,aset,natoms,t,2,pi)
         call equiv(maxap3,b,aset,atmchg,natoms,itst)
         if (itst .eq. 0) goto 240
         call move(maxap3,b,a,numatm)
         return
 240  continue
      itst = 0
      return
      end
