*deck @(#)orkey.f	5.1  11/6/94
      integer function orkey(maxap3,natoms,a,atmchg,nset,npop,aset)
      implicit real*8(a-h,o-z)
c
c     the "key atom" in a symmetric top molecule is hereby defined to
c     be the lowest numbered atom in the first circular-set.  the
c     first circular-set is that set which is nearest the cartesian
c     xy plane.  if two sets are in the same plane then the inner
c     one takes precedence and the lower atomic numbered one takes
c     precedence next.  if two sets are equidistant from the xy
c     plane, the one with a positive projection on the z-axis takes
c     precedence.
c
      common/tol/ toler,tol2
      dimension a(maxap3,3), atmchg(1), npop(1), nset(1), aset(maxap3,3)
      data zero/0.0d+00/
      save zero
c
c     call rtrace(6horkey ,1)
      call cirset(maxap3,natoms,a,atmchg,3,nset,npop,aset,numset)
      iset = 99
      do 100 jat=1,natoms
         jset = nset(jat)
         p = aset(jat,2)
         if (jset.eq.0 .or. jset.eq.iset) goto 100
         if (iset .ne. 99) goto 20
         iset = jset
         iat = jat
         small = p
         goto 100
c
 20      continue
         test = abs(small) - abs(p)
         if (abs(test) .gt. toler) goto 40
         test = p - small
         if (abs(test) .gt. toler) goto 40
         test = aset(iat,3) - aset(jat,3)
         if (abs(test) .gt. toler) goto 40
         test = aset(iat,1) - aset(jat,1)
 40      if (test .lt. zero) goto 100
 60      iset = jset
         small = p
         iat = jat
 100  continue
      orkey = iat
      return
      end
