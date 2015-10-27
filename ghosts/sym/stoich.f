*deck stoich
      subroutine stoich(maxap3,iatflg,ian,a,natoms,length,ia,ixyz,lin)
      implicit real*8(a-h,o-z)
c
c     this routine determines the stoichiometric formula for the
c     atoms in the current symmetric sub-space (i.e. those atoms
c     whose flags are set to 2) and places it in the framework
c     group buffer.
c     lin is used to indicate that:
c     lin=0  the current symmetric sub-space is not linear.
c     lin=1  the current symmetric sub-space is linear.
c     lin=2  the current symmetric sub-space is linear and
c     contains a unique central point.
c     lin=3  initialize idx.
c
      common/tol/ toler,tol2
      dimension iatflg(1), ian(1), a(maxap3,3), ia(1), idx(92), jdx(92)
      character*(*) ia
      character*2 el(92), itoc
      data jdx/ 6,  1, 89, 47, 13, 18, 33, 85, 79,  5,
     $     56, 4, 83, 35, 20, 48, 58, 17, 27, 24,
     $     55, 29, 66, 68, 63,  9, 26, 87, 31, 64,
     $     32,  2, 72, 80, 67, 53, 49, 77, 19, 36,
     $     57,  3, 71, 12, 25, 42,  7, 11, 41, 60,
     $     10, 28,  8, 76, 15, 91, 82, 46, 61, 59,
     $     84, 78, 88, 37, 75, 45, 86, 44, 16, 51,
     $     21, 34, 14, 62, 50, 38, 73, 65, 43, 52,
     $     90, 22, 81, 69, 92, 23, 74, 54, 39, 70,
     $     30, 40/
      data zero,one,big/0.d0,1.d0,1.0d20/
c
c     call rtrace(6hstoich,1)
      call fillel(1,92,el)
      if (lin .eq. 3) goto 100
      if (lin .gt. 0) goto 200
      call append('(',ia,length)
      do 40 i=1,92
         iatnum = idx(i)
         num = 0
         do 20 iat=1,natoms
            if (iatflg(iat) .ne. 2) goto 20
            if (ian(iat) .ne. iatnum) goto 20
            num = num + 1
            iatflg(iat) = 1
 20      continue
         if (num .eq. 0) goto 40
         call append(el(iatnum),ia,length)
         call append(itoc(num),ia,length)
 40   continue
      call append(')',ia,length)
      call append(',',ia,length)
      return
c
c     initialize idx which governs the order of the atomic symbols
c     in a stoichiometric formula.  the hill system is used:  for
c     carbon containing compounds the order is carbon, hydrogen,
c     with remaining elements in alphabetical order.  for non-
c     carbonaceous species the order is strictly alphabetical.
c     the order given in the data statement assumes the prescence
c     of carbon.
c
 100  continue
      do 110 i=1,92
         idx(i) = jdx(i)
 110  continue
      do 120 iat=1,natoms
         if (ian(iat) .eq. 6) return
 120  continue
      do 130 i=1,12
 130     idx(i) = idx(i+2)
         idx(13) = 6
         do 140 i=14,30
 140        idx(i) = idx(i+1)
            idx(31) = 1
            return
c
c     for linear subspaces the atoms are listed in sequential order
c     along the axis.  if the molecule has a center of symmetry
c     (lin=2) then the central point will be represented by a dot.
c
 200        continue
            call append('(',ia,length)
c
c     count the flagged atoms.
c
            numfla = 0
            do 220 iat=1,natoms
               if (iatflg(iat) .ne. 2) goto 220
               numfla = numfla + 1
 220        continue
c
c     scan the axis from negative to positive.
c
 240        continue
            oldz = -big
            zmin = big
 250        iflag = 0
            do 260 iat=1,natoms
               if (iatflg(iat) .ne. 2) goto 260
               iflag = 1
               curz = a(iat,ixyz)
               zmin = min(zmin,curz)
               if (abs(curz-zmin) .lt. toler) imin = iat
 260        continue
            if (iflag .eq. 0) goto 280
            if (lin.eq.2 .and. numfla.gt.1 .and.
     $           sign(one,oldz).ne.sign(one,zmin))
     $           call append('.',ia,length)
            iatflg(imin) = 1
            call append(itoc(ian(imin)),ia,length)
            call append('1',ia,length)
            oldz = zmin
            zmin = big
            goto 250
c
 280        continue
            if (lin.eq.2 .and. numfla.gt.1 .and.
     $           sign(one,oldz).lt.zero)  call append('.',ia,length)
            call append(')',ia,length)
            call append(',',ia,length)
            return
            end
