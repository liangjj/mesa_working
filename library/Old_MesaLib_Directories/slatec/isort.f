*deck isort
      subroutine isort (ix, iy, n, kflag)
c***begin prologue  isort
c***purpose  sort an array and optionally make the same interchanges in
c            an auxiliary array.  the array may be sorted in increasing
c            or decreasing order.  a slightly modified quicksort
c            algorithm is used.
c***library   slatec
c***category  n6a2a
c***type      integer (ssort-s, dsort-d, isort-i)
c***keywords  singleton quicksort, sort, sorting
c***author  jones, r. e., (snla)
c           kahaner, d. k., (nbs)
c           wisniewski, j. a., (snla)
c***description
c
c   isort sorts array ix and optionally makes the same interchanges in
c   array iy.  the array ix may be sorted in increasing order or
c   decreasing order.  a slightly modified quicksort algorithm is used.
c
c   description of parameters
c      ix - integer array of values to be sorted
c      iy - integer array to be (optionally) carried along
c      n  - number of values in integer array ix to be sorted
c      kflag - control parameter
c            =  2  means sort ix in increasing order and carry iy along.
c            =  1  means sort ix in increasing order (ignoring iy)
c            = -1  means sort ix in decreasing order (ignoring iy)
c            = -2  means sort ix in decreasing order and carry iy along.
c
c***references  r. c. singleton, algorithm 347, an efficient algorithm
c                 for sorting with minimal storage, communications of
c                 the acm, 12, 3 (1969), pp. 185-187.
c***routines called  xermsg
c***revision history  (yymmdd)
c   761118  date written
c   810801  modified by david k. kahaner.
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891009  removed unreferenced statement labels.  (wrb)
c   891009  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   901012  declared all variables; changed x,y to ix,iy. (m. mcclain)
c   920501  reformatted the references section.  (dwl, wrb)
c   920519  clarified error messages.  (dwl)
c   920801  declarations section rebuilt and code restructured to use
c           if-then-else-endif.  (rwc, wrb)
c***end prologue  isort
c     .. scalar arguments ..
      integer kflag, n
c     .. array arguments ..
      integer ix(*), iy(*)
c     .. local scalars ..
      real r
      integer i, ij, j, k, kk, l, m, nn, t, tt, tty, ty
c     .. local arrays ..
      integer il(21), iu(21)
c     .. external subroutines ..
      external xermsg
c     .. intrinsic functions ..
      intrinsic abs, int
c***first executable statement  isort
      nn = n
      if (nn .lt. 1) then
         call xermsg ('slatec', 'isort',
     +      'the number of values to be sorted is not positive.', 1, 1)
         return
      endif
c
      kk = abs(kflag)
      if (kk.ne.1 .and. kk.ne.2) then
         call xermsg ('slatec', 'isort',
     +      'the sort control parameter, k, is not 2, 1, -1, or -2.', 2,
     +      1)
         return
      endif
c
c     alter array ix to get decreasing order if needed
c
      if (kflag .le. -1) then
         do 10 i=1,nn
            ix(i) = -ix(i)
   10    continue
      endif
c
      if (kk .eq. 2) go to 100
c
c     sort ix only
c
      m = 1
      i = 1
      j = nn
      r = 0.375e0
c
   20 if (i .eq. j) go to 60
      if (r .le. 0.5898437e0) then
         r = r+3.90625e-2
      else
         r = r-0.21875e0
      endif
c
   30 k = i
c
c     select a central element of the array and save it in location t
c
      ij = i + int((j-i)*r)
      t = ix(ij)
c
c     if first element of array is greater than t, interchange with t
c
      if (ix(i) .gt. t) then
         ix(ij) = ix(i)
         ix(i) = t
         t = ix(ij)
      endif
      l = j
c
c     if last element of array is less than than t, interchange with t
c
      if (ix(j) .lt. t) then
         ix(ij) = ix(j)
         ix(j) = t
         t = ix(ij)
c
c        if first element of array is greater than t, interchange with t
c
         if (ix(i) .gt. t) then
            ix(ij) = ix(i)
            ix(i) = t
            t = ix(ij)
         endif
      endif
c
c     find an element in the second half of the array which is smaller
c     than t
c
   40 l = l-1
      if (ix(l) .gt. t) go to 40
c
c     find an element in the first half of the array which is greater
c     than t
c
   50 k = k+1
      if (ix(k) .lt. t) go to 50
c
c     interchange these elements
c
      if (k .le. l) then
         tt = ix(l)
         ix(l) = ix(k)
         ix(k) = tt
         go to 40
      endif
c
c     save upper and lower subscripts of the array yet to be sorted
c
      if (l-i .gt. j-k) then
         il(m) = i
         iu(m) = l
         i = k
         m = m+1
      else
         il(m) = k
         iu(m) = j
         j = l
         m = m+1
      endif
      go to 70
c
c     begin again on another portion of the unsorted array
c
   60 m = m-1
      if (m .eq. 0) go to 190
      i = il(m)
      j = iu(m)
c
   70 if (j-i .ge. 1) go to 30
      if (i .eq. 1) go to 20
      i = i-1
c
   80 i = i+1
      if (i .eq. j) go to 60
      t = ix(i+1)
      if (ix(i) .le. t) go to 80
      k = i
c
   90 ix(k+1) = ix(k)
      k = k-1
      if (t .lt. ix(k)) go to 90
      ix(k+1) = t
      go to 80
c
c     sort ix and carry iy along
c
  100 m = 1
      i = 1
      j = nn
      r = 0.375e0
c
  110 if (i .eq. j) go to 150
      if (r .le. 0.5898437e0) then
         r = r+3.90625e-2
      else
         r = r-0.21875e0
      endif
c
  120 k = i
c
c     select a central element of the array and save it in location t
c
      ij = i + int((j-i)*r)
      t = ix(ij)
      ty = iy(ij)
c
c     if first element of array is greater than t, interchange with t
c
      if (ix(i) .gt. t) then
         ix(ij) = ix(i)
         ix(i) = t
         t = ix(ij)
         iy(ij) = iy(i)
         iy(i) = ty
         ty = iy(ij)
      endif
      l = j
c
c     if last element of array is less than t, interchange with t
c
      if (ix(j) .lt. t) then
         ix(ij) = ix(j)
         ix(j) = t
         t = ix(ij)
         iy(ij) = iy(j)
         iy(j) = ty
         ty = iy(ij)
c
c        if first element of array is greater than t, interchange with t
c
         if (ix(i) .gt. t) then
            ix(ij) = ix(i)
            ix(i) = t
            t = ix(ij)
            iy(ij) = iy(i)
            iy(i) = ty
            ty = iy(ij)
         endif
      endif
c
c     find an element in the second half of the array which is smaller
c     than t
c
  130 l = l-1
      if (ix(l) .gt. t) go to 130
c
c     find an element in the first half of the array which is greater
c     than t
c
  140 k = k+1
      if (ix(k) .lt. t) go to 140
c
c     interchange these elements
c
      if (k .le. l) then
         tt = ix(l)
         ix(l) = ix(k)
         ix(k) = tt
         tty = iy(l)
         iy(l) = iy(k)
         iy(k) = tty
         go to 130
      endif
c
c     save upper and lower subscripts of the array yet to be sorted
c
      if (l-i .gt. j-k) then
         il(m) = i
         iu(m) = l
         i = k
         m = m+1
      else
         il(m) = k
         iu(m) = j
         j = l
         m = m+1
      endif
      go to 160
c
c     begin again on another portion of the unsorted array
c
  150 m = m-1
      if (m .eq. 0) go to 190
      i = il(m)
      j = iu(m)
c
  160 if (j-i .ge. 1) go to 120
      if (i .eq. 1) go to 110
      i = i-1
c
  170 i = i+1
      if (i .eq. j) go to 150
      t = ix(i+1)
      ty = iy(i+1)
      if (ix(i) .le. t) go to 170
      k = i
c
  180 ix(k+1) = ix(k)
      iy(k+1) = iy(k)
      k = k-1
      if (t .lt. ix(k)) go to 180
      ix(k+1) = t
      iy(k+1) = ty
      go to 170
c
c     clean up
c
  190 if (kflag .le. -1) then
         do 200 i=1,nn
            ix(i) = -ix(i)
  200    continue
      endif
      return
      end
