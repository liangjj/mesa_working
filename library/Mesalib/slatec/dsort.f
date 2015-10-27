*deck dsort
      subroutine dsort (dx, dy, n, kflag)
c***begin prologue  dsort
c***purpose  sort an array and optionally make the same interchanges in
c            an auxiliary array.  the array may be sorted in increasing
c            or decreasing order.  a slightly modified quicksort
c            algorithm is used.
c***library   slatec
c***category  n6a2b
c***type      double precision (ssort-s, dsort-d, isort-i)
c***keywords  singleton quicksort, sort, sorting
c***author  jones, r. e., (snla)
c           wisniewski, j. a., (snla)
c***description
c
c   dsort sorts array dx and optionally makes the same interchanges in
c   array dy.  the array dx may be sorted in increasing order or
c   decreasing order.  a slightly modified quicksort algorithm is used.
c
c   description of parameters
c      dx - array of values to be sorted   (usually abscissas)
c      dy - array to be (optionally) carried along
c      n  - number of values in array dx to be sorted
c      kflag - control parameter
c            =  2  means sort dx in increasing order and carry dy along.
c            =  1  means sort dx in increasing order (ignoring dy)
c            = -1  means sort dx in decreasing order (ignoring dy)
c            = -2  means sort dx in decreasing order and carry dy along.
c
c***references  r. c. singleton, algorithm 347, an efficient algorithm
c                 for sorting with minimal storage, communications of
c                 the acm, 12, 3 (1969), pp. 185-187.
c***routines called  xermsg
c***revision history  (yymmdd)
c   761101  date written
c   761118  modified to use the singleton quicksort algorithm.  (jaw)
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891009  removed unreferenced statement labels.  (wrb)
c   891024  changed category.  (wrb)
c   891024  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   901012  declared all variables; changed x,y to dx,dy; changed
c           code to parallel ssort. (m. mcclain)
c   920501  reformatted the references section.  (dwl, wrb)
c   920519  clarified error messages.  (dwl)
c   920801  declarations section rebuilt and code restructured to use
c           if-then-else-endif.  (rwc, wrb)
c***end prologue  dsort
c     .. scalar arguments ..
      integer kflag, n
c     .. array arguments ..
      double precision dx(*), dy(*)
c     .. local scalars ..
      double precision r, t, tt, tty, ty
      integer i, ij, j, k, kk, l, m, nn
c     .. local arrays ..
      integer il(21), iu(21)
c     .. external subroutines ..
      external xermsg
c     .. intrinsic functions ..
      intrinsic abs, int
c***first executable statement  dsort
      nn = n
      if (nn .lt. 1) then
         call xermsg ('slatec', 'dsort',
     +      'the number of values to be sorted is not positive.', 1, 1)
         return
      endif
c
      kk = abs(kflag)
      if (kk.ne.1 .and. kk.ne.2) then
         call xermsg ('slatec', 'dsort',
     +      'the sort control parameter, k, is not 2, 1, -1, or -2.', 2,
     +      1)
         return
      endif
c
c     alter array dx to get decreasing order if needed
c
      if (kflag .le. -1) then
         do 10 i=1,nn
            dx(i) = -dx(i)
   10    continue
      endif
c
      if (kk .eq. 2) go to 100
c
c     sort dx only
c
      m = 1
      i = 1
      j = nn
      r = 0.375d0
c
   20 if (i .eq. j) go to 60
      if (r .le. 0.5898437d0) then
         r = r+3.90625d-2
      else
         r = r-0.21875d0
      endif
c
   30 k = i
c
c     select a central element of the array and save it in location t
c
      ij = i + int((j-i)*r)
      t = dx(ij)
c
c     if first element of array is greater than t, interchange with t
c
      if (dx(i) .gt. t) then
         dx(ij) = dx(i)
         dx(i) = t
         t = dx(ij)
      endif
      l = j
c
c     if last element of array is less than than t, interchange with t
c
      if (dx(j) .lt. t) then
         dx(ij) = dx(j)
         dx(j) = t
         t = dx(ij)
c
c        if first element of array is greater than t, interchange with t
c
         if (dx(i) .gt. t) then
            dx(ij) = dx(i)
            dx(i) = t
            t = dx(ij)
         endif
      endif
c
c     find an element in the second half of the array which is smaller
c     than t
c
   40 l = l-1
      if (dx(l) .gt. t) go to 40
c
c     find an element in the first half of the array which is greater
c     than t
c
   50 k = k+1
      if (dx(k) .lt. t) go to 50
c
c     interchange these elements
c
      if (k .le. l) then
         tt = dx(l)
         dx(l) = dx(k)
         dx(k) = tt
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
      t = dx(i+1)
      if (dx(i) .le. t) go to 80
      k = i
c
   90 dx(k+1) = dx(k)
      k = k-1
      if (t .lt. dx(k)) go to 90
      dx(k+1) = t
      go to 80
c
c     sort dx and carry dy along
c
  100 m = 1
      i = 1
      j = nn
      r = 0.375d0
c
  110 if (i .eq. j) go to 150
      if (r .le. 0.5898437d0) then
         r = r+3.90625d-2
      else
         r = r-0.21875d0
      endif
c
  120 k = i
c
c     select a central element of the array and save it in location t
c
      ij = i + int((j-i)*r)
      t = dx(ij)
      ty = dy(ij)
c
c     if first element of array is greater than t, interchange with t
c
      if (dx(i) .gt. t) then
         dx(ij) = dx(i)
         dx(i) = t
         t = dx(ij)
         dy(ij) = dy(i)
         dy(i) = ty
         ty = dy(ij)
      endif
      l = j
c
c     if last element of array is less than t, interchange with t
c
      if (dx(j) .lt. t) then
         dx(ij) = dx(j)
         dx(j) = t
         t = dx(ij)
         dy(ij) = dy(j)
         dy(j) = ty
         ty = dy(ij)
c
c        if first element of array is greater than t, interchange with t
c
         if (dx(i) .gt. t) then
            dx(ij) = dx(i)
            dx(i) = t
            t = dx(ij)
            dy(ij) = dy(i)
            dy(i) = ty
            ty = dy(ij)
         endif
      endif
c
c     find an element in the second half of the array which is smaller
c     than t
c
  130 l = l-1
      if (dx(l) .gt. t) go to 130
c
c     find an element in the first half of the array which is greater
c     than t
c
  140 k = k+1
      if (dx(k) .lt. t) go to 140
c
c     interchange these elements
c
      if (k .le. l) then
         tt = dx(l)
         dx(l) = dx(k)
         dx(k) = tt
         tty = dy(l)
         dy(l) = dy(k)
         dy(k) = tty
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
      t = dx(i+1)
      ty = dy(i+1)
      if (dx(i) .le. t) go to 170
      k = i
c
  180 dx(k+1) = dx(k)
      dy(k+1) = dy(k)
      k = k-1
      if (t .lt. dx(k)) go to 180
      dx(k+1) = t
      dy(k+1) = ty
      go to 170
c
c     clean up
c
  190 if (kflag .le. -1) then
         do 200 i=1,nn
            dx(i) = -dx(i)
  200    continue
      endif
      return
      end
