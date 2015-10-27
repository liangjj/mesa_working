*deck ssort
      subroutine ssort (x, y, n, kflag)
c***begin prologue  ssort
c***purpose  sort an array and optionally make the same interchanges in
c            an auxiliary array.  the array may be sorted in increasing
c            or decreasing order.  a slightly modified quicksort
c            algorithm is used.
c***library   slatec
c***category  n6a2b
c***type      single precision (ssort-s, dsort-d, isort-i)
c***keywords  singleton quicksort, sort, sorting
c***author  jones, r. e., (snla)
c           wisniewski, j. a., (snla)
c***description
c
c   ssort sorts array x and optionally makes the same interchanges in
c   array y.  the array x may be sorted in increasing order or
c   decreasing order.  a slightly modified quicksort algorithm is used.
c
c   description of parameters
c      x - array of values to be sorted   (usually abscissas)
c      y - array to be (optionally) carried along
c      n - number of values in array x to be sorted
c      kflag - control parameter
c            =  2  means sort x in increasing order and carry y along.
c            =  1  means sort x in increasing order (ignoring y)
c            = -1  means sort x in decreasing order (ignoring y)
c            = -2  means sort x in decreasing order and carry y along.
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
c   901012  declared all variables; changed x,y to sx,sy. (m. mcclain)
c   920501  reformatted the references section.  (dwl, wrb)
c   920519  clarified error messages.  (dwl)
c   920801  declarations section rebuilt and code restructured to use
c           if-then-else-endif.  (rwc, wrb)
c***end prologue  ssort
c     .. scalar arguments ..
      integer kflag, n
c     .. array arguments ..
      real x(*), y(*)
c     .. local scalars ..
      real r, t, tt, tty, ty
      integer i, ij, j, k, kk, l, m, nn
c     .. local arrays ..
      integer il(21), iu(21)
c     .. external subroutines ..
      external xermsg
c     .. intrinsic functions ..
      intrinsic abs, int
c***first executable statement  ssort
      nn = n
      if (nn .lt. 1) then
         call xermsg ('slatec', 'ssort',
     +      'the number of values to be sorted is not positive.', 1, 1)
         return
      endif
c
      kk = abs(kflag)
      if (kk.ne.1 .and. kk.ne.2) then
         call xermsg ('slatec', 'ssort',
     +      'the sort control parameter, k, is not 2, 1, -1, or -2.', 2,
     +      1)
         return
      endif
c
c     alter array x to get decreasing order if needed
c
      if (kflag .le. -1) then
         do 10 i=1,nn
            x(i) = -x(i)
   10    continue
      endif
c
      if (kk .eq. 2) go to 100
c
c     sort x only
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
      t = x(ij)
c
c     if first element of array is greater than t, interchange with t
c
      if (x(i) .gt. t) then
         x(ij) = x(i)
         x(i) = t
         t = x(ij)
      endif
      l = j
c
c     if last element of array is less than than t, interchange with t
c
      if (x(j) .lt. t) then
         x(ij) = x(j)
         x(j) = t
         t = x(ij)
c
c        if first element of array is greater than t, interchange with t
c
         if (x(i) .gt. t) then
            x(ij) = x(i)
            x(i) = t
            t = x(ij)
         endif
      endif
c
c     find an element in the second half of the array which is smaller
c     than t
c
   40 l = l-1
      if (x(l) .gt. t) go to 40
c
c     find an element in the first half of the array which is greater
c     than t
c
   50 k = k+1
      if (x(k) .lt. t) go to 50
c
c     interchange these elements
c
      if (k .le. l) then
         tt = x(l)
         x(l) = x(k)
         x(k) = tt
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
      t = x(i+1)
      if (x(i) .le. t) go to 80
      k = i
c
   90 x(k+1) = x(k)
      k = k-1
      if (t .lt. x(k)) go to 90
      x(k+1) = t
      go to 80
c
c     sort x and carry y along
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
      t = x(ij)
      ty = y(ij)
c
c     if first element of array is greater than t, interchange with t
c
      if (x(i) .gt. t) then
         x(ij) = x(i)
         x(i) = t
         t = x(ij)
         y(ij) = y(i)
         y(i) = ty
         ty = y(ij)
      endif
      l = j
c
c     if last element of array is less than t, interchange with t
c
      if (x(j) .lt. t) then
         x(ij) = x(j)
         x(j) = t
         t = x(ij)
         y(ij) = y(j)
         y(j) = ty
         ty = y(ij)
c
c        if first element of array is greater than t, interchange with t
c
         if (x(i) .gt. t) then
            x(ij) = x(i)
            x(i) = t
            t = x(ij)
            y(ij) = y(i)
            y(i) = ty
            ty = y(ij)
         endif
      endif
c
c     find an element in the second half of the array which is smaller
c     than t
c
  130 l = l-1
      if (x(l) .gt. t) go to 130
c
c     find an element in the first half of the array which is greater
c     than t
c
  140 k = k+1
      if (x(k) .lt. t) go to 140
c
c     interchange these elements
c
      if (k .le. l) then
         tt = x(l)
         x(l) = x(k)
         x(k) = tt
         tty = y(l)
         y(l) = y(k)
         y(k) = tty
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
      t = x(i+1)
      ty = y(i+1)
      if (x(i) .le. t) go to 170
      k = i
c
  180 x(k+1) = x(k)
      y(k+1) = y(k)
      k = k-1
      if (t .lt. x(k)) go to 180
      x(k+1) = t
      y(k+1) = ty
      go to 170
c
c     clean up
c
  190 if (kflag .le. -1) then
         do 200 i=1,nn
            x(i) = -x(i)
  200    continue
      endif
      return
      end
