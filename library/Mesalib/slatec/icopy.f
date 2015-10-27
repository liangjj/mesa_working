*deck icopy
      subroutine icopy (n, ix, incx, iy, incy)
c***begin prologue  icopy
c***purpose  copy a vector.
c***library   slatec (blas)
c***category  d1a5
c***type      integer (icopy-s, dcopy-d, ccopy-c, icopy-i)
c***keywords  blas, copy, linear algebra, vector
c***author  boland, w. robert, (lanl)
c           clemens, reginald, (plk)
c***description
c
c                b l a s  subprogram
c    description of parameters
c
c     --input--
c        n  number of elements in input vector(s)
c       ix  integer vector with n elements
c     incx  storage spacing between elements of ix
c       iy  integer vector with n elements
c     incy  storage spacing between elements of iy
c
c     --output--
c       iy  copy of vector ix (unchanged if n .le. 0)
c
c     copy integer ix to integer iy.
c     for i = 0 to n-1, copy  ix(lx+i*incx) to iy(ly+i*incy),
c     where lx = 1 if incx .ge. 0, else lx = 1+(1-n)*incx, and ly is
c     defined in a similar way using incy.
c
c***references  c. l. lawson, r. j. hanson, d. r. kincaid and f. t.
c                 krogh, basic linear algebra subprograms for fortran
c                 usage, algorithm no. 539, transactions on mathematical
c                 software 5, 3 (september 1979), pp. 308-323.
c***routines called  (none)
c***revision history  (yymmdd)
c   930201  date written
c***end prologue  icopy
      integer ix(*), iy(*)
c***first executable statement  icopy
      if (n .le. 0) return
      if (incx .eq. incy) if (incx-1) 5,20,60
c
c     code for unequal or nonpositive increments.
c
    5 iix = 1
      iiy = 1
      if (incx .lt. 0) iix = (-n+1)*incx + 1
      if (incy .lt. 0) iiy = (-n+1)*incy + 1
      do 10 i = 1,n
        iy(iiy) = ix(iix)
        iix = iix + incx
        iiy = iiy + incy
   10 continue
      return
c
c     code for both increments equal to 1.
c
c     clean-up loop so remaining vector length is a multiple of 7.
c
   20 m = mod(n,7)
      if (m .eq. 0) go to 40
      do 30 i = 1,m
        iy(i) = ix(i)
   30 continue
      if (n .lt. 7) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        iy(i) = ix(i)
        iy(i+1) = ix(i+1)
        iy(i+2) = ix(i+2)
        iy(i+3) = ix(i+3)
        iy(i+4) = ix(i+4)
        iy(i+5) = ix(i+5)
        iy(i+6) = ix(i+6)
   50 continue
      return
c
c     code for equal, positive, non-unit increments.
c
   60 ns = n*incx
      do 70 i = 1,ns,incx
        iy(i) = ix(i)
   70 continue
      return
      end
