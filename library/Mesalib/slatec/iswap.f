*deck iswap
      subroutine iswap (n, ix, incx, iy, incy)
c***begin prologue  iswap
c***purpose  interchange two vectors.
c***library   slatec (blas)
c***category  d1a5
c***type      integer (sswap-s, dswap-d, cswap-c, iswap-i)
c***keywords  blas, interchange, linear algebra, vector
c***author  vandevender, w. h., (snla)
c***description
c
c                extended b l a s  subprogram
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
c       ix  input vector iy (unchanged if n .le. 0)
c       iy  input vector ix (unchanged if n .le. 0)
c
c     interchange integer ix and integer iy.
c     for i = 0 to n-1, interchange  ix(lx+i*incx) and iy(ly+i*incy),
c     where lx = 1 if incx .ge. 0, else lx = 1+(1-n)*incx, and ly is
c     defined in a similar way using incy.
c
c***references  c. l. lawson, r. j. hanson, d. r. kincaid and f. t.
c                 krogh, basic linear algebra subprograms for fortran
c                 usage, algorithm no. 539, transactions on mathematical
c                 software 5, 3 (september 1979), pp. 308-323.
c***routines called  (none)
c***revision history  (yymmdd)
c   850601  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920310  corrected definition of lx in description.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  iswap
      integer ix(*), iy(*), itemp1, itemp2, itemp3
c***first executable statement  iswap
      if (n .le. 0) return
      if (incx .ne. incy) go to 5
      if (incx-1) 5,20,60
c
c     code for unequal or nonpositive increments.
c
    5 iix = 1
      iiy = 1
      if (incx .lt. 0) iix = (1-n)*incx + 1
      if (incy .lt. 0) iiy = (1-n)*incy + 1
      do 10 i = 1,n
        itemp1 = ix(iix)
        ix(iix) = iy(iiy)
        iy(iiy) = itemp1
        iix = iix + incx
        iiy = iiy + incy
   10 continue
      return
c
c     code for both increments equal to 1.
c
c     clean-up loop so remaining vector length is a multiple of 3.
c
   20 m = mod(n,3)
      if (m .eq. 0) go to 40
      do 30 i = 1,m
        itemp1 = ix(i)
        ix(i) = iy(i)
        iy(i) = itemp1
   30 continue
      if (n .lt. 3) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        itemp1 = ix(i)
        itemp2 = ix(i+1)
        itemp3 = ix(i+2)
        ix(i) = iy(i)
        ix(i+1) = iy(i+1)
        ix(i+2) = iy(i+2)
        iy(i) = itemp1
        iy(i+1) = itemp2
        iy(i+2) = itemp3
   50 continue
      return
c
c     code for equal, positive, non-unit increments.
c
   60 ns = n*incx
      do 70 i = 1,ns,incx
        itemp1 = ix(i)
        ix(i) = iy(i)
        iy(i) = itemp1
   70 continue
      return
      end
