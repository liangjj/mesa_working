*deck scopy
      subroutine scopy (n, sx, incx, sy, incy)
c***begin prologue  scopy
c***purpose  copy a vector.
c***library   slatec (blas)
c***category  d1a5
c***type      single precision (scopy-s, dcopy-d, ccopy-c, icopy-i)
c***keywords  blas, copy, linear algebra, vector
c***author  lawson, c. l., (jpl)
c           hanson, r. j., (snla)
c           kincaid, d. r., (u. of texas)
c           krogh, f. t., (jpl)
c***description
c
c                b l a s  subprogram
c    description of parameters
c
c     --input--
c        n  number of elements in input vector(s)
c       sx  single precision vector with n elements
c     incx  storage spacing between elements of sx
c       sy  single precision vector with n elements
c     incy  storage spacing between elements of sy
c
c     --output--
c       sy  copy of vector sx (unchanged if n .le. 0)
c
c     copy single precision sx to single precision sy.
c     for i = 0 to n-1, copy  sx(lx+i*incx) to sy(ly+i*incy),
c     where lx = 1 if incx .ge. 0, else lx = 1+(1-n)*incx, and ly is
c     defined in a similar way using incy.
c
c***references  c. l. lawson, r. j. hanson, d. r. kincaid and f. t.
c                 krogh, basic linear algebra subprograms for fortran
c                 usage, algorithm no. 539, transactions on mathematical
c                 software 5, 3 (september 1979), pp. 308-323.
c***routines called  (none)
c***revision history  (yymmdd)
c   791001  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920310  corrected definition of lx in description.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  scopy
      real sx(*), sy(*)
c***first executable statement  scopy
      if (n .le. 0) return
      if (incx .eq. incy) if (incx-1) 5,20,60
c
c     code for unequal or nonpositive increments.
c
    5 ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      if (incy .lt. 0) iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
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
        sy(i) = sx(i)
   30 continue
      if (n .lt. 7) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = sx(i)
        sy(i+1) = sx(i+1)
        sy(i+2) = sx(i+2)
        sy(i+3) = sx(i+3)
        sy(i+4) = sx(i+4)
        sy(i+5) = sx(i+5)
        sy(i+6) = sx(i+6)
   50 continue
      return
c
c     code for equal, positive, non-unit increments.
c
   60 ns = n*incx
      do 70 i = 1,ns,incx
        sy(i) = sx(i)
   70 continue
      return
      end