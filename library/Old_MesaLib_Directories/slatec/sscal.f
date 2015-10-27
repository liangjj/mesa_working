*deck sscal
      subroutine sscal (n, sa, sx, incx)
c***begin prologue  sscal
c***purpose  multiply a vector by a constant.
c***library   slatec (blas)
c***category  d1a6
c***type      single precision (sscal-s, dscal-d, cscal-c)
c***keywords  blas, linear algebra, scale, vector
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
c       sa  single precision scale factor
c       sx  single precision vector with n elements
c     incx  storage spacing between elements of sx
c
c     --output--
c       sx  single precision result (unchanged if n .le. 0)
c
c     replace single precision sx by single precision sa*sx.
c     for i = 0 to n-1, replace sx(ix+i*incx) with  sa * sx(ix+i*incx),
c     where ix = 1 if incx .ge. 0, else ix = 1+(1-n)*incx.
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
c   900821  modified to correct problem with a negative increment.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  sscal
      real sa, sx(*)
      integer i, incx, ix, m, mp1, n
c***first executable statement  sscal
      if (n .le. 0) return
      if (incx .eq. 1) goto 20
c
c     code for increment not equal to 1.
c
      ix = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      do 10 i = 1,n
        sx(ix) = sa*sx(ix)
        ix = ix + incx
   10 continue
      return
c
c     code for increment equal to 1.
c
c     clean-up loop so remaining vector length is a multiple of 5.
c
   20 m = mod(n,5)
      if (m .eq. 0) goto 40
      do 30 i = 1,m
        sx(i) = sa*sx(i)
   30 continue
      if (n .lt. 5) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sx(i) = sa*sx(i)
        sx(i+1) = sa*sx(i+1)
        sx(i+2) = sa*sx(i+2)
        sx(i+3) = sa*sx(i+3)
        sx(i+4) = sa*sx(i+4)
   50 continue
      return
      end
