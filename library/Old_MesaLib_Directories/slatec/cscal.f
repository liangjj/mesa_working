*deck cscal
      subroutine cscal (n, ca, cx, incx)
c***begin prologue  cscal
c***purpose  multiply a vector by a constant.
c***library   slatec (blas)
c***category  d1a6
c***type      complex (sscal-s, dscal-d, cscal-c)
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
c       ca  complex scale factor
c       cx  complex vector with n elements
c     incx  storage spacing between elements of cx
c
c     --output--
c       cx  complex result (unchanged if n .le. 0)
c
c     replace complex cx by complex ca*cx.
c     for i = 0 to n-1, replace cx(ix+i*incx) with ca*cx(ix+i*incx),
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
c***end prologue  cscal
      complex ca, cx(*)
      integer i, incx, ix, n
c***first executable statement  cscal
      if (n .le. 0) return
c
      if (incx .eq. 1) goto 20
c
c     code for increment not equal to 1.
c
      ix = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      do 10 i = 1,n
        cx(ix) = ca*cx(ix)
        ix = ix + incx
   10 continue
      return
c
c     code for increment equal to 1.
c
   20 do 30 i = 1,n
        cx(i) = ca*cx(i)
   30 continue
      return
      end
