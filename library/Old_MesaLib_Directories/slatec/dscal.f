*deck dscal
      subroutine dscal (n, da, dx, incx)
c***begin prologue  dscal
c***purpose  multiply a vector by a constant.
c***library   slatec (blas)
c***category  d1a6
c***type      double precision (sscal-s, dscal-d, cscal-c)
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
c       da  double precision scale factor
c       dx  double precision vector with n elements
c     incx  storage spacing between elements of dx
c
c     --output--
c       dx  double precision result (unchanged if n.le.0)
c
c     replace double precision dx by double precision da*dx.
c     for i = 0 to n-1, replace dx(ix+i*incx) with  da * dx(ix+i*incx),
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
c***end prologue  dscal
      double precision da, dx(*)
      integer i, incx, ix, m, mp1, n
c***first executable statement  dscal
      if (n .le. 0) return
      if (incx .eq. 1) goto 20
c
c     code for increment not equal to 1.
c
      ix = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      do 10 i = 1,n
        dx(ix) = da*dx(ix)
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
        dx(i) = da*dx(i)
   30 continue
      if (n .lt. 5) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i+1) = da*dx(i+1)
        dx(i+2) = da*dx(i+2)
        dx(i+3) = da*dx(i+3)
        dx(i+4) = da*dx(i+4)
   50 continue
      return
      end
