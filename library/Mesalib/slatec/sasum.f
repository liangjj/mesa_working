*deck sasum
      real function sasum (n, sx, incx)
c***begin prologue  sasum
c***purpose  compute the sum of the magnitudes of the elements of a
c            vector.
c***library   slatec (blas)
c***category  d1a3a
c***type      single precision (sasum-s, dasum-d, scasum-c)
c***keywords  blas, linear algebra, sum of magnitudes of a vector
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
c
c     --output--
c    sasum  single precision result (zero if n .le. 0)
c
c     returns sum of magnitudes of single precision sx.
c     sasum = sum from 0 to n-1 of abs(sx(ix+i*incx)),
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
c***end prologue  sasum
      real sx(*)
      integer i, incx, ix, m, mp1, n
c***first executable statement  sasum
      sasum = 0.0e0
      if (n .le. 0) return
c
      if (incx .eq. 1) goto 20
c
c     code for increment not equal to 1.
c
      ix = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      do 10 i = 1,n
        sasum = sasum + abs(sx(ix))
        ix = ix + incx
   10 continue
      return
c
c     code for increment equal to 1.
c
c     clean-up loop so remaining vector length is a multiple of 6.
c
   20 m = mod(n,6)
      if (m .eq. 0) goto 40
      do 30 i = 1,m
        sasum = sasum + abs(sx(i))
   30 continue
      if (n .lt. 6) return
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        sasum = sasum + abs(sx(i)) + abs(sx(i+1)) + abs(sx(i+2)) +
     1          abs(sx(i+3)) + abs(sx(i+4)) + abs(sx(i+5))
   50 continue
      return
      end
