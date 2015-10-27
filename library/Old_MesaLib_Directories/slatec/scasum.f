*deck scasum
      function scasum (n, cx, incx)
c***begin prologue  scasum
c***purpose  compute the sum of the magnitudes of the real and
c            imaginary elements of a complex vector.
c***library   slatec (blas)
c***category  d1a3a
c***type      complex (sasum-s, dasum-d, scasum-c)
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
c       cx  complex vector with n elements
c     incx  storage spacing between elements of cx
c
c     --output--
c   scasum  single precision result (zero if n .le. 0)
c
c     returns sums of magnitudes of real and imaginary parts of
c     components of cx.  note that this is not the l1 norm of cx.
c     casum = sum from 0 to n-1 of abs(real(cx(ix+i*incx))) +
c             abs(imag(cx(ix+i*incx))),
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
c***end prologue  scasum
      complex cx(*)
      integer i, incx, ix, n
c***first executable statement  scasum
      scasum = 0.0e0
      if (n .le. 0) return
c
      if (incx .eq. 1) goto 20
c
c     code for increment not equal to 1.
c
      ix = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      do 10 i = 1,n
        scasum = scasum + abs(real(cx(ix))) + abs(aimag(cx(ix)))
        ix = ix + incx
   10 continue
      return
c
c     code for increment equal to 1.
c
   20 do 30 i = 1,n
        scasum = scasum + abs(real(cx(i))) + abs(aimag(cx(i)))
   30 continue
      return
      end
