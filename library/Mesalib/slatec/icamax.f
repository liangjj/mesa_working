*deck icamax
      integer function icamax (n, cx, incx)
c***begin prologue  icamax
c***purpose  find the smallest index of the component of a complex
c            vector having the maximum sum of magnitudes of real
c            and imaginary parts.
c***library   slatec (blas)
c***category  d1a2
c***type      complex (isamax-s, idamax-d, icamax-c)
c***keywords  blas, linear algebra, maximum component, vector
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
c   icamax  smallest index (zero if n .le. 0)
c
c     returns the smallest index of the component of cx having the
c     largest sum of magnitudes of real and imaginary parts.
c     icamax = first i, i = 1 to n, to maximize
c     abs(real(cx(ix+(i-1)*incx))) + abs(imag(cx(ix+(i-1)*incx))),
c     where ix = 1 if incx .ge. 0, else ix = 1+(1-n)*incx.
c
c***references  c. l. lawson, r. j. hanson, d. r. kincaid and f. t.
c                 krogh, basic linear algebra subprograms for fortran
c                 usage, algorithm no. 539, transactions on mathematical
c                 software 5, 3 (september 1979), pp. 308-323.
c***routines called  (none)
c***revision history  (yymmdd)
c   791001  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900821  modified to correct problem with a negative increment.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  icamax
      complex cx(*)
      real smax, xmag
      integer i, incx, ix, n
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c***first executable statement  icamax
      icamax = 0
      if (n .le. 0) return
      icamax = 1
      if (n .eq. 1) return
c
      if (incx .eq. 1) goto 20
c
c     code for increment not equal to 1.
c
      ix = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      smax = cabs1(cx(ix))
      ix = ix + incx
      do 10 i = 2,n
        xmag = cabs1(cx(ix))
        if (xmag .gt. smax) then
          icamax = i
          smax = xmag
        endif
        ix = ix + incx
   10 continue
      return
c
c     code for increment equal to 1.
c
   20 smax = cabs1(cx(1))
      do 30 i = 2,n
        xmag = cabs1(cx(i))
        if (xmag .gt. smax) then
          icamax = i
          smax = xmag
        endif
   30 continue
      return
      end
