*deck cdotc
      complex function cdotc (n, cx, incx, cy, incy)
c***begin prologue  cdotc
c***purpose  dot product of two complex vectors using the complex
c            conjugate of the first vector.
c***library   slatec (blas)
c***category  d1a4
c***type      complex (cdotc-c)
c***keywords  blas, inner product, linear algebra, vector
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
c       cy  complex vector with n elements
c     incy  storage spacing between elements of cy
c
c     --output--
c    cdotc  complex result (zero if n .le. 0)
c
c     returns the dot product of complex cx and cy, using conjugate(cx)
c     cdotc = sum for i = 0 to n-1 of conj(cx(lx+i*incx))*cy(ly+i*incy),
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
c***end prologue  cdotc
      complex cx(*),cy(*)
c***first executable statement  cdotc
      cdotc = (0.0,0.0)
      if (n .le. 0) return
      if (incx.eq.incy .and. incx.gt.0) go to 20
c
c     code for unequal or nonpositive increments.
c
      kx = 1
      ky = 1
      if (incx .lt. 0) kx = 1+(1-n)*incx
      if (incy .lt. 0) ky = 1+(1-n)*incy
      do 10 i = 1,n
        cdotc = cdotc + conjg(cx(kx))*cy(ky)
        kx = kx + incx
        ky = ky + incy
   10 continue
      return
c
c     code for equal, positive increments.
c
   20 ns = n*incx
      do 30 i = 1,ns,incx
      cdotc = cdotc + conjg(cx(i))*cy(i)
   30 continue
      return
      end
