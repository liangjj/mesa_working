      function cdotc(n,cx,incx,cy,incy)
c***begin prologue  cdotc
c***date written   791001   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  d1a4
c***keywords  blas,complex,inner product,linear algebra,vector
c***author  lawson, c. l., (jpl)
c           hanson, r. j., (snla)
c           kincaid, d. r., (u. of texas)
c           krogh, f. t., (jpl)
c***purpose  dot product of complex vectors, uses complx conjugate of
c            first vector
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
c     returns the dot product for complex cx and cy, uses conjugate(cx)
c     cdotc = sum for i = 0 to n-1 of conj(cx(lx+i*incx))*cy(ly+i*incy)
c     where lx = 1 if incx .ge. 0, else lx = (-incx)*n, and ly is
c     defined in a similar way using incy.
c***references  lawson c.l., hanson r.j., kincaid d.r., krogh f.t.,
c                 *basic linear algebra subprograms for fortran usage*,
c                 algorithm no. 539, transactions on mathematical
c                 software, volume 5, number 3, september 1979, 308-323
c***routines called  (none)
c***end prologue  cdotc
c
      implicit real *8 (a-h,o-z)
      complex*16 cx(*), cy(*), cdotc, zdotc
c***first executable statement  cdotc
      cdotc=zdotc(n,cx,incx,cy,incy)
      return
      end
