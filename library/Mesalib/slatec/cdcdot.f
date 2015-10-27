*deck cdcdot
      complex function cdcdot (n, cb, cx, incx, cy, incy)
c***begin prologue  cdcdot
c***purpose  compute the inner product of two vectors with extended
c            precision accumulation.
c***library   slatec (blas)
c***category  d1a4
c***type      complex (sdsdot-s, cdcdot-c)
c***keywords  blas, dot product, inner product, linear algebra, vector
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
c       cb  complex scalar to be added to inner product
c       cx  complex vector with n elements
c     incx  storage spacing between elements of cx
c       cy  complex vector with n elements
c     incy  storage spacing between elements of cy
c
c     --output--
c   cdcdot  complex dot product (cb if n .le. 0)
c
c     returns complex result with dot product accumulated in d.p.
c     cdcdot = cb + sum for i = 0 to n-1 of cx(lx+i*incy)*cy(ly+i*incy)
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
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920310  corrected definition of lx in description.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  cdcdot
      integer n, incx, incy, i, kx, ky
      complex cx(*), cy(*), cb
      double precision dsdotr, dsdoti, dt1, dt2, dt3, dt4
c***first executable statement  cdcdot
      dsdotr = dble(real(cb))
      dsdoti = dble(aimag(cb))
      if (n .le. 0) go to 10
      kx = 1
      ky = 1
      if(incx.lt.0) kx = 1+(1-n)*incx
      if(incy.lt.0) ky = 1+(1-n)*incy
      do 5 i = 1,n
        dt1 = dble(real(cx(kx)))
        dt2 = dble(real(cy(ky)))
        dt3 = dble(aimag(cx(kx)))
        dt4 = dble(aimag(cy(ky)))
        dsdotr = dsdotr+(dt1*dt2)-(dt3*dt4)
        dsdoti = dsdoti+(dt1*dt4)+(dt3*dt2)
        kx = kx+incx
        ky = ky+incy
    5 continue
   10 cdcdot = cmplx(real(dsdotr),real(dsdoti))
      return
      end
