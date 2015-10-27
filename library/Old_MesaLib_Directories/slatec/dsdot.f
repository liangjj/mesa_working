*deck dsdot
      double precision function dsdot (n, sx, incx, sy, incy)
c***begin prologue  dsdot
c***purpose  compute the inner product of two vectors with extended
c            precision accumulation and result.
c***library   slatec (blas)
c***category  d1a4
c***type      double precision (dsdot-d, dcdot-c)
c***keywords  blas, complex vectors, dot product, inner product,
c             linear algebra, vector
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
c    dsdot  double precision dot product (zero if n.le.0)
c
c     returns d.p. dot product accumulated in d.p., for s.p. sx and sy
c     dsdot = sum for i = 0 to n-1 of  sx(lx+i*incx) * sy(ly+i*incy),
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
c***end prologue  dsdot
      real sx(*),sy(*)
c***first executable statement  dsdot
      dsdot = 0.0d0
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
        dsdot = dsdot + dble(sx(kx))*dble(sy(ky))
        kx = kx + incx
        ky = ky + incy
   10 continue
      return
c
c     code for equal, positive, non-unit increments.
c
   20 ns = n*incx
      do 30 i = 1,ns,incx
        dsdot = dsdot + dble(sx(i))*dble(sy(i))
   30 continue
      return
      end
