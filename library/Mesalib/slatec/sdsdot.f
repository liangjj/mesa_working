*deck sdsdot
      real function sdsdot (n, sb, sx, incx, sy, incy)
c***begin prologue  sdsdot
c***purpose  compute the inner product of two vectors with extended
c            precision accumulation.
c***library   slatec (blas)
c***category  d1a4
c***type      single precision (sdsdot-s, cdcdot-c)
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
c       sb  single precision scalar to be added to inner product
c       sx  single precision vector with n elements
c     incx  storage spacing between elements of sx
c       sy  single precision vector with n elements
c     incy  storage spacing between elements of sy
c
c     --output--
c   sdsdot  single precision dot product (sb if n .le. 0)
c
c     returns s.p. result with dot product accumulated in d.p.
c     sdsdot = sb + sum for i = 0 to n-1 of sx(lx+i*incx)*sy(ly+i*incy),
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
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920310  corrected definition of lx in description.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  sdsdot
      real sx(*), sy(*), sb
      double precision dsdot
c***first executable statement  sdsdot
      dsdot = sb
      if (n .le. 0) go to 30
      if (incx.eq.incy .and. incx.gt.0) go to 40
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
   30 sdsdot = dsdot
      return
c
c     code for equal and positive increments.
c
   40 ns = n*incx
      do 50 i = 1,ns,incx
        dsdot = dsdot + dble(sx(i))*dble(sy(i))
   50 continue
      sdsdot = dsdot
      return
      end
