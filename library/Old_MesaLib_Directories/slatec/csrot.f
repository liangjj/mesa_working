*deck csrot
      subroutine csrot (n, cx, incx, cy, incy, c, s)
c***begin prologue  csrot
c***purpose  apply a plane givens rotation.
c***library   slatec (blas)
c***category  d1b10
c***type      complex (srot-s, drot-d, csrot-c)
c***keywords  blas, givens rotation, givens transformation,
c             linear algebra, plane rotation, vector
c***author  dongarra, j., (anl)
c***description
c
c     csrot applies the complex givens rotation
c
c          (x)   ( c s)(x)
c          (y) = (-s c)(y)
c
c     n times where for i = 0,...,n-1
c
c          x = cx(lx+i*incx)
c          y = cy(ly+i*incy),
c
c     where lx = 1 if incx .ge. 0, else lx = 1+(1-n)*incx, and ly is
c     defined in a similar way using incy.
c
c     argument description
c
c        n      (integer)  number of elements in each vector
c
c        cx     (complex array)  beginning of one vector
c
c        incx   (integer)  memory spacing of successive elements
c               of vector cx
c
c        cy     (complex array)  beginning of the other vector
c
c        incy   (integer)  memory spacing of successive elements
c               of vector cy
c
c        c      (real)  cosine term of the rotation
c
c        s      (real)  sine term of the rotation.
c
c***references  j. j. dongarra, j. r. bunch, c. b. moler, and g. w.
c                 stewart, linpack users' guide, siam, 1979.
c***routines called  (none)
c***revision history  (yymmdd)
c   810223  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920310  corrected definition of lx in description.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  csrot
      complex cx(*), cy(*), ctemp
      real c, s
      integer i, incx, incy, ix, iy, n
c***first executable statement  csrot
      if (n .le. 0) return
      if (incx.eq.1 .and. incy.eq.1)go to 20
c
c     code for unequal increments or equal increments not equal to 1.
c
      ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      if (incy .lt. 0) iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ctemp = c*cx(ix) + s*cy(iy)
        cy(iy) = c*cy(iy) - s*cx(ix)
        cx(ix) = ctemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c     code for both increments equal to 1.
c
   20 do 30 i = 1,n
        ctemp = c*cx(i) + s*cy(i)
        cy(i) = c*cy(i) - s*cx(i)
        cx(i) = ctemp
   30 continue
      return
      end
