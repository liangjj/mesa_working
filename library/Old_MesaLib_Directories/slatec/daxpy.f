*deck daxpy
      subroutine daxpy (n, da, dx, incx, dy, incy)
c***begin prologue  daxpy
c***purpose  compute a constant times a vector plus a vector.
c***library   slatec (blas)
c***category  d1a7
c***type      double precision (saxpy-s, daxpy-d, caxpy-c)
c***keywords  blas, linear algebra, triad, vector
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
c       da  double precision scalar multiplier
c       dx  double precision vector with n elements
c     incx  storage spacing between elements of dx
c       dy  double precision vector with n elements
c     incy  storage spacing between elements of dy
c
c     --output--
c       dy  double precision result (unchanged if n .le. 0)
c
c     overwrite double precision dy with double precision da*dx + dy.
c     for i = 0 to n-1, replace  dy(ly+i*incy) with da*dx(lx+i*incx) +
c       dy(ly+i*incy),
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
c***end prologue  daxpy
      double precision dx(*), dy(*), da
c***first executable statement  daxpy
      if (n.le.0 .or. da.eq.0.0d0) return
      if (incx .eq. incy) if (incx-1) 5,20,60
c
c     code for unequal or nonpositive increments.
c
    5 ix = 1
      iy = 1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      if (incy .lt. 0) iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c     code for both increments equal to 1.
c
c     clean-up loop so remaining vector length is a multiple of 4.
c
   20 m = mod(n,4)
      if (m .eq. 0) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if (n .lt. 4) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i+1) = dy(i+1) + da*dx(i+1)
        dy(i+2) = dy(i+2) + da*dx(i+2)
        dy(i+3) = dy(i+3) + da*dx(i+3)
   50 continue
      return
c
c     code for equal, positive, non-unit increments.
c
   60 ns = n*incx
      do 70 i = 1,ns,incx
        dy(i) = da*dx(i) + dy(i)
   70 continue
      return
      end