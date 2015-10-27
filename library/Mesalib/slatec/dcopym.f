*deck dcopym
      subroutine dcopym (n, dx, incx, dy, incy)
c***begin prologue  dcopym
c***purpose  copy the negative of a vector to a vector.
c***library   slatec (blas)
c***category  d1a5
c***type      double precision (scopym-s, dcopym-d)
c***keywords  blas, copy, vector
c***author  kahaner, d. k., (nbs)
c***description
c
c       description of parameters
c           the * flags output variables
c
c       n   number of elements in vector(s)
c      dx   double precision vector with n elements
c    incx   storage spacing between elements of dx
c      dy*  double precision negative copy of dx
c    incy   storage spacing between elements of dy
c
c      ***  note that dy = -dx  ***
c
c     copy negative of d.p. dx to d.p. dy.  for i=0 to n-1,
c     copy  -dx(lx+i*incx) to dy(ly+i*incy), where lx=1 if
c     incx .ge. 0, else lx = 1+(1-n)*incx, and ly is defined
c     in a similar way using incy.
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   801001  date written
c   861211  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920310  corrected definition of lx in description.  (wrb)
c***end prologue  dcopym
      double precision dx(*), dy(*)
c***first executable statement  dcopym
      if (n .le. 0) return
      if (incx .eq. incy) if (incx-1) 5,20,60
c
c     code for unequal or nonpositive increments.
c
   5  ix=1
      iy=1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      if (incy .lt. 0) iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = -dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c     code for both increments equal to 1.
c
c     clean-up loop so remaining vector length is a multiple of 7.
c
   20 m = mod(n,7)
      if (m .eq. 0) go to 40
      do 30 i = 1,m
        dy(i) = -dx(i)
   30 continue
      if (n .lt. 7) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = -dx(i)
        dy(i+1) = -dx(i+1)
        dy(i+2) = -dx(i+2)
        dy(i+3) = -dx(i+3)
        dy(i+4) = -dx(i+4)
        dy(i+5) = -dx(i+5)
        dy(i+6) = -dx(i+6)
   50 continue
      return
c
c     code for equal, positive, non-unit increments.
c
   60 ns = n*incx
      do 70 i = 1,ns,incx
        dy(i) = -dx(i)
   70 continue
      return
      end
