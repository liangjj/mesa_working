*deck scopym
      subroutine scopym (n, sx, incx, sy, incy)
c***begin prologue  scopym
c***purpose  copy the negative of a vector to a vector.
c***library   slatec (blas)
c***category  d1a5
c***type      single precision (scopym-s, dcopym-d)
c***keywords  blas, copy, vector
c***author  kahaner, d. k., (nbs)
c***description
c
c       description of parameters
c           the * flags output variables
c
c       n   number of elements in vector(s)
c      sx   real vector with n elements
c    incx   storage spacing between elements of sx
c      sy*  real negative copy of sx
c    incy   storage spacing between elements of sy
c
c      ***  note that sy = -sx  ***
c
c     copy negative of real sx to real sy.  for i=0 to n-1,
c     copy  -sx(lx+i*incx) to sy(ly+i*incy), where lx=1 if
c     incx .ge. 0, else lx = 1+(1-n)*incx, and ly is defined
c     in a similar way using incy.
c
c***references  (none)
c***routines called  (none)
c***revision history  (yymmdd)
c   801001  date written
c   890831  modified array declarations.  (wrb)
c   890831  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920310  corrected definition of lx in description.  (wrb)
c***end prologue  scopym
      real sx(*),sy(*)
c***first executable statement  scopym
      if (n .le. 0) return
      if (incx .eq. incy) if (incx-1) 5,20,60
c
c     code for unequal or nonpositive increments.
c
    5 ix=1
      iy=1
      if (incx .lt. 0) ix = (-n+1)*incx + 1
      if (incy .lt. 0) iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = -sx(ix)
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
        sy(i) = -sx(i)
   30 continue
      if (n .lt. 7) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        sy(i) = -sx(i)
        sy(i+1) = -sx(i+1)
        sy(i+2) = -sx(i+2)
        sy(i+3) = -sx(i+3)
        sy(i+4) = -sx(i+4)
        sy(i+5) = -sx(i+5)
        sy(i+6) = -sx(i+6)
   50 continue
      return
c
c     code for equal, positive, non-unit increments.
c
   60 ns = n*incx
      do 70 i = 1,ns,incx
        sy(i) = -sx(i)
   70 continue
      return
      end
