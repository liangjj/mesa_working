*deck csroot
      subroutine csroot (xr, xi, yr, yi)
c***begin prologue  csroot
c***subsidiary
c***purpose  compute the complex square root of a complex number.
c***library   slatec
c***type      single precision (csroot-s)
c***author  (unknown)
c***description
c
c     (yr,yi) = complex sqrt(xr,xi)
c
c***see also  eisdoc
c***routines called  pythag
c***revision history  (yymmdd)
c   811101  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900402  added type section.  (wrb)
c***end prologue  csroot
      real xr,xi,yr,yi,s,tr,ti,pythag
c
c     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)
c***first executable statement  csroot
      tr = xr
      ti = xi
      s = sqrt(0.5e0*(pythag(tr,ti) + abs(tr)))
      if (tr .ge. 0.0e0) yr = s
      if (ti .lt. 0.0e0) s = -s
      if (tr .le. 0.0e0) yi = s
      if (tr .lt. 0.0e0) yr = 0.5e0*(ti/yi)
      if (tr .gt. 0.0e0) yi = 0.5e0*(ti/yr)
      return
      end
