      subroutine csroot(xr,xi,yr,yi)
c***begin prologue  csroot
c***refer to  eisdoc
c
c     (yr,yi) = complex sqrt(xr,xi)
c***routines called  pythag
c***end prologue  csroot
      implicit real *8 (a-h,o-z)
c
c     branch chosen so that yr .ge. 0.0 and sign(yi) .eq. sign(xi)
c***first executable statement  csroot
      tr = xr
      ti = xi
      s = sqrt(0.5d0*(pythag(tr,ti) + abs(tr)))
      if (tr .ge. 0.0d0) yr = s
      if (ti .lt. 0.0d0) s = -s
      if (tr .le. 0.0d0) yi = s
      if (tr .lt. 0.0d0) yr = 0.5d0*(ti/yi)
      if (tr .gt. 0.0d0) yi = 0.5d0*(ti/yr)
      return
      end
