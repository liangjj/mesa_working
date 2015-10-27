*deck zsqrt
      subroutine zsqrt (ar, ai, br, bi)
c***begin prologue  zsqrt
c***subsidiary
c***purpose  subsidiary to zbesh, zbesi, zbesj, zbesk, zbesy, zairy and
c            zbiry
c***library   slatec
c***type      all (zsqrt-a)
c***author  amos, d. e., (snl)
c***description
c
c     double precision complex square root, b=csqrt(a)
c
c***see also  zairy, zbesh, zbesi, zbesj, zbesk, zbesy, zbiry
c***routines called  zabs
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zsqrt
      double precision ar, ai, br, bi, zm, dtheta, dpi, drt
      double precision zabs
      external zabs
      data drt , dpi / 7.071067811865475244008443621d-1,
     1                 3.141592653589793238462643383d+0/
c***first executable statement  zsqrt
      zm = zabs(ar,ai)
      zm = sqrt(zm)
      if (ar.eq.0.0d+0) go to 10
      if (ai.eq.0.0d+0) go to 20
      dtheta = datan(ai/ar)
      if (dtheta.le.0.0d+0) go to 40
      if (ar.lt.0.0d+0) dtheta = dtheta - dpi
      go to 50
   10 if (ai.gt.0.0d+0) go to 60
      if (ai.lt.0.0d+0) go to 70
      br = 0.0d+0
      bi = 0.0d+0
      return
   20 if (ar.gt.0.0d+0) go to 30
      br = 0.0d+0
      bi = sqrt(abs(ar))
      return
   30 br = sqrt(ar)
      bi = 0.0d+0
      return
   40 if (ar.lt.0.0d+0) dtheta = dtheta + dpi
   50 dtheta = dtheta*0.5d+0
      br = zm*cos(dtheta)
      bi = zm*sin(dtheta)
      return
   60 br = zm*drt
      bi = zm*drt
      return
   70 br = zm*drt
      bi = -zm*drt
      return
      end
