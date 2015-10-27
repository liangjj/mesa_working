*deck zlog
      subroutine zlog (ar, ai, br, bi, ierr)
c***begin prologue  zlog
c***subsidiary
c***purpose  subsidiary to zbesh, zbesi, zbesj, zbesk, zbesy, zairy and
c            zbiry
c***library   slatec
c***type      all (zlog-a)
c***author  amos, d. e., (snl)
c***description
c
c     double precision complex logarithm b=clog(a)
c     ierr=0,normal return      ierr=1, z=cmplx(0.0,0.0)
c***see also  zairy, zbesh, zbesi, zbesj, zbesk, zbesy, zbiry
c***routines called  zabs
c***revision history  (yymmdd)
c   830501  date written
c   910415  prologue converted to version 4.0 format.  (bab)
c***end prologue  zlog
      double precision ar, ai, br, bi, zm, dtheta, dpi, dhpi
      double precision zabs
      integer ierr
      external zabs
      data dpi , dhpi  / 3.141592653589793238462643383d+0,
     1                   1.570796326794896619231321696d+0/
c***first executable statement  zlog
      ierr=0
      if (ar.eq.0.0d+0) go to 10
      if (ai.eq.0.0d+0) go to 20
      dtheta = datan(ai/ar)
      if (dtheta.le.0.0d+0) go to 40
      if (ar.lt.0.0d+0) dtheta = dtheta - dpi
      go to 50
   10 if (ai.eq.0.0d+0) go to 60
      bi = dhpi
      br = log(abs(ai))
      if (ai.lt.0.0d+0) bi = -bi
      return
   20 if (ar.gt.0.0d+0) go to 30
      br = log(abs(ar))
      bi = dpi
      return
   30 br = log(ar)
      bi = 0.0d+0
      return
   40 if (ar.lt.0.0d+0) dtheta = dtheta + dpi
   50 zm = zabs(ar,ai)
      br = log(zm)
      bi = dtheta
      return
   60 continue
      ierr=1
      return
      end
