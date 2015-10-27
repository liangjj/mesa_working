*deck dxred
      subroutine dxred (x, ix, ierror)
c***begin prologue  dxred
c***purpose  to provide double-precision floating-point arithmetic
c            with an extended exponent range.
c***library   slatec
c***category  a3d
c***type      double precision (xred-s, dxred-d)
c***keywords  extended-range double-precision arithmetic
c***author  lozier, daniel w., (national bureau of standards)
c           smith, john m., (nbs and george mason university)
c***description
c     double precision x
c     integer ix
c
c                  if
c                  radix**(-2l) .le. (abs(x),ix) .le. radix**(2l)
c                  then dxred transforms (x,ix) so that ix=0.
c                  if (x,ix) is outside the above range,
c                  then dxred takes no action.
c                  this subroutine is useful if the
c                  results of extended-range calculations
c                  are to be used in subsequent ordinary
c                  double-precision calculations.
c
c***see also  dxset
c***references  (none)
c***routines called  (none)
c***common blocks    dxblk2
c***revision history  (yymmdd)
c   820712  date written
c   881020  revised to meet slatec cml recommendations.  (dwl and jms)
c   901019  revisions to prologue.  (dwl and wrb)
c   901106  changed all specific intrinsics to generic.  (wrb)
c           corrected order of sections in prologue and added type
c           section.  (wrb)
c   920127  revised purpose section of prologue.  (dwl)
c***end prologue  dxred
      double precision x
      integer ix
      double precision radix, radixl, rad2l, dlg10r, xa
      integer l, l2, kmax
      common /dxblk2/ radix, radixl, rad2l, dlg10r, l, l2, kmax
      save /dxblk2/
c
c***first executable statement  dxred
      ierror=0
      if (x.eq.0.0d0) go to 90
      xa = abs(x)
      if (ix.eq.0) go to 70
      ixa = abs(ix)
      ixa1 = ixa/l2
      ixa2 = mod(ixa,l2)
      if (ix.gt.0) go to 40
   10 continue
      if (xa.gt.1.0d0) go to 20
      xa = xa*rad2l
      ixa1 = ixa1 + 1
      go to 10
   20 xa = xa/radix**ixa2
      if (ixa1.eq.0) go to 70
      do 30 i=1,ixa1
        if (xa.lt.1.0d0) go to 100
        xa = xa/rad2l
   30 continue
      go to 70
c
   40 continue
      if (xa.lt.1.0d0) go to 50
      xa = xa/rad2l
      ixa1 = ixa1 + 1
      go to 40
   50 xa = xa*radix**ixa2
      if (ixa1.eq.0) go to 70
      do 60 i=1,ixa1
        if (xa.gt.1.0d0) go to 100
        xa = xa*rad2l
   60 continue
   70 if (xa.gt.rad2l) go to 100
      if (xa.gt.1.0d0) go to 80
      if (rad2l*xa.lt.1.0d0) go to 100
   80 x = sign(xa,x)
   90 ix = 0
  100 return
      end
