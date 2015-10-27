*deck dxadd
      subroutine dxadd (x, ix, y, iy, z, iz, ierror)
c***begin prologue  dxadd
c***purpose  to provide double-precision floating-point arithmetic
c            with an extended exponent range.
c***library   slatec
c***category  a3d
c***type      double precision (xadd-s, dxadd-d)
c***keywords  extended-range double-precision arithmetic
c***author  lozier, daniel w., (national bureau of standards)
c           smith, john m., (nbs and george mason university)
c***description
c     double precision x, y, z
c     integer ix, iy, iz
c
c                  forms the extended-range sum  (z,iz) =
c                  (x,ix) + (y,iy).  (z,iz) is adjusted
c                  before returning. the input operands
c                  need not be in adjusted form, but their
c                  principal parts must satisfy
c                  radix**(-2l).le.abs(x).le.radix**(2l),
c                  radix**(-2l).le.abs(y).le.radix**(2l).
c
c***see also  dxset
c***references  (none)
c***routines called  dxadj
c***common blocks    dxblk2
c***revision history  (yymmdd)
c   820712  date written
c   881020  revised to meet slatec cml recommendations.  (dwl and jms)
c   901019  revisions to prologue.  (dwl and wrb)
c   901106  changed all specific intrinsics to generic.  (wrb)
c           corrected order of sections in prologue and added type
c           section.  (wrb)
c   920127  revised purpose section of prologue.  (dwl)
c***end prologue  dxadd
      double precision x, y, z
      integer ix, iy, iz
      double precision radix, radixl, rad2l, dlg10r
      integer l, l2, kmax
      common /dxblk2/ radix, radixl, rad2l, dlg10r, l, l2, kmax
      save /dxblk2/
      double precision s, t
c
c   the conditions imposed on l and kmax by this subroutine
c are
c     (1) 1 .lt. l .le. 0.5d0*logr(0.5d0*dzero)
c
c     (2) nradpl .lt. l .le. kmax/6
c
c     (3) kmax .le. (2**nbits - 4*l - 1)/2
c
c these conditions must be met by appropriate coding
c in subroutine dxset.
c
c***first executable statement  dxadd
      ierror=0
      if (x.ne.0.0d0) go to 10
      z = y
      iz = iy
      go to 220
   10 if (y.ne.0.0d0) go to 20
      z = x
      iz = ix
      go to 220
   20 continue
      if (ix.ge.0 .and. iy.ge.0) go to 40
      if (ix.lt.0 .and. iy.lt.0) go to 40
      if (abs(ix).le.6*l .and. abs(iy).le.6*l) go to 40
      if (ix.ge.0) go to 30
      z = y
      iz = iy
      go to 220
   30 continue
      z = x
      iz = ix
      go to 220
   40 i = ix - iy
      if (i) 80, 50, 90
   50 if (abs(x).gt.1.0d0 .and. abs(y).gt.1.0d0) go to 60
      if (abs(x).lt.1.0d0 .and. abs(y).lt.1.0d0) go to 70
      z = x + y
      iz = ix
      go to 220
   60 s = x/radixl
      t = y/radixl
      z = s + t
      iz = ix + l
      go to 220
   70 s = x*radixl
      t = y*radixl
      z = s + t
      iz = ix - l
      go to 220
   80 s = y
      is = iy
      t = x
      go to 100
   90 s = x
      is = ix
      t = y
  100 continue
c
c  at this point, the one of (x,ix) or (y,iy) that has the
c larger auxiliary index is stored in (s,is). the principal
c part of the other input is stored in t.
c
      i1 = abs(i)/l
      i2 = mod(abs(i),l)
      if (abs(t).ge.radixl) go to 130
      if (abs(t).ge.1.0d0) go to 120
      if (radixl*abs(t).ge.1.0d0) go to 110
      j = i1 + 1
      t = t*radix**(l-i2)
      go to 140
  110 j = i1
      t = t*radix**(-i2)
      go to 140
  120 j = i1 - 1
      if (j.lt.0) go to 110
      t = t*radix**(-i2)/radixl
      go to 140
  130 j = i1 - 2
      if (j.lt.0) go to 120
      t = t*radix**(-i2)/rad2l
  140 continue
c
c  at this point, some or all of the difference in the
c auxiliary indices has been used to effect a left shift
c of t.  the shifted value of t satisfies
c
c       radix**(-2*l) .le. abs(t) .le. 1.0d0
c
c and, if j=0, no further shifting remains to be done.
c
      if (j.eq.0) go to 190
      if (abs(s).ge.radixl .or. j.gt.3) go to 150
      if (abs(s).ge.1.0d0) go to (180, 150, 150), j
      if (radixl*abs(s).ge.1.0d0) go to (180, 170, 150), j
      go to (180, 170, 160), j
  150 z = s
      iz = is
      go to 220
  160 s = s*radixl
  170 s = s*radixl
  180 s = s*radixl
  190 continue
c
c   at this point, the remaining difference in the
c auxiliary indices has been used to effect a right shift
c of s.  if the shifted value of s would have exceeded
c radix**l, then (s,is) is returned as the value of the
c sum.
c
      if (abs(s).gt.1.0d0 .and. abs(t).gt.1.0d0) go to 200
      if (abs(s).lt.1.0d0 .and. abs(t).lt.1.0d0) go to 210
      z = s + t
      iz = is - j*l
      go to 220
  200 s = s/radixl
      t = t/radixl
      z = s + t
      iz = is - j*l + l
      go to 220
  210 s = s*radixl
      t = t*radixl
      z = s + t
      iz = is - j*l - l
  220 call dxadj(z, iz,ierror)
      return
      end
