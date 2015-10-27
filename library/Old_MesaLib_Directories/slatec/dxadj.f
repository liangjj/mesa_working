*deck dxadj
      subroutine dxadj (x, ix, ierror)
c***begin prologue  dxadj
c***purpose  to provide double-precision floating-point arithmetic
c            with an extended exponent range.
c***library   slatec
c***category  a3d
c***type      double precision (xadj-s, dxadj-d)
c***keywords  extended-range double-precision arithmetic
c***author  lozier, daniel w., (national bureau of standards)
c           smith, john m., (nbs and george mason university)
c***description
c     double precision x
c     integer ix
c
c                  transforms (x,ix) so that
c                  radix**(-l) .le. abs(x) .lt. radix**l.
c                  on most computers this transformation does
c                  not change the mantissa of x provided radix is
c                  the number base of double-precision arithmetic.
c
c***see also  dxset
c***references  (none)
c***routines called  xermsg
c***common blocks    dxblk2
c***revision history  (yymmdd)
c   820712  date written
c   881020  revised to meet slatec cml recommendations.  (dwl and jms)
c   901019  revisions to prologue.  (dwl and wrb)
c   901106  changed all specific intrinsics to generic.  (wrb)
c           corrected order of sections in prologue and added type
c           section.  (wrb)
c           calls to xerror changed to calls to xermsg.  (wrb)
c   920127  revised purpose section of prologue.  (dwl)
c***end prologue  dxadj
      double precision x
      integer ix
      double precision radix, radixl, rad2l, dlg10r
      integer l, l2, kmax
      common /dxblk2/ radix, radixl, rad2l, dlg10r, l, l2, kmax
      save /dxblk2/
c
c   the condition imposed on l and kmax by this subroutine
c is
c     2*l .le. kmax
c
c this condition must be met by appropriate coding
c in subroutine dxset.
c
c***first executable statement  dxadj
      ierror=0
      if (x.eq.0.0d0) go to 50
      if (abs(x).ge.1.0d0) go to 20
      if (radixl*abs(x).ge.1.0d0) go to 60
      x = x*rad2l
      if (ix.lt.0) go to 10
      ix = ix - l2
      go to 70
   10 if (ix.lt.-kmax+l2) go to 40
      ix = ix - l2
      go to 70
   20 if (abs(x).lt.radixl) go to 60
      x = x/rad2l
      if (ix.gt.0) go to 30
      ix = ix + l2
      go to 70
   30 if (ix.gt.kmax-l2) go to 40
      ix = ix + l2
      go to 70
   40 call xermsg ('slatec', 'dxadj', 'overflow in auxiliary index',
     +             207, 1)
      ierror=207
      return
   50 ix = 0
   60 if (abs(ix).gt.kmax) go to 40
   70 return
      end
