*deck xadj
      subroutine xadj (x, ix, ierror)
c***begin prologue  xadj
c***purpose  to provide single-precision floating-point arithmetic
c            with an extended exponent range.
c***library   slatec
c***category  a3d
c***type      single precision (xadj-s, dxadj-d)
c***keywords  extended-range single-precision arithmetic
c***author  lozier, daniel w., (national bureau of standards)
c           smith, john m., (nbs and george mason university)
c***description
c     real x
c     integer ix
c
c                  transforms (x,ix) so that
c                  radix**(-l) .le. abs(x) .lt. radix**l.
c                  on most computers this transformation does
c                  not change the mantissa of x provided radix is
c                  the number base of single-precision arithmetic.
c
c***see also  xset
c***references  (none)
c***routines called  xermsg
c***common blocks    xblk2
c***revision history  (yymmdd)
c   820712  date written
c   881020  revised to meet slatec cml recommendations.  (dwl and jms)
c   901019  revisions to prologue.  (dwl and wrb)
c   901106  changed all specific intrinsics to generic.  (wrb)
c           corrected order of sections in prologue and added type
c           section.  (wrb)
c           calls to xerror changed to calls to xermsg.  (wrb)
c   920127  revised purpose section of prologue.  (dwl)
c***end prologue  xadj
      real x
      integer ix
      real radix, radixl, rad2l, dlg10r
      integer l, l2, kmax
      common /xblk2/ radix, radixl, rad2l, dlg10r, l, l2, kmax
      save /xblk2/
c
c   the condition imposed on l and kmax by this subroutine
c is
c     2*l .le. kmax
c
c this condition must be met by appropriate coding
c in subroutine xset.
c
c***first executable statement  xadj
      ierror=0
      if (x.eq.0.0) go to 50
      if (abs(x).ge.1.0) go to 20
      if (radixl*abs(x).ge.1.0) go to 60
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
   40 call xermsg ('slatec', 'xadj', 'overflow in auxiliary index', 107,
     +             1)
      ierror=107
      return
   50 ix = 0
   60 if (abs(ix).gt.kmax) go to 40
   70 return
      end
