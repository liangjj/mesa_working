*deck dxc210
      subroutine dxc210 (k, z, j, ierror)
c***begin prologue  dxc210
c***purpose  to provide double-precision floating-point arithmetic
c            with an extended exponent range.
c***library   slatec
c***category  a3d
c***type      double precision (xc210-s, dxc210-d)
c***keywords  extended-range double-precision arithmetic
c***author  lozier, daniel w., (national bureau of standards)
c           smith, john m., (nbs and george mason university)
c***description
c     integer k, j
c     double precision z
c
c                  given k this subroutine computes j and z
c                  such that  radix**k = z*10**j, where z is in
c                  the range 1/10 .le. z .lt. 1.
c                  the value of z will be accurate to full
c                  double-precision provided the number
c                  of decimal places in the largest
c                  integer plus the number of decimal
c                  places carried in double-precision does not
c                  exceed 60. dxc210 is called by subroutine
c                  dxcon when necessary. the user should
c                  never need to call dxc210 directly.
c
c***see also  dxset
c***references  (none)
c***routines called  xermsg
c***common blocks    dxblk3
c***revision history  (yymmdd)
c   820712  date written
c   890126  revised to meet slatec cml recommendations.  (dwl and jms)
c   901019  revisions to prologue.  (dwl and wrb)
c   901106  changed all specific intrinsics to generic.  (wrb)
c           corrected order of sections in prologue and added type
c           section.  (wrb)
c           calls to xerror changed to calls to xermsg.  (wrb)
c   920127  revised purpose section of prologue.  (dwl)
c***end prologue  dxc210
      double precision z
      integer k, j
      integer nlg102, mlg102, lg102
      common /dxblk3/ nlg102, mlg102, lg102(21)
      save /dxblk3/
c
c   the conditions imposed on nlg102, mlg102, and lg102 by
c this subroutine are
c
c     (1) nlg102 .ge. 2
c
c     (2) mlg102 .ge. 1
c
c     (3) 2*mlg102*(mlg102 - 1) .le. 2**nbits - 1
c
c these conditions must be met by appropriate coding
c in subroutine dxset.
c
c***first executable statement  dxc210
      ierror=0
      if (k.eq.0) go to 70
      m = mlg102
      ka = abs(k)
      ka1 = ka/m
      ka2 = mod(ka,m)
      if (ka1.ge.m) go to 60
      nm1 = nlg102 - 1
      np1 = nlg102 + 1
      it = ka2*lg102(np1)
      ic = it/m
      id = mod(it,m)
      z = id
      if (ka1.gt.0) go to 20
      do 10 ii=1,nm1
        i = np1 - ii
        it = ka2*lg102(i) + ic
        ic = it/m
        id = mod(it,m)
        z = z/m + id
   10 continue
      ja = ka*lg102(1) + ic
      go to 40
   20 continue
      do 30 ii=1,nm1
        i = np1 - ii
        it = ka2*lg102(i) + ka1*lg102(i+1) + ic
        ic = it/m
        id = mod(it,m)
        z = z/m + id
   30 continue
      ja = ka*lg102(1) + ka1*lg102(2) + ic
   40 continue
      z = z/m
      if (k.gt.0) go to 50
      j = -ja
      z = 10.0d0**(-z)
      go to 80
   50 continue
      j = ja + 1
      z = 10.0d0**(z-1.0d0)
      go to 80
   60 continue
c   this error occurs if k exceeds  mlg102**2 - 1  in magnitude.
c
      call xermsg ('slatec', 'dxc210', 'k too large', 208, 1)
      ierror=208
      return
   70 continue
      j = 0
      z = 1.0d0
   80 return
      end
