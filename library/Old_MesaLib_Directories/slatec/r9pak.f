*deck r9pak
      function r9pak (y, n)
c***begin prologue  r9pak
c***purpose  pack a base 2 exponent into a floating point number.
c***library   slatec (fnlib)
c***category  a6b
c***type      single precision (r9pak-s, d9pak-d)
c***keywords  fnlib, pack
c***author  fullerton, w., (lanl)
c***description
c
c pack a base 2 exponent into floating point number y.  this
c routine is almost the inverse of r9upak.  it is not exactly
c the inverse, because abs(x) need not be between 0.5 and
c 1.0.  if both r9pak and 2.0**n were known to be in range, we
c could compute
c       r9pak = y * 2.0**n .
c
c***references  (none)
c***routines called  i1mach, r1mach, r9upak, xermsg
c***revision history  (yymmdd)
c   790801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   901009  routine used i1mach(7) where it should use i1mach(10),
c           corrected (rwc)
c***end prologue  r9pak
      logical first
      save nmin, nmax, a1n210, first
      data a1n210 / 3.321928094 887362 e0/
      data first /.true./
c***first executable statement  r9pak
      if (first) then
         a1n2b = 1.0
         if (i1mach(10).ne.2) a1n2b = r1mach(5)*a1n210
         nmin = a1n2b*i1mach(12)
         nmax = a1n2b*i1mach(13)
      endif
      first = .false.
c
      call r9upak(y,r9pak,ny)
c
      nsum = n + ny
      if (nsum.lt.nmin) go to 40
      if (nsum .gt. nmax) call xermsg ('slatec', 'r9pak',
     +   'packed number overflows', 2, 2)
c
      if (nsum.eq.0) return
      if (nsum.gt.0) go to 30
c
 20   r9pak = 0.5*r9pak
      nsum = nsum + 1
      if(nsum.ne.0) go to 20
      return
c
30    r9pak = 2.0*r9pak
      nsum = nsum - 1
      if(nsum.ne.0) go to 30
      return
c
40    call xermsg ('slatec', 'r9pak', 'packed number underflows', 1, 1)
      r9pak = 0.0
      return
c
      end
