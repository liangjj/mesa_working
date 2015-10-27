*deck clnrel
      complex function clnrel (z)
c***begin prologue  clnrel
c***purpose  evaluate ln(1+x) accurate in the sense of relative error.
c***library   slatec (fnlib)
c***category  c4b
c***type      complex (alnrel-s, dlnrel-d, clnrel-c)
c***keywords  elementary functions, fnlib, logarithm
c***author  fullerton, w., (lanl)
c***description
c
c clnrel(z) = log(1+z) with relative error accuracy near z = 0.
c let   rho = abs(z)  and
c       r**2 = abs(1+z)**2 = (1+x)**2 + y**2 = 1 + 2*x + rho**2 .
c now if rho is small we may evaluate clnrel(z) accurately by
c       log(1+z) = cmplx  (log(r), carg(1+z))
c                 = cmplx  (0.5*log(r**2), carg(1+z))
c                 = cmplx  (0.5*alnrel(2*x+rho**2), carg(1+z))
c
c***references  (none)
c***routines called  alnrel, carg, r1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  clnrel
      complex z
      save sqeps
      data sqeps /0.0/
c***first executable statement  clnrel
      if (sqeps.eq.0.) sqeps = sqrt (r1mach(4))
c
      if (abs(1.+z) .lt. sqeps) call xermsg ('slatec', 'clnrel',
     +   'answer lt half precision because z too near -1', 1, 1)
c
      rho = abs(z)
      if (rho.gt.0.375) clnrel = log (1.0+z)
      if (rho.gt.0.375) return
c
      x = real(z)
      clnrel = cmplx (0.5*alnrel(2.*x+rho**2), carg(1.0+z))
c
      return
      end
