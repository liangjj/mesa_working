      function clnrel(z)
c***begin prologue  clnrel
c***date written   770401   (yymmdd)
c***revision date  820801   (yymmdd)
c***category no.  c4b
c***keywords  complex,elementary function,logarithm,relative error
c***author  fullerton, w., (lanl)
c***purpose  computes the principal value of the complex natural
c            logarithm of 1+z with relative error accuracy for small
c            cabs(z).
c***description
c
c clnrel(z) = clog(1+z) with relative error accuracy near z = 0.
c let   rho = cabs(z)  and
c       r**2 = cabs(1+z)**2 = (1+x)**2 + y**2 = 1 + 2*x + rho**2 .
c now if rho is small we may evaluate clnrel(z) accurately by
c       clog(1+z) = cmplx  (alog(r), carg(1+z))
c                 = cmplx  (0.5*alog(r**2), carg(1+z))
c                 = cmplx  (0.5*alnrel(2*x+rho**2), carg(1+z))
c***references  (none)
c***routines called  alnrel,carg,r1mach,xerror
c***end prologue  clnrel
      implicit real*8(a-h,o-z)
      complex*16 z
      data sqeps /0.0d0/
c***first executable statement  clnrel
      if (sqeps.eq.0.d0) sqeps = sqrt (r1mach(4))
c
      if (abs(1.d0+z).lt.sqeps) call lnkerr ('clnrel  answer lt half pre
     1cision because z too near -1')
c
      rho = abs(z)
      if (rho.gt.0.375d0) clnrel = log (1.0d0+z)
      if (rho.gt.0.375d0) return
c
      x = real(z)
      clnrel = dcmplx (0.5d0*alnrel(2.d0*x+rho**2), carg(1.0d0+z))
c
      return
      end
