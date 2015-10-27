*deck clngam
      complex function clngam (zin)
c***begin prologue  clngam
c***purpose  compute the logarithm of the absolute value of the gamma
c            function.
c***library   slatec (fnlib)
c***category  c7a
c***type      complex (alngam-s, dlngam-d, clngam-c)
c***keywords  absolute value, complete gamma function, fnlib, logarithm,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c clngam computes the natural log of the complex valued gamma function
c at zin, where zin is a complex number.  this is a preliminary version,
c which is not accurate.
c
c***references  (none)
c***routines called  c9lgmc, carg, clnrel, r1mach, xermsg
c***revision history  (yymmdd)
c   780401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  clngam
      complex zin, z, corr, clnrel, c9lgmc
      logical first
      save pi, sq2pil, bound, dxrel, first
      data pi / 3.1415926535 8979324e0 /
      data sq2pil / 0.9189385332 0467274e0 /
      data first /.true./
c***first executable statement  clngam
      if (first) then
         n = -0.30*log(r1mach(3))
c bound = n*(0.1*eps)**(-1/(2*n-1))/(pi*exp(1))
         bound = 0.1171*n*(0.1*r1mach(3))**(-1./(2*n-1))
         dxrel = sqrt (r1mach(4))
      endif
      first = .false.
c
      z = zin
      x = real(zin)
      y = aimag(zin)
c
      corr = (0.0, 0.0)
      cabsz = abs(z)
      if (x.ge.0.0 .and. cabsz.gt.bound) go to 50
      if (x.lt.0.0 .and. abs(y).gt.bound) go to 50
c
      if (cabsz.lt.bound) go to 20
c
c use the reflection formula for real(z) negative, abs(z) large, and
c abs(aimag(y)) small.
c
      if (y.gt.0.0) z = conjg (z)
      corr = exp (-cmplx(0.0,2.0*pi)*z)
      if (real(corr) .eq. 1.0 .and. aimag(corr) .eq. 0.0) call xermsg
     +   ('slatec', 'clngam', 'z is a negative integer', 3, 2)
c
      clngam = sq2pil + 1.0 - cmplx(0.0,pi)*(z-0.5) - clnrel(-corr)
     1  + (z-0.5)*log(1.0-z) - z - c9lgmc(1.0-z)
      if (y.gt.0.0) clngam = conjg (clngam)
      return
c
c use the recursion relation for abs(z) small.
c
 20   if (x.ge.(-0.5) .or. abs(y).gt.dxrel) go to 30
      if (abs((z-aint(x-0.5))/x) .lt. dxrel) call xermsg ('slatec',
     +   'clngam',
     +   'answer lt half precision because z too near negative integer',
     +   1, 1)
c
 30   n = sqrt (bound**2 - y**2) - x + 1.0
      argsum = 0.0
      corr = (1.0, 0.0)
      do 40 i=1,n
        argsum = argsum + carg(z)
        corr = z*corr
        z = 1.0 + z
 40   continue
c
      if (real(corr) .eq. 0.0 .and. aimag(corr) .eq. 0.0) call xermsg
     +   ('slatec', 'clngam', 'z is a negative integer', 3, 2)
      corr = -cmplx (log(abs(corr)), argsum)
c
c use stirling-s approximation for large z.
c
 50   clngam = sq2pil + (z-0.5)*log(z) - z + corr + c9lgmc(z)
      return
c
      end
