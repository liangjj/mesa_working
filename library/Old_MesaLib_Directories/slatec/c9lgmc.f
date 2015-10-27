*deck c9lgmc
      complex function c9lgmc (zin)
c***begin prologue  c9lgmc
c***subsidiary
c***purpose  compute the log gamma correction factor so that
c            log(cgamma(z)) = 0.5*log(2.*pi) + (z-0.5)*log(z) - z
c            + c9lgmc(z).
c***library   slatec (fnlib)
c***category  c7a
c***type      complex (r9lgmc-s, d9lgmc-d, c9lgmc-c)
c***keywords  complete gamma function, correction term, fnlib,
c             log gamma, logarithm, special functions
c***author  fullerton, w., (lanl)
c***description
c
c compute the log gamma correction term for large abs(z) when real(z)
c .ge. 0.0 and for large abs(aimag(y)) when real(z) .lt. 0.0.  we find
c c9lgmc so that
c   log(z) = 0.5*log(2.*pi) + (z-0.5)*log(z) - z + c9lgmc(z)
c
c***references  (none)
c***routines called  r1mach, xermsg
c***revision history  (yymmdd)
c   780401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   900720  routine changed from user-callable to subsidiary.  (wrb)
c***end prologue  c9lgmc
      complex zin, z, z2inv
      dimension bern(11)
      logical first
      save bern, nterm, bound, xbig, xmax, first
      data bern( 1) /    .08333333333 3333333e0   /
      data bern( 2) /   -.002777777777 7777778e0  /
      data bern( 3) /    .0007936507936 5079365e0 /
      data bern( 4) /   -.0005952380952 3809524e0 /
      data bern( 5) /    .0008417508417 5084175e0 /
      data bern( 6) /   -.001917526917 5269175e0  /
      data bern( 7) /    .006410256410 2564103e0  /
      data bern( 8) /   -.02955065359 4771242e0   /
      data bern( 9) /    .1796443723 6883057e0    /
      data bern(10) /  -1.392432216 9059011e0     /
      data bern(11) /  13.40286404 4168392e0      /
      data first /.true./
c***first executable statement  c9lgmc
      if (first) then
         nterm = -0.30*log(r1mach(3))
         bound = 0.1170*nterm*(0.1*r1mach(3))**(-1./(2*nterm-1))
         xbig = 1.0/sqrt(r1mach(3))
         xmax = exp (min(log(r1mach(2)/12.0), -log(12.*r1mach(1))) )
      endif
      first = .false.
c
      z = zin
      x = real (z)
      y = aimag(z)
      cabsz = abs(z)
c
      if (x .lt. 0.0 .and. abs(y) .lt. bound) call xermsg ('slatec',
     +   'c9lgmc', 'not valid for negative real(z) and small ' //
     +   'abs(aimag(z))', 2, 2)
      if (cabsz .lt. bound) call xermsg ('slatec', 'c9lgmc',
     +   'not valid for small abs(z)', 3, 2)
c
      if (cabsz.ge.xmax) go to 50
c
      if (cabsz.ge.xbig) c9lgmc = 1.0/(12.0*z)
      if (cabsz.ge.xbig) return
c
      z2inv = 1.0/z**2
      c9lgmc = (0.0, 0.0)
      do 40 i=1,nterm
        ndx = nterm + 1 - i
        c9lgmc = bern(ndx) + c9lgmc*z2inv
 40   continue
c
      c9lgmc = c9lgmc/z
      return
c
 50   c9lgmc = (0.0, 0.0)
      call xermsg ('slatec', 'c9lgmc', 'z so big c9lgmc underflows', 1,
     +   1)
      return
c
      end
