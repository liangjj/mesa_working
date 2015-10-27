*deck cpsi
      complex function cpsi (zin)
c***begin prologue  cpsi
c***purpose  compute the psi (or digamma) function.
c***library   slatec (fnlib)
c***category  c7c
c***type      complex (psi-s, dpsi-d, cpsi-c)
c***keywords  digamma function, fnlib, psi function, special functions
c***author  fullerton, w., (lanl)
c***description
c
c psi(x) calculates the psi (or digamma) function of x.  psi(x)
c is the logarithmic derivative of the gamma function of x.
c
c***references  (none)
c***routines called  ccot, r1mach, xermsg
c***revision history  (yymmdd)
c   780501  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900727  added external statement.  (wrb)
c***end prologue  cpsi
      complex zin, z, z2inv, corr, ccot
      dimension bern(13)
      logical first
      external ccot
      save bern, pi, nterm, bound, dxrel, rmin, rbig, first
      data bern( 1) /   .8333333333 3333333 e-1 /
      data bern( 2) /  -.8333333333 3333333 e-2 /
      data bern( 3) /   .3968253968 2539683 e-2 /
      data bern( 4) /  -.4166666666 6666667 e-2 /
      data bern( 5) /   .7575757575 7575758 e-2 /
      data bern( 6) /  -.2109279609 2796093 e-1 /
      data bern( 7) /   .8333333333 3333333 e-1 /
      data bern( 8) /  -.4432598039 2156863 e0 /
      data bern( 9) /   .3053954330 2701197 e1 /
      data bern(10) /  -.2645621212 1212121 e2 /
      data bern(11) /   .2814601449 2753623 e3 /
      data bern(12) /  -.3454885393 7728938 e4 /
      data bern(13) /   .5482758333 3333333 e5 /
      data pi / 3.141592653 589793 e0 /
      data first /.true./
c***first executable statement  cpsi
      if (first) then
         nterm = -0.30*log(r1mach(3))
c maybe bound = n*(0.1*eps)**(-1/(2*n-1)) / (pi*exp(1))
         bound = 0.1171*nterm*(0.1*r1mach(3))**(-1.0/(2*nterm-1))
         dxrel = sqrt(r1mach(4))
         rmin = exp (max (log(r1mach(1)), -log(r1mach(2))) + 0.011 )
         rbig = 1.0/r1mach(3)
      endif
      first = .false.
c
      z = zin
      x = real(z)
      y = aimag(z)
      if (y.lt.0.0) z = conjg(z)
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
      corr = -pi*ccot(pi*z)
      z = 1.0 - z
      go to 50
c
c use the recursion relation for abs(z) small.
c
 20   if (cabsz .lt. rmin) call xermsg ('slatec', 'cpsi',
     +   'cpsi called with z so near 0 that cpsi overflows', 2, 2)
c
      if (x.ge.(-0.5) .or. abs(y).gt.dxrel) go to 30
      if (abs((z-aint(x-0.5))/x) .lt. dxrel) call xermsg ('slatec',
     +   'cpsi',
     +   'answer lt half precision because z too near negative integer',
     +   1, 1)
      if (y .eq. 0.0 .and. x .eq. aint(x)) call xermsg ('slatec',
     +   'cpsi', 'z is a negative integer', 3, 2)
c
 30   n = sqrt(bound**2-y**2) - x + 1.0
      do 40 i=1,n
        corr = corr - 1.0/z
        z = z + 1.0
 40   continue
c
c now evaluate the asymptotic series for suitably large z.
c
 50   if (cabsz.gt.rbig) cpsi = log(z) + corr
      if (cabsz.gt.rbig) go to 70
c
      cpsi = (0.0, 0.0)
      z2inv = 1.0/z**2
      do 60 i=1,nterm
        ndx = nterm + 1 - i
        cpsi = bern(ndx) + z2inv*cpsi
 60   continue
      cpsi = log(z) - 0.5/z - cpsi*z2inv + corr
c
 70   if (y.lt.0.0) cpsi = conjg(cpsi)
c
      return
      end
