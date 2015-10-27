*deck besi1
      function besi1 (x)
c***begin prologue  besi1
c***purpose  compute the modified (hyperbolic) bessel function of the
c            first kind of order one.
c***library   slatec (fnlib)
c***category  c10b1
c***type      single precision (besi1-s, dbesi1-d)
c***keywords  first kind, fnlib, hyperbolic bessel function,
c             modified bessel function, order one, special functions
c***author  fullerton, w., (lanl)
c***description
c
c besi1(x) calculates the modified (hyperbolic) bessel function
c of the first kind of order one for real argument x.
c
c series for bi1        on the interval  0.          to  9.00000d+00
c                                        with weighted error   2.40e-17
c                                         log weighted error  16.62
c                               significant figures required  16.23
c                                    decimal places required  17.14
c
c***references  (none)
c***routines called  besi1e, csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  besi1
      dimension bi1cs(11)
      logical first
      save bi1cs, nti1, xmin, xsml, xmax, first
      data bi1cs( 1) /   -.0019717132 61099859e0 /
      data bi1cs( 2) /    .4073488766 7546481e0 /
      data bi1cs( 3) /    .0348389942 99959456e0 /
      data bi1cs( 4) /    .0015453945 56300123e0 /
      data bi1cs( 5) /    .0000418885 21098377e0 /
      data bi1cs( 6) /    .0000007649 02676483e0 /
      data bi1cs( 7) /    .0000000100 42493924e0 /
      data bi1cs( 8) /    .0000000000 99322077e0 /
      data bi1cs( 9) /    .0000000000 00766380e0 /
      data bi1cs(10) /    .0000000000 00004741e0 /
      data bi1cs(11) /    .0000000000 00000024e0 /
      data first /.true./
c***first executable statement  besi1
      if (first) then
         nti1 = inits (bi1cs, 11, 0.1*r1mach(3))
         xmin = 2.0*r1mach(1)
         xsml = sqrt (4.5*r1mach(3))
         xmax = log (r1mach(2))
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.3.0) go to 20
c
      besi1 = 0.0
      if (y.eq.0.0)  return
c
      if (y .le. xmin) call xermsg ('slatec', 'besi1',
     +   'abs(x) so small i1 underflows', 1, 1)
      if (y.gt.xmin) besi1 = 0.5*x
      if (y.gt.xsml) besi1 = x * (.875 + csevl(y*y/4.5-1., bi1cs, nti1))
      return
c
 20   if (y .gt. xmax) call xermsg ('slatec', 'besi1',
     +   'abs(x) so big i1 overflows', 2, 2)
c
      besi1 = exp(y) * besi1e(x)
c
      return
      end
