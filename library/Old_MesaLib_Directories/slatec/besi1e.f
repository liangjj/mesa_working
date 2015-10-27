*deck besi1e
      function besi1e (x)
c***begin prologue  besi1e
c***purpose  compute the exponentially scaled modified (hyperbolic)
c            bessel function of the first kind of order one.
c***library   slatec (fnlib)
c***category  c10b1
c***type      single precision (besi1e-s, dbsi1e-d)
c***keywords  exponentially scaled, first kind, fnlib,
c             hyperbolic bessel function, modified bessel function,
c             order one, special functions
c***author  fullerton, w., (lanl)
c***description
c
c besi1e(x) calculates the exponentially scaled modified (hyperbolic)
c bessel function of the first kind of order one for real argument x;
c i.e., exp(-abs(x))*i1(x).
c
c series for bi1        on the interval  0.          to  9.00000d+00
c                                        with weighted error   2.40e-17
c                                         log weighted error  16.62
c                               significant figures required  16.23
c                                    decimal places required  17.14
c
c series for ai1        on the interval  1.25000d-01 to  3.33333d-01
c                                        with weighted error   6.98e-17
c                                         log weighted error  16.16
c                               significant figures required  14.53
c                                    decimal places required  16.82
c
c series for ai12       on the interval  0.          to  1.25000d-01
c                                        with weighted error   3.55e-17
c                                         log weighted error  16.45
c                               significant figures required  14.69
c                                    decimal places required  17.12
c
c***references  (none)
c***routines called  csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890210  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920618  removed space from variable names.  (rwc, wrb)
c***end prologue  besi1e
      dimension bi1cs(11), ai1cs(21), ai12cs(22)
      logical first
      save bi1cs, ai1cs, ai12cs, nti1, ntai1, ntai12, xmin, xsml, first
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
      data ai1cs( 1) /   -.0284674418 1881479e0 /
      data ai1cs( 2) /   -.0192295323 1443221e0 /
      data ai1cs( 3) /   -.0006115185 8579437e0 /
      data ai1cs( 4) /   -.0000206997 1253350e0 /
      data ai1cs( 5) /    .0000085856 1914581e0 /
      data ai1cs( 6) /    .0000010494 9824671e0 /
      data ai1cs( 7) /   -.0000002918 3389184e0 /
      data ai1cs( 8) /   -.0000000155 9378146e0 /
      data ai1cs( 9) /    .0000000131 8012367e0 /
      data ai1cs(10) /   -.0000000014 4842341e0 /
      data ai1cs(11) /   -.0000000002 9085122e0 /
      data ai1cs(12) /    .0000000001 2663889e0 /
      data ai1cs(13) /   -.0000000000 1664947e0 /
      data ai1cs(14) /   -.0000000000 0166665e0 /
      data ai1cs(15) /    .0000000000 0124260e0 /
      data ai1cs(16) /   -.0000000000 0027315e0 /
      data ai1cs(17) /    .0000000000 0002023e0 /
      data ai1cs(18) /    .0000000000 0000730e0 /
      data ai1cs(19) /   -.0000000000 0000333e0 /
      data ai1cs(20) /    .0000000000 0000071e0 /
      data ai1cs(21) /   -.0000000000 0000006e0 /
      data ai12cs( 1) /    .0285762350 1828014e0 /
      data ai12cs( 2) /   -.0097610974 9136147e0 /
      data ai12cs( 3) /   -.0001105889 3876263e0 /
      data ai12cs( 4) /   -.0000038825 6480887e0 /
      data ai12cs( 5) /   -.0000002512 2362377e0 /
      data ai12cs( 6) /   -.0000000263 1468847e0 /
      data ai12cs( 7) /   -.0000000038 3538039e0 /
      data ai12cs( 8) /   -.0000000005 5897433e0 /
      data ai12cs( 9) /   -.0000000000 1897495e0 /
      data ai12cs(10) /    .0000000000 3252602e0 /
      data ai12cs(11) /    .0000000000 1412580e0 /
      data ai12cs(12) /    .0000000000 0203564e0 /
      data ai12cs(13) /   -.0000000000 0071985e0 /
      data ai12cs(14) /   -.0000000000 0040836e0 /
      data ai12cs(15) /   -.0000000000 0002101e0 /
      data ai12cs(16) /    .0000000000 0004273e0 /
      data ai12cs(17) /    .0000000000 0001041e0 /
      data ai12cs(18) /   -.0000000000 0000382e0 /
      data ai12cs(19) /   -.0000000000 0000186e0 /
      data ai12cs(20) /    .0000000000 0000033e0 /
      data ai12cs(21) /    .0000000000 0000028e0 /
      data ai12cs(22) /   -.0000000000 0000003e0 /
      data first /.true./
c***first executable statement  besi1e
      if (first) then
         nti1 = inits (bi1cs, 11, 0.1*r1mach(3))
         ntai1 = inits (ai1cs, 21, 0.1*r1mach(3))
         ntai12 = inits (ai12cs, 22, 0.1*r1mach(3))
c
         xmin = 2.0*r1mach(1)
         xsml = sqrt (4.5*r1mach(3))
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.3.0) go to 20
c
      besi1e = 0.0
      if (y.eq.0.0)  return
c
      if (y .le. xmin) call xermsg ('slatec', 'besi1e',
     +   'abs(x) so small i1 underflows', 1, 1)
      if (y.gt.xmin) besi1e = 0.5*x
      if (y.gt.xsml) besi1e = x * (.875 + csevl(y*y/4.5-1., bi1cs,nti1))
      besi1e = exp(-y) * besi1e
      return
c
 20   if (y.le.8.) besi1e = (.375 + csevl ((48./y-11.)/5., ai1cs, ntai1)
     1  ) / sqrt(y)
      if (y.gt.8.) besi1e = (.375 + csevl (16./y-1.0, ai12cs, ntai12))
     1  / sqrt(y)
      besi1e = sign (besi1e, x)
c
      return
      end
