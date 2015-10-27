*deck atanh
      function atanh (x)
c***begin prologue  atanh
c***purpose  compute the arc hyperbolic tangent.
c***library   slatec (fnlib)
c***category  c4c
c***type      single precision (atanh-s, datanh-d, catanh-c)
c***keywords  arc hyperbolic tangent, atanh, elementary functions,
c             fnlib, inverse hyperbolic tangent
c***author  fullerton, w., (lanl)
c***description
c
c atanh(x) computes the arc hyperbolic tangent of x.
c
c series for atnh       on the interval  0.          to  2.50000d-01
c                                        with weighted error   6.70e-18
c                                         log weighted error  17.17
c                               significant figures required  16.01
c                                    decimal places required  17.76
c
c***references  (none)
c***routines called  csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  atanh
      dimension atnhcs(15)
      logical first
      save atnhcs, nterms, dxrel, sqeps, first
      data atnhcs( 1) /    .0943951023 93195492e0 /
      data atnhcs( 2) /    .0491984370 55786159e0 /
      data atnhcs( 3) /    .0021025935 22455432e0 /
      data atnhcs( 4) /    .0001073554 44977611e0 /
      data atnhcs( 5) /    .0000059782 67249293e0 /
      data atnhcs( 6) /    .0000003505 06203088e0 /
      data atnhcs( 7) /    .0000000212 63743437e0 /
      data atnhcs( 8) /    .0000000013 21694535e0 /
      data atnhcs( 9) /    .0000000000 83658755e0 /
      data atnhcs(10) /    .0000000000 05370503e0 /
      data atnhcs(11) /    .0000000000 00348665e0 /
      data atnhcs(12) /    .0000000000 00022845e0 /
      data atnhcs(13) /    .0000000000 00001508e0 /
      data atnhcs(14) /    .0000000000 00000100e0 /
      data atnhcs(15) /    .0000000000 00000006e0 /
      data first /.true./
c***first executable statement  atanh
      if (first) then
         nterms = inits (atnhcs, 15, 0.1*r1mach(3))
         dxrel = sqrt (r1mach(4))
         sqeps = sqrt (3.0*r1mach(3))
      endif
      first = .false.
c
      y = abs(x)
      if (y .ge. 1.0) call xermsg ('slatec', 'atanh', 'abs(x) ge 1', 2,
     +   2)
c
      if (1.0-y .lt. dxrel) call xermsg ('slatec', 'atanh',
     +   'answer lt half precision because abs(x) too near 1', 1, 1)
c
      atanh = x
      if (y.gt.sqeps .and. y.le.0.5) atanh = x*(1.0 + csevl (8.*x*x-1.,
     1  atnhcs, nterms))
      if (y.gt.0.5) atanh = 0.5*log((1.0+x)/(1.0-x))
c
      return
      end
