*deck ai
      function ai (x)
c***begin prologue  ai
c***purpose  evaluate the airy function.
c***library   slatec (fnlib)
c***category  c10d
c***type      single precision (ai-s, dai-d)
c***keywords  airy function, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c ai(x) computes the airy function ai(x)
c series for aif        on the interval -1.00000d+00 to  1.00000d+00
c                                        with weighted error   1.09e-19
c                                         log weighted error  18.96
c                               significant figures required  17.76
c                                    decimal places required  19.44
c
c series for aig        on the interval -1.00000d+00 to  1.00000d+00
c                                        with weighted error   1.51e-17
c                                         log weighted error  16.82
c                               significant figures required  15.19
c                                    decimal places required  17.27
c
c***references  (none)
c***routines called  aie, csevl, inits, r1mach, r9aimp, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920618  removed space from variable names.  (rwc, wrb)
c***end prologue  ai
      dimension aifcs(9), aigcs(8)
      logical first
      save aifcs, aigcs, naif, naig, x3sml, xmax, first
      data aifcs( 1) /   -.0379713584 9666999750e0 /
      data aifcs( 2) /    .0591918885 3726363857e0 /
      data aifcs( 3) /    .0009862928 0577279975e0 /
      data aifcs( 4) /    .0000068488 4381907656e0 /
      data aifcs( 5) /    .0000000259 4202596219e0 /
      data aifcs( 6) /    .0000000000 6176612774e0 /
      data aifcs( 7) /    .0000000000 0010092454e0 /
      data aifcs( 8) /    .0000000000 0000012014e0 /
      data aifcs( 9) /    .0000000000 0000000010e0 /
      data aigcs( 1) /    .0181523655 8116127e0 /
      data aigcs( 2) /    .0215725631 6601076e0 /
      data aigcs( 3) /    .0002567835 6987483e0 /
      data aigcs( 4) /    .0000014265 2141197e0 /
      data aigcs( 5) /    .0000000045 7211492e0 /
      data aigcs( 6) /    .0000000000 0952517e0 /
      data aigcs( 7) /    .0000000000 0001392e0 /
      data aigcs( 8) /    .0000000000 0000001e0 /
      data first /.true./
c***first executable statement  ai
      if (first) then
         naif = inits (aifcs, 9, 0.1*r1mach(3))
         naig = inits (aigcs, 8, 0.1*r1mach(3))
c
         x3sml = r1mach(3)**0.3334
         xmaxt = (-1.5*log(r1mach(1)))**0.6667
         xmax = xmaxt - xmaxt*log(xmaxt)/
     *                   (4.0*sqrt(xmaxt)+1.0) - 0.01
      endif
      first = .false.
c
      if (x.ge.(-1.0)) go to 20
      call r9aimp (x, xm, theta)
      ai = xm * cos(theta)
      return
c
 20   if (x.gt.1.0) go to 30
      z = 0.0
      if (abs(x).gt.x3sml) z = x**3
      ai = 0.375 + (csevl (z, aifcs, naif) - x*(0.25 +
     1  csevl (z, aigcs, naig)) )
      return
c
 30   if (x.gt.xmax) go to 40
      ai = aie(x) * exp(-2.0*x*sqrt(x)/3.0)
      return
c
 40   ai = 0.0
      call xermsg ('slatec', 'ai', 'x so big ai underflows', 1, 1)
      return
c
      end
