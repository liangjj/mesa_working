*deck besy0
      function besy0 (x)
c***begin prologue  besy0
c***purpose  compute the bessel function of the second kind of order
c            zero.
c***library   slatec (fnlib)
c***category  c10a1
c***type      single precision (besy0-s, dbesy0-d)
c***keywords  bessel function, fnlib, order zero, second kind,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c besy0(x) calculates the bessel function of the second kind
c of order zero for real argument x.
c
c series for by0        on the interval  0.          to  1.60000d+01
c                                        with weighted error   1.20e-17
c                                         log weighted error  16.92
c                               significant figures required  16.15
c                                    decimal places required  17.48
c
c series for bm0        on the interval  0.          to  6.25000d-02
c                                        with weighted error   4.98e-17
c                                         log weighted error  16.30
c                               significant figures required  14.97
c                                    decimal places required  16.96
c
c series for bth0       on the interval  0.          to  6.25000d-02
c                                        with weighted error   3.67e-17
c                                         log weighted error  16.44
c                               significant figures required  15.53
c                                    decimal places required  17.13
c
c***references  (none)
c***routines called  besj0, csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c***end prologue  besy0
      dimension by0cs(13), bm0cs(21), bth0cs(24)
      logical first
      save by0cs, bm0cs, bth0cs, twodpi, pi4,
     1 nty0, ntm0, ntth0, xsml, xmax, first
      data by0cs( 1) /   -.0112778393 92865573e0 /
      data by0cs( 2) /   -.1283452375 6042035e0 /
      data by0cs( 3) /   -.1043788479 9794249e0 /
      data by0cs( 4) /    .0236627491 83969695e0 /
      data by0cs( 5) /   -.0020903916 47700486e0 /
      data by0cs( 6) /    .0001039754 53939057e0 /
      data by0cs( 7) /   -.0000033697 47162423e0 /
      data by0cs( 8) /    .0000000772 93842676e0 /
      data by0cs( 9) /   -.0000000013 24976772e0 /
      data by0cs(10) /    .0000000000 17648232e0 /
      data by0cs(11) /   -.0000000000 00188105e0 /
      data by0cs(12) /    .0000000000 00001641e0 /
      data by0cs(13) /   -.0000000000 00000011e0 /
      data bm0cs( 1) /    .0928496163 7381644e0 /
      data bm0cs( 2) /   -.0014298770 7403484e0 /
      data bm0cs( 3) /    .0000283057 9271257e0 /
      data bm0cs( 4) /   -.0000014330 0611424e0 /
      data bm0cs( 5) /    .0000001202 8628046e0 /
      data bm0cs( 6) /   -.0000000139 7113013e0 /
      data bm0cs( 7) /    .0000000020 4076188e0 /
      data bm0cs( 8) /   -.0000000003 5399669e0 /
      data bm0cs( 9) /    .0000000000 7024759e0 /
      data bm0cs(10) /   -.0000000000 1554107e0 /
      data bm0cs(11) /    .0000000000 0376226e0 /
      data bm0cs(12) /   -.0000000000 0098282e0 /
      data bm0cs(13) /    .0000000000 0027408e0 /
      data bm0cs(14) /   -.0000000000 0008091e0 /
      data bm0cs(15) /    .0000000000 0002511e0 /
      data bm0cs(16) /   -.0000000000 0000814e0 /
      data bm0cs(17) /    .0000000000 0000275e0 /
      data bm0cs(18) /   -.0000000000 0000096e0 /
      data bm0cs(19) /    .0000000000 0000034e0 /
      data bm0cs(20) /   -.0000000000 0000012e0 /
      data bm0cs(21) /    .0000000000 0000004e0 /
      data bth0cs( 1) /   -.2463916377 4300119e0 /
      data bth0cs( 2) /    .0017370983 07508963e0 /
      data bth0cs( 3) /   -.0000621836 33402968e0 /
      data bth0cs( 4) /    .0000043680 50165742e0 /
      data bth0cs( 5) /   -.0000004560 93019869e0 /
      data bth0cs( 6) /    .0000000621 97400101e0 /
      data bth0cs( 7) /   -.0000000103 00442889e0 /
      data bth0cs( 8) /    .0000000019 79526776e0 /
      data bth0cs( 9) /   -.0000000004 28198396e0 /
      data bth0cs(10) /    .0000000001 02035840e0 /
      data bth0cs(11) /   -.0000000000 26363898e0 /
      data bth0cs(12) /    .0000000000 07297935e0 /
      data bth0cs(13) /   -.0000000000 02144188e0 /
      data bth0cs(14) /    .0000000000 00663693e0 /
      data bth0cs(15) /   -.0000000000 00215126e0 /
      data bth0cs(16) /    .0000000000 00072659e0 /
      data bth0cs(17) /   -.0000000000 00025465e0 /
      data bth0cs(18) /    .0000000000 00009229e0 /
      data bth0cs(19) /   -.0000000000 00003448e0 /
      data bth0cs(20) /    .0000000000 00001325e0 /
      data bth0cs(21) /   -.0000000000 00000522e0 /
      data bth0cs(22) /    .0000000000 00000210e0 /
      data bth0cs(23) /   -.0000000000 00000087e0 /
      data bth0cs(24) /    .0000000000 00000036e0 /
      data twodpi / 0.6366197723 6758134e0 /
      data pi4 / 0.7853981633 9744831e0 /
      data first /.true./
c***first executable statement  besy0
      if (first) then
         nty0 = inits (by0cs, 13, 0.1*r1mach(3))
         ntm0 = inits (bm0cs, 21, 0.1*r1mach(3))
         ntth0 = inits (bth0cs, 24, 0.1*r1mach(3))
c
         xsml = sqrt (4.0*r1mach(3))
         xmax = 1.0/r1mach(4)
      endif
      first = .false.
c
      if (x .le. 0.) call xermsg ('slatec', 'besy0',
     +   'x is zero or negative', 1, 2)
      if (x.gt.4.0) go to 20
c
      y = 0.
      if (x.gt.xsml) y = x*x
      besy0 = twodpi*log(0.5*x)*besj0(x) + .375 + csevl (.125*y-1.,
     1  by0cs, nty0)
      return
c
 20   if (x .gt. xmax) call xermsg ('slatec', 'besy0',
     +   'no precision because x is big', 2, 2)
c
      z = 32.0/x**2 - 1.0
      ampl = (0.75 + csevl (z, bm0cs, ntm0)) / sqrt(x)
      theta = x - pi4 + csevl (z, bth0cs, ntth0) / x
      besy0 = ampl * sin (theta)
c
      return
      end
