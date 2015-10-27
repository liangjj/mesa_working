*deck daws
      function daws (x)
c***begin prologue  daws
c***purpose  compute dawson's function.
c***library   slatec (fnlib)
c***category  c8c
c***type      single precision (daws-s, ddaws-d)
c***keywords  dawson's function, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c daws(x) calculates dawson's integral for real argument x.
c
c series for daw        on the interval  0.          to  1.00000d+00
c                                        with weighted error   3.83e-17
c                                         log weighted error  16.42
c                               significant figures required  15.78
c                                    decimal places required  16.97
c
c series for daw2       on the interval  0.          to  1.60000d+01
c                                        with weighted error   5.17e-17
c                                         log weighted error  16.29
c                               significant figures required  15.90
c                                    decimal places required  17.02
c
c series for dawa       on the interval  0.          to  6.25000d-02
c                                        with weighted error   2.24e-17
c                                         log weighted error  16.65
c                               significant figures required  14.73
c                                    decimal places required  17.36
c
c***references  (none)
c***routines called  csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   780401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920618  removed space from variable names.  (rwc, wrb)
c***end prologue  daws
      dimension dawcs(13), daw2cs(29), dawacs(26)
      logical first
      save dawcs, daw2cs, dawacs, ntdaw, ntdaw2, ntdawa,
     1 xsml, xbig, xmax, first
      data dawcs( 1) /   -.0063517343 75145949e0 /
      data dawcs( 2) /   -.2294071479 6773869e0 /
      data dawcs( 3) /    .0221305009 39084764e0 /
      data dawcs( 4) /   -.0015492654 53892985e0 /
      data dawcs( 5) /    .0000849732 77156849e0 /
      data dawcs( 6) /   -.0000038282 66270972e0 /
      data dawcs( 7) /    .0000001462 85480625e0 /
      data dawcs( 8) /   -.0000000048 51982381e0 /
      data dawcs( 9) /    .0000000001 42146357e0 /
      data dawcs(10) /   -.0000000000 03728836e0 /
      data dawcs(11) /    .0000000000 00088549e0 /
      data dawcs(12) /   -.0000000000 00001920e0 /
      data dawcs(13) /    .0000000000 00000038e0 /
      data daw2cs( 1) /   -.0568865441 05215527e0 /
      data daw2cs( 2) /   -.3181134699 6168131e0 /
      data daw2cs( 3) /    .2087384541 3642237e0 /
      data daw2cs( 4) /   -.1247540991 3779131e0 /
      data daw2cs( 5) /    .0678693051 86676777e0 /
      data daw2cs( 6) /   -.0336591448 95270940e0 /
      data daw2cs( 7) /    .0152607812 71987972e0 /
      data daw2cs( 8) /   -.0063483709 62596214e0 /
      data daw2cs( 9) /    .0024326740 92074852e0 /
      data daw2cs(10) /   -.0008621954 14910650e0 /
      data daw2cs(11) /    .0002837657 33363216e0 /
      data daw2cs(12) /   -.0000870575 49874170e0 /
      data daw2cs(13) /    .0000249868 49985481e0 /
      data daw2cs(14) /   -.0000067319 28676416e0 /
      data daw2cs(15) /    .0000017078 57878557e0 /
      data daw2cs(16) /   -.0000004091 75512264e0 /
      data daw2cs(17) /    .0000000928 28292216e0 /
      data daw2cs(18) /   -.0000000199 91403610e0 /
      data daw2cs(19) /    .0000000040 96349064e0 /
      data daw2cs(20) /   -.0000000008 00324095e0 /
      data daw2cs(21) /    .0000000001 49385031e0 /
      data daw2cs(22) /   -.0000000000 26687999e0 /
      data daw2cs(23) /    .0000000000 04571221e0 /
      data daw2cs(24) /   -.0000000000 00751873e0 /
      data daw2cs(25) /    .0000000000 00118931e0 /
      data daw2cs(26) /   -.0000000000 00018116e0 /
      data daw2cs(27) /    .0000000000 00002661e0 /
      data daw2cs(28) /   -.0000000000 00000377e0 /
      data daw2cs(29) /    .0000000000 00000051e0 /
      data dawacs( 1) /    .0169048563 7765704e0 /
      data dawacs( 2) /    .0086832522 7840695e0 /
      data dawacs( 3) /    .0002424864 0424177e0 /
      data dawacs( 4) /    .0000126118 2399572e0 /
      data dawacs( 5) /    .0000010664 5331463e0 /
      data dawacs( 6) /    .0000001358 1597947e0 /
      data dawacs( 7) /    .0000000217 1042356e0 /
      data dawacs( 8) /    .0000000028 6701050e0 /
      data dawacs( 9) /   -.0000000001 9013363e0 /
      data dawacs(10) /   -.0000000003 0977804e0 /
      data dawacs(11) /   -.0000000001 0294148e0 /
      data dawacs(12) /   -.0000000000 0626035e0 /
      data dawacs(13) /    .0000000000 0856313e0 /
      data dawacs(14) /    .0000000000 0303304e0 /
      data dawacs(15) /   -.0000000000 0025236e0 /
      data dawacs(16) /   -.0000000000 0042106e0 /
      data dawacs(17) /   -.0000000000 0004431e0 /
      data dawacs(18) /    .0000000000 0004911e0 /
      data dawacs(19) /    .0000000000 0001235e0 /
      data dawacs(20) /   -.0000000000 0000578e0 /
      data dawacs(21) /   -.0000000000 0000228e0 /
      data dawacs(22) /    .0000000000 0000076e0 /
      data dawacs(23) /    .0000000000 0000038e0 /
      data dawacs(24) /   -.0000000000 0000011e0 /
      data dawacs(25) /   -.0000000000 0000006e0 /
      data dawacs(26) /    .0000000000 0000002e0 /
      data first /.true./
c***first executable statement  daws
      if (first) then
         eps = r1mach(3)
         ntdaw  = inits (dawcs,  13, 0.1*eps)
         ntdaw2 = inits (daw2cs, 29, 0.1*eps)
         ntdawa = inits (dawacs, 26, 0.1*eps)
c
         xsml = sqrt (1.5*eps)
         xbig = sqrt (0.5/eps)
         xmax = exp (min (-log(2.*r1mach(1)), log(r1mach(2))) - 1.0)
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.1.0) go to 20
c
      daws = x
      if (y.le.xsml) return
c
      daws = x * (0.75 + csevl (2.0*y*y-1.0, dawcs, ntdaw))
      return
c
 20   if (y.gt.4.0) go to 30
      daws = x * (0.25 + csevl (0.125*y*y-1.0, daw2cs, ntdaw2))
      return
c
 30   if (y.gt.xmax) go to 40
      daws = 0.5/x
      if (y.gt.xbig) return
c
      daws = (0.5 + csevl (32.0/y**2-1.0, dawacs, ntdawa)) / x
      return
c
 40   call xermsg ('slatec', 'daws', 'abs(x) so large daws underflows',
     +   1, 1)
      daws = 0.0
      return
c
      end
