*deck r9ln2r
      function r9ln2r (x)
c***begin prologue  r9ln2r
c***subsidiary
c***purpose  evaluate log(1+x) from second order relative accuracy so
c            that log(1+x) = x - x**2/2 + x**3*r9ln2r(x).
c***library   slatec (fnlib)
c***category  c4b
c***type      single precision (r9ln2r-s, d9ln2r-d, c9ln2r-c)
c***keywords  elementary functions, fnlib, logarithm, second order
c***author  fullerton, w., (lanl)
c***description
c
c evaluate  log(1+x)  from 2-nd order with relative error accuracy so
c that    log(1+x) = x - x**2/2 + x**3*r9ln2r(x)
c
c series for ln21       on the interval -6.25000d-01 to  0.
c                                        with weighted error   2.49e-17
c                                         log weighted error  16.60
c                               significant figures required  15.87
c                                    decimal places required  17.31
c
c series for ln22       on the interval  0.          to  8.12500d-01
c                                        with weighted error   1.42e-17
c                                         log weighted error  16.85
c                               significant figures required  15.95
c                                    decimal places required  17.50
c
c***references  (none)
c***routines called  csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   780401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900720  routine changed from user-callable to subsidiary.  (wrb)
c***end prologue  r9ln2r
      real ln21cs(26), ln22cs(20)
      logical first
      save ln21cs, ln22cs, ntln21, ntln22, xmin, xbig, xmax, first
      data ln21cs( 1) /    .1811196251 3478810e0 /
      data ln21cs( 2) /   -.1562712319 2872463e0 /
      data ln21cs( 3) /    .0286763053 61557275e0 /
      data ln21cs( 4) /   -.0055586996 55948139e0 /
      data ln21cs( 5) /    .0011178976 65229983e0 /
      data ln21cs( 6) /   -.0002308050 89823279e0 /
      data ln21cs( 7) /    .0000485988 53341100e0 /
      data ln21cs( 8) /   -.0000103901 27388903e0 /
      data ln21cs( 9) /    .0000022484 56370739e0 /
      data ln21cs(10) /   -.0000004914 05927392e0 /
      data ln21cs(11) /    .0000001082 82565070e0 /
      data ln21cs(12) /   -.0000000240 25872763e0 /
      data ln21cs(13) /    .0000000053 62460047e0 /
      data ln21cs(14) /   -.0000000012 02995136e0 /
      data ln21cs(15) /    .0000000002 71078892e0 /
      data ln21cs(16) /   -.0000000000 61323562e0 /
      data ln21cs(17) /    .0000000000 13920858e0 /
      data ln21cs(18) /   -.0000000000 03169930e0 /
      data ln21cs(19) /    .0000000000 00723837e0 /
      data ln21cs(20) /   -.0000000000 00165700e0 /
      data ln21cs(21) /    .0000000000 00038018e0 /
      data ln21cs(22) /   -.0000000000 00008741e0 /
      data ln21cs(23) /    .0000000000 00002013e0 /
      data ln21cs(24) /   -.0000000000 00000464e0 /
      data ln21cs(25) /    .0000000000 00000107e0 /
      data ln21cs(26) /   -.0000000000 00000024e0 /
      data ln22cs( 1) /   -.2224253253 5020461e0 /
      data ln22cs( 2) /   -.0610471001 08078624e0 /
      data ln22cs( 3) /    .0074272350 09750394e0 /
      data ln22cs( 4) /   -.0009335018 26163697e0 /
      data ln22cs( 5) /    .0001200499 07687260e0 /
      data ln22cs( 6) /   -.0000157047 22952820e0 /
      data ln22cs( 7) /    .0000020818 74781051e0 /
      data ln22cs( 8) /   -.0000002789 19557764e0 /
      data ln22cs( 9) /    .0000000376 93558237e0 /
      data ln22cs(10) /   -.0000000051 30902896e0 /
      data ln22cs(11) /    .0000000007 02714117e0 /
      data ln22cs(12) /   -.0000000000 96748595e0 /
      data ln22cs(13) /    .0000000000 13381046e0 /
      data ln22cs(14) /   -.0000000000 01858102e0 /
      data ln22cs(15) /    .0000000000 00258929e0 /
      data ln22cs(16) /   -.0000000000 00036195e0 /
      data ln22cs(17) /    .0000000000 00005074e0 /
      data ln22cs(18) /   -.0000000000 00000713e0 /
      data ln22cs(19) /    .0000000000 00000100e0 /
      data ln22cs(20) /   -.0000000000 00000014e0 /
      data first /.true./
c***first executable statement  r9ln2r
      if (first) then
         eps = r1mach(3)
         ntln21 = inits (ln21cs, 26, 0.1*eps)
         ntln22 = inits (ln22cs, 20, 0.1*eps)
c
         xmin = -1.0 + sqrt(r1mach(4))
         sqeps = sqrt(eps)
         txmax = 6.0/sqeps
         xmax = txmax - (eps*txmax**2 - 2.0*log(txmax)) /
     1                                              (2.0*eps*txmax)
         txbig = 4.0/sqrt(sqeps)
         xbig = txbig - (sqeps*txbig**2 - 2.0*log(txbig)) /
     1                                                (2.*sqeps*txbig)
      endif
      first = .false.
c
      if (x.lt.(-0.625) .or. x.gt.0.8125) go to 20
c
      if (x.lt.0.0) r9ln2r = 0.375 + csevl (16.*x/5.+1.0, ln21cs,
     1  ntln21)
      if (x.ge.0.0) r9ln2r = 0.375 + csevl (32.*x/13.-1.0, ln22cs,
     1  ntln22)
      return
c
 20   if (x .lt. xmin) call xermsg ('slatec', 'r9ln2r',
     +   'answer lt half precision because x is too near -1', 1, 1)
      if (x .gt. xmax) call xermsg ('slatec', 'r9ln2r',
     +   'no precision in answer because x is too big', 3, 2)
      if (x .gt. xbig) call xermsg ('slatec', 'r9ln2r',
     +   'answer lt half precision because x is too big', 2, 1)
c
      r9ln2r = (log(1.0+x) - x*(1.0-0.5*x) ) / x**3
      return
c
      end
