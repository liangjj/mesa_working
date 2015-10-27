*deck dbesk0
      double precision function dbesk0 (x)
c***begin prologue  dbesk0
c***purpose  compute the modified (hyperbolic) bessel function of the
c            third kind of order zero.
c***library   slatec (fnlib)
c***category  c10b1
c***type      double precision (besk0-s, dbesk0-d)
c***keywords  fnlib, hyperbolic bessel function,
c             modified bessel function, order zero, special functions,
c             third kind
c***author  fullerton, w., (lanl)
c***description
c
c dbesk0(x) calculates the double precision modified (hyperbolic)
c bessel function of the third kind of order zero for double
c precision argument x.  the argument must be greater than zero
c but not so large that the result underflows.
c
c series for bk0        on the interval  0.          to  4.00000e+00
c                                        with weighted error   3.08e-33
c                                         log weighted error  32.51
c                               significant figures required  32.05
c                                    decimal places required  33.11
c
c***references  (none)
c***routines called  d1mach, dbesi0, dbsk0e, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dbesk0
      double precision x, bk0cs(16), xmax, xmaxt, xsml, y,
     1  d1mach, dcsevl, dbesi0, dbsk0e
      logical first
      save bk0cs, ntk0, xsml, xmax, first
      data bk0cs(  1) / -.3532739323 3902768720 1140060063 153 d-1    /
      data bk0cs(  2) / +.3442898999 2462848688 6344927529 213 d+0    /
      data bk0cs(  3) / +.3597993651 5361501626 5721303687 231 d-1    /
      data bk0cs(  4) / +.1264615411 4469259233 8479508673 447 d-2    /
      data bk0cs(  5) / +.2286212103 1194517860 8269830297 585 d-4    /
      data bk0cs(  6) / +.2534791079 0261494573 0790013428 354 d-6    /
      data bk0cs(  7) / +.1904516377 2202088589 7214059381 366 d-8    /
      data bk0cs(  8) / +.1034969525 7633624585 1008317853 089 d-10   /
      data bk0cs(  9) / +.4259816142 7910825765 2445327170 133 d-13   /
      data bk0cs( 10) / +.1374465435 8807508969 4238325440 000 d-15   /
      data bk0cs( 11) / +.3570896528 5083735909 9688597333 333 d-18   /
      data bk0cs( 12) / +.7631643660 1164373766 7498666666 666 d-21   /
      data bk0cs( 13) / +.1365424988 4407818590 8053333333 333 d-23   /
      data bk0cs( 14) / +.2075275266 9066680831 9999999999 999 d-26   /
      data bk0cs( 15) / +.2712814218 0729856000 0000000000 000 d-29   /
      data bk0cs( 16) / +.3082593887 9146666666 6666666666 666 d-32   /
      data first /.true./
c***first executable statement  dbesk0
      if (first) then
         ntk0 = initds (bk0cs, 16, 0.1*real(d1mach(3)))
         xsml = sqrt(4.0d0*d1mach(3))
         xmaxt = -log(d1mach(1))
         xmax = xmaxt - 0.5d0*xmaxt*log(xmaxt)/(xmaxt+0.5d0)
      endif
      first = .false.
c
      if (x .le. 0.d0) call xermsg ('slatec', 'dbesk0',
     +   'x is zero or negative', 2, 2)
      if (x.gt.2.0d0) go to 20
c
      y = 0.d0
      if (x.gt.xsml) y = x*x
      dbesk0 = -log(0.5d0*x)*dbesi0(x) - 0.25d0 + dcsevl (.5d0*y-1.d0,
     1  bk0cs, ntk0)
      return
c
 20   dbesk0 = 0.d0
      if (x .gt. xmax) call xermsg ('slatec', 'dbesk0',
     +   'x so big k0 underflows', 1, 1)
      if (x.gt.xmax) return
c
      dbesk0 = exp(-x) * dbsk0e(x)
c
      return
      end
