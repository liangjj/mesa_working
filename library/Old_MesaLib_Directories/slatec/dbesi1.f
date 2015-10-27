*deck dbesi1
      double precision function dbesi1 (x)
c***begin prologue  dbesi1
c***purpose  compute the modified (hyperbolic) bessel function of the
c            first kind of order one.
c***library   slatec (fnlib)
c***category  c10b1
c***type      double precision (besi1-s, dbesi1-d)
c***keywords  first kind, fnlib, hyperbolic bessel function,
c             modified bessel function, order one, special functions
c***author  fullerton, w., (lanl)
c***description
c
c dbesi1(x) calculates the double precision modified (hyperbolic)
c bessel function of the first kind of order one and double precision
c argument x.
c
c series for bi1        on the interval  0.          to  9.00000e+00
c                                        with weighted error   1.44e-32
c                                         log weighted error  31.84
c                               significant figures required  31.45
c                                    decimal places required  32.46
c
c***references  (none)
c***routines called  d1mach, dbsi1e, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dbesi1
      double precision x, bi1cs(17), xmax, xmin, xsml, y, d1mach,
     1  dcsevl, dbsi1e
      logical first
      save bi1cs, nti1, xmin, xsml, xmax, first
      data bi1cs(  1) / -.1971713261 0998597316 1385032181 49 d-2     /
      data bi1cs(  2) / +.4073488766 7546480608 1553936520 14 d+0     /
      data bi1cs(  3) / +.3483899429 9959455866 2450377837 87 d-1     /
      data bi1cs(  4) / +.1545394556 3001236038 5984010584 89 d-2     /
      data bi1cs(  5) / +.4188852109 8377784129 4588320041 20 d-4     /
      data bi1cs(  6) / +.7649026764 8362114741 9597039660 69 d-6     /
      data bi1cs(  7) / +.1004249392 4741178689 1798080372 38 d-7     /
      data bi1cs(  8) / +.9932207791 9238106481 3712980548 63 d-10    /
      data bi1cs(  9) / +.7663801791 8447637275 2001716813 49 d-12    /
      data bi1cs( 10) / +.4741418923 8167394980 3880919481 60 d-14    /
      data bi1cs( 11) / +.2404114404 0745181799 8631720320 00 d-16    /
      data bi1cs( 12) / +.1017150500 7093713649 1211007999 99 d-18    /
      data bi1cs( 13) / +.3645093565 7866949458 4917333333 33 d-21    /
      data bi1cs( 14) / +.1120574950 2562039344 8106666666 66 d-23    /
      data bi1cs( 15) / +.2987544193 4468088832 0000000000 00 d-26    /
      data bi1cs( 16) / +.6973231093 9194709333 3333333333 33 d-29    /
      data bi1cs( 17) / +.1436794822 0620800000 0000000000 00 d-31    /
      data first /.true./
c***first executable statement  dbesi1
      if (first) then
         nti1 = initds (bi1cs, 17, 0.1*real(d1mach(3)))
         xmin = 2.0d0*d1mach(1)
         xsml = sqrt(4.5d0*d1mach(3))
         xmax = log (d1mach(2))
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.3.0d0) go to 20
c
      dbesi1 = 0.d0
      if (y.eq.0.d0)  return
c
      if (y .le. xmin) call xermsg ('slatec', 'dbesi1',
     +   'abs(x) so small i1 underflows', 1, 1)
      if (y.gt.xmin) dbesi1 = 0.5d0*x
      if (y.gt.xsml) dbesi1 = x*(0.875d0 + dcsevl (y*y/4.5d0-1.d0,
     1  bi1cs, nti1))
      return
c
 20   if (y .gt. xmax) call xermsg ('slatec', 'dbesi1',
     +   'abs(x) so big i1 overflows', 2, 2)
c
      dbesi1 = exp(y) * dbsi1e(x)
c
      return
      end
