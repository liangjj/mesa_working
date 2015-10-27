*deck dbesy1
      double precision function dbesy1 (x)
c***begin prologue  dbesy1
c***purpose  compute the bessel function of the second kind of order
c            one.
c***library   slatec (fnlib)
c***category  c10a1
c***type      double precision (besy1-s, dbesy1-d)
c***keywords  bessel function, fnlib, order one, second kind,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c dbesy1(x) calculates the double precision bessel function of the
c second kind of order for double precision argument x.
c
c series for by1        on the interval  0.          to  1.60000e+01
c                                        with weighted error   8.65e-33
c                                         log weighted error  32.06
c                               significant figures required  32.17
c                                    decimal places required  32.71
c
c***references  (none)
c***routines called  d1mach, d9b1mp, dbesj1, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dbesy1
      double precision x, by1cs(20), ampl, theta, twodpi, xmin, xsml,
     1  y, d1mach, dcsevl, dbesj1
      logical first
      save by1cs, twodpi, nty1, xmin, xsml, first
      data by1cs(  1) / +.3208047100 6119086293 2352018628 015 d-1    /
      data by1cs(  2) / +.1262707897 4335004495 3431725999 727 d+1    /
      data by1cs(  3) / +.6499961899 9231750009 7490637314 144 d-2    /
      data by1cs(  4) / -.8936164528 8605041165 3144160009 712 d-1    /
      data by1cs(  5) / +.1325088122 1757095451 2375510370 043 d-1    /
      data by1cs(  6) / -.8979059119 6483523775 3039508298 105 d-3    /
      data by1cs(  7) / +.3647361487 9583067824 2287368165 349 d-4    /
      data by1cs(  8) / -.1001374381 6660005554 9075523845 295 d-5    /
      data by1cs(  9) / +.1994539657 3901739703 1159372421 243 d-7    /
      data by1cs( 10) / -.3023065601 8033816728 4799332520 743 d-9    /
      data by1cs( 11) / +.3609878156 9478119611 6252914242 474 d-11   /
      data by1cs( 12) / -.3487488297 2875824241 4552947409 066 d-13   /
      data by1cs( 13) / +.2783878971 5591766581 3507698517 333 d-15   /
      data by1cs( 14) / -.1867870968 6194876876 6825352533 333 d-17   /
      data by1cs( 15) / +.1068531533 9116825975 7070336000 000 d-19   /
      data by1cs( 16) / -.5274721956 6844822894 3872000000 000 d-22   /
      data by1cs( 17) / +.2270199403 1556641437 0133333333 333 d-24   /
      data by1cs( 18) / -.8595390353 9452310869 3333333333 333 d-27   /
      data by1cs( 19) / +.2885404379 8337945600 0000000000 000 d-29   /
      data by1cs( 20) / -.8647541138 9371733333 3333333333 333 d-32   /
      data twodpi / 0.6366197723 6758134307 5535053490 057 d0 /
      data first /.true./
c***first executable statement  dbesy1
      if (first) then
         nty1 = initds (by1cs, 20, 0.1*real(d1mach(3)))
c
         xmin = 1.571d0 * exp (max(log(d1mach(1)), -log(d1mach(2))) +
     1     0.01d0)
         xsml = sqrt(4.0d0*d1mach(3))
      endif
      first = .false.
c
      if (x .le. 0.d0) call xermsg ('slatec', 'dbesy1',
     +   'x is zero or negative', 1, 2)
      if (x.gt.4.0d0) go to 20
c
      if (x .lt. xmin) call xermsg ('slatec', 'dbesy1',
     +   'x so small y1 overflows', 3, 2)
      y = 0.d0
      if (x.gt.xsml) y = x*x
      dbesy1 = twodpi * log(0.5d0*x)*dbesj1(x) + (0.5d0 +
     1  dcsevl (.125d0*y-1.d0, by1cs, nty1))/x
      return
c
 20   call d9b1mp (x, ampl, theta)
      dbesy1 = ampl * sin(theta)
      return
c
      end
