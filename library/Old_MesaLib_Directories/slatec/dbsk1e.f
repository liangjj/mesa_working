*deck dbsk1e
      double precision function dbsk1e (x)
c***begin prologue  dbsk1e
c***purpose  compute the exponentially scaled modified (hyperbolic)
c            bessel function of the third kind of order one.
c***library   slatec (fnlib)
c***category  c10b1
c***type      double precision (besk1e-s, dbsk1e-d)
c***keywords  exponentially scaled, fnlib, hyperbolic bessel function,
c             modified bessel function, order one, special functions,
c             third kind
c***author  fullerton, w., (lanl)
c***description
c
c dbsk1e(s) computes the double precision exponentially scaled
c modified (hyperbolic) bessel function of the third kind of order
c one for positive double precision argument x.
c
c series for bk1        on the interval  0.          to  4.00000e+00
c                                        with weighted error   9.16e-32
c                                         log weighted error  31.04
c                               significant figures required  30.61
c                                    decimal places required  31.64
c
c series for ak1        on the interval  1.25000e-01 to  5.00000e-01
c                                        with weighted error   3.07e-32
c                                         log weighted error  31.51
c                               significant figures required  30.71
c                                    decimal places required  32.30
c
c series for ak12       on the interval  0.          to  1.25000e-01
c                                        with weighted error   2.41e-32
c                                         log weighted error  31.62
c                               significant figures required  30.25
c                                    decimal places required  32.38
c
c***references  (none)
c***routines called  d1mach, dbesi1, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dbsk1e
      double precision x, bk1cs(16), ak1cs(38), ak12cs(33), xmin,
     1  xsml, y, d1mach, dcsevl, dbesi1
      logical first
      save bk1cs, ak1cs, ak12cs, ntk1, ntak1, ntak12, xmin, xsml,
     1  first
      data bk1cs(  1) / +.2530022733 8947770532 5311208685 33 d-1     /
      data bk1cs(  2) / -.3531559607 7654487566 7238316918 01 d+0     /
      data bk1cs(  3) / -.1226111808 2265714823 4790679300 42 d+0     /
      data bk1cs(  4) / -.6975723859 6398643501 8129202960 83 d-2     /
      data bk1cs(  5) / -.1730288957 5130520630 1765073689 79 d-3     /
      data bk1cs(  6) / -.2433406141 5659682349 6007350301 64 d-5     /
      data bk1cs(  7) / -.2213387630 7347258558 3152525451 26 d-7     /
      data bk1cs(  8) / -.1411488392 6335277610 9583302126 08 d-9     /
      data bk1cs(  9) / -.6666901694 1993290060 8537512643 73 d-12    /
      data bk1cs( 10) / -.2427449850 5193659339 2631968648 53 d-14    /
      data bk1cs( 11) / -.7023863479 3862875971 7837971200 00 d-17    /
      data bk1cs( 12) / -.1654327515 5100994675 4910293333 33 d-19    /
      data bk1cs( 13) / -.3233834745 9944491991 8933333333 33 d-22    /
      data bk1cs( 14) / -.5331275052 9265274999 4666666666 66 d-25    /
      data bk1cs( 15) / -.7513040716 2157226666 6666666666 66 d-28    /
      data bk1cs( 16) / -.9155085717 6541866666 6666666666 66 d-31    /
      data ak1cs(  1) / +.2744313406 9738829695 2576662272 66 d+0     /
      data ak1cs(  2) / +.7571989953 1993678170 8923781492 90 d-1     /
      data ak1cs(  3) / -.1441051556 4754061229 8531161756 25 d-2     /
      data ak1cs(  4) / +.6650116955 1257479394 2513854770 36 d-4     /
      data ak1cs(  5) / -.4369984709 5201407660 5808450891 67 d-5     /
      data ak1cs(  6) / +.3540277499 7630526799 4171390085 34 d-6     /
      data ak1cs(  7) / -.3311163779 2932920208 9826882457 04 d-7     /
      data ak1cs(  8) / +.3445977581 9010534532 3114997709 92 d-8     /
      data ak1cs(  9) / -.3898932347 4754271048 9819374927 58 d-9     /
      data ak1cs( 10) / +.4720819750 4658356400 9474493390 05 d-10    /
      data ak1cs( 11) / -.6047835662 8753562345 3735915628 90 d-11    /
      data ak1cs( 12) / +.8128494874 8658747888 1938379856 63 d-12    /
      data ak1cs( 13) / -.1138694574 7147891428 9239159510 42 d-12    /
      data ak1cs( 14) / +.1654035840 8462282325 9729482050 90 d-13    /
      data ak1cs( 15) / -.2480902567 7068848221 5160104405 33 d-14    /
      data ak1cs( 16) / +.3829237890 7024096948 4292272991 57 d-15    /
      data ak1cs( 17) / -.6064734104 0012418187 7682103773 86 d-16    /
      data ak1cs( 18) / +.9832425623 2648616038 1940046506 66 d-17    /
      data ak1cs( 19) / -.1628416873 8284380035 6666201156 26 d-17    /
      data ak1cs( 20) / +.2750153649 6752623718 2841203370 66 d-18    /
      data ak1cs( 21) / -.4728966646 3953250924 2810695680 00 d-19    /
      data ak1cs( 22) / +.8268150002 8109932722 3920503466 66 d-20    /
      data ak1cs( 23) / -.1468140513 6624956337 1939648853 33 d-20    /
      data ak1cs( 24) / +.2644763926 9208245978 0858948266 66 d-21    /
      data ak1cs( 25) / -.4829015756 4856387897 9698688000 00 d-22    /
      data ak1cs( 26) / +.8929302074 3610130180 6563327999 99 d-23    /
      data ak1cs( 27) / -.1670839716 8972517176 9977514666 66 d-23    /
      data ak1cs( 28) / +.3161645603 4040694931 3686186666 66 d-24    /
      data ak1cs( 29) / -.6046205531 2274989106 5064106666 66 d-25    /
      data ak1cs( 30) / +.1167879894 2042732700 7184213333 33 d-25    /
      data ak1cs( 31) / -.2277374158 2653996232 8678400000 00 d-26    /
      data ak1cs( 32) / +.4481109730 0773675795 3058133333 33 d-27    /
      data ak1cs( 33) / -.8893288476 9020194062 3360000000 00 d-28    /
      data ak1cs( 34) / +.1779468001 8850275131 3920000000 00 d-28    /
      data ak1cs( 35) / -.3588455596 7329095821 9946666666 66 d-29    /
      data ak1cs( 36) / +.7290629049 2694257991 6799999999 99 d-30    /
      data ak1cs( 37) / -.1491844984 5546227073 0240000000 00 d-30    /
      data ak1cs( 38) / +.3073657387 2934276300 7999999999 99 d-31    /
      data ak12cs(  1) / +.6379308343 7390010366 0048853410 2 d-1      /
      data ak12cs(  2) / +.2832887813 0497209358 3503028470 8 d-1      /
      data ak12cs(  3) / -.2475370673 9052503454 1454556673 2 d-3      /
      data ak12cs(  4) / +.5771972451 6072488204 7097662576 3 d-5      /
      data ak12cs(  5) / -.2068939219 5365483027 4553319655 2 d-6      /
      data ak12cs(  6) / +.9739983441 3818041803 0921309788 7 d-8      /
      data ak12cs(  7) / -.5585336140 3806249846 8889551112 9 d-9      /
      data ak12cs(  8) / +.3732996634 0461852402 2121285473 1 d-10     /
      data ak12cs(  9) / -.2825051961 0232254451 3506575492 8 d-11     /
      data ak12cs( 10) / +.2372019002 4841441736 4349695548 6 d-12     /
      data ak12cs( 11) / -.2176677387 9917539792 6830166793 8 d-13     /
      data ak12cs( 12) / +.2157914161 6160324539 3956268970 6 d-14     /
      data ak12cs( 13) / -.2290196930 7182692759 9155133815 4 d-15     /
      data ak12cs( 14) / +.2582885729 8232749619 1993956522 6 d-16     /
      data ak12cs( 15) / -.3076752641 2684631876 2109817344 0 d-17     /
      data ak12cs( 16) / +.3851487721 2804915970 9489684479 9 d-18     /
      data ak12cs( 17) / -.5044794897 6415289771 1728250880 0 d-19     /
      data ak12cs( 18) / +.6888673850 4185442370 1829222399 9 d-20     /
      data ak12cs( 19) / -.9775041541 9501183030 0213248000 0 d-21     /
      data ak12cs( 20) / +.1437416218 5238364610 0165973333 3 d-21     /
      data ak12cs( 21) / -.2185059497 3443473734 9973333333 3 d-22     /
      data ak12cs( 22) / +.3426245621 8092206316 4538880000 0 d-23     /
      data ak12cs( 23) / -.5531064394 2464082325 0124800000 0 d-24     /
      data ak12cs( 24) / +.9176601505 6859954037 8282666666 6 d-25     /
      data ak12cs( 25) / -.1562287203 6180249114 4874666666 6 d-25     /
      data ak12cs( 26) / +.2725419375 4843331323 4943999999 9 d-26     /
      data ak12cs( 27) / -.4865674910 0748279923 7802666666 6 d-27     /
      data ak12cs( 28) / +.8879388552 7235025873 5786666666 6 d-28     /
      data ak12cs( 29) / -.1654585918 0392575489 3653333333 3 d-28     /
      data ak12cs( 30) / +.3145111321 3578486743 0399999999 9 d-29     /
      data ak12cs( 31) / -.6092998312 1931276124 1600000000 0 d-30     /
      data ak12cs( 32) / +.1202021939 3698158346 2399999999 9 d-30     /
      data ak12cs( 33) / -.2412930801 4594088413 8666666666 6 d-31     /
      data first /.true./
c***first executable statement  dbsk1e
      if (first) then
         eta = 0.1*real(d1mach(3))
         ntk1 = initds (bk1cs, 16, eta)
         ntak1 = initds (ak1cs, 38, eta)
         ntak12 = initds (ak12cs, 33, eta)
c
         xmin = exp (max(log(d1mach(1)), -log(d1mach(2))) + 0.01d0)
         xsml = sqrt(4.0d0*d1mach(3))
      endif
      first = .false.
c
      if (x .le. 0.d0) call xermsg ('slatec', 'dbsk1e',
     +   'x is zero or negative', 2, 2)
      if (x.gt.2.0d0) go to 20
c
      if (x .lt. xmin) call xermsg ('slatec', 'dbsk1e',
     +   'x so small k1 overflows', 3, 2)
      y = 0.d0
      if (x.gt.xsml) y = x*x
      dbsk1e = exp(x)*(log(0.5d0*x)*dbesi1(x) + (0.75d0 +
     1  dcsevl (0.5d0*y-1.d0, bk1cs, ntk1))/x )
      return
c
 20   if (x.le.8.d0) dbsk1e = (1.25d0 + dcsevl ((16.d0/x-5.d0)/3.d0,
     1  ak1cs, ntak1))/sqrt(x)
      if (x.gt.8.d0) dbsk1e = (1.25d0 +
     1  dcsevl (16.d0/x-1.d0, ak12cs, ntak12))/sqrt(x)
c
      return
      end
