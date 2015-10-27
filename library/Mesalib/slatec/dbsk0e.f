*deck dbsk0e
      double precision function dbsk0e (x)
c***begin prologue  dbsk0e
c***purpose  compute the exponentially scaled modified (hyperbolic)
c            bessel function of the third kind of order zero.
c***library   slatec (fnlib)
c***category  c10b1
c***type      double precision (besk0e-s, dbsk0e-d)
c***keywords  exponentially scaled, fnlib, hyperbolic bessel function,
c             modified bessel function, order zero, special functions,
c             third kind
c***author  fullerton, w., (lanl)
c***description
c
c dbsk0e(x) computes the double precision exponentially scaled
c modified (hyperbolic) bessel function of the third kind of
c order zero for positive double precision argument x.
c
c series for bk0        on the interval  0.          to  4.00000e+00
c                                        with weighted error   3.08e-33
c                                         log weighted error  32.51
c                               significant figures required  32.05
c                                    decimal places required  33.11
c
c series for ak0        on the interval  1.25000e-01 to  5.00000e-01
c                                        with weighted error   2.85e-32
c                                         log weighted error  31.54
c                               significant figures required  30.19
c                                    decimal places required  32.33
c
c series for ak02       on the interval  0.          to  1.25000e-01
c                                        with weighted error   2.30e-32
c                                         log weighted error  31.64
c                               significant figures required  29.68
c                                    decimal places required  32.40
c
c***references  (none)
c***routines called  d1mach, dbesi0, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dbsk0e
      double precision x, bk0cs(16), ak0cs(38), ak02cs(33),
     1  xsml, y, d1mach, dcsevl, dbesi0
      logical first
      save bk0cs, ak0cs, ak02cs, ntk0, ntak0, ntak02, xsml, first
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
      data ak0cs(  1) / -.7643947903 3279414240 8297827008 8 d-1      /
      data ak0cs(  2) / -.2235652605 6998190520 2309555079 1 d-1      /
      data ak0cs(  3) / +.7734181154 6938582353 0061817404 7 d-3      /
      data ak0cs(  4) / -.4281006688 8860994644 5214643541 6 d-4      /
      data ak0cs(  5) / +.3081700173 8629747436 5001482666 0 d-5      /
      data ak0cs(  6) / -.2639367222 0096649740 6744889272 3 d-6      /
      data ak0cs(  7) / +.2563713036 4034692062 9408826574 2 d-7      /
      data ak0cs(  8) / -.2742705549 9002012638 5721191524 4 d-8      /
      data ak0cs(  9) / +.3169429658 0974995920 8083287340 3 d-9      /
      data ak0cs( 10) / -.3902353286 9621841416 0106571796 2 d-10     /
      data ak0cs( 11) / +.5068040698 1885754020 5009212728 6 d-11     /
      data ak0cs( 12) / -.6889574741 0078706795 4171355798 4 d-12     /
      data ak0cs( 13) / +.9744978497 8259176913 8820133683 1 d-13     /
      data ak0cs( 14) / -.1427332841 8845485053 8985534012 2 d-13     /
      data ak0cs( 15) / +.2156412571 0214630395 5806297652 7 d-14     /
      data ak0cs( 16) / -.3349654255 1495627721 8878205853 0 d-15     /
      data ak0cs( 17) / +.5335260216 9529116921 4528039260 1 d-16     /
      data ak0cs( 18) / -.8693669980 8907538076 3962237883 7 d-17     /
      data ak0cs( 19) / +.1446404347 8622122278 8776344234 6 d-17     /
      data ak0cs( 20) / -.2452889825 5001296824 0467875157 3 d-18     /
      data ak0cs( 21) / +.4233754526 2321715728 2170634240 0 d-19     /
      data ak0cs( 22) / -.7427946526 4544641956 9534129493 3 d-20     /
      data ak0cs( 23) / +.1323150529 3926668662 7796746240 0 d-20     /
      data ak0cs( 24) / -.2390587164 7396494513 3598146559 9 d-21     /
      data ak0cs( 25) / +.4376827585 9232261401 6571255466 6 d-22     /
      data ak0cs( 26) / -.8113700607 3451180593 3901141333 3 d-23     /
      data ak0cs( 27) / +.1521819913 8321729583 1037815466 6 d-23     /
      data ak0cs( 28) / -.2886041941 4833977702 3595861333 3 d-24     /
      data ak0cs( 29) / +.5530620667 0547179799 9261013333 3 d-25     /
      data ak0cs( 30) / -.1070377329 2498987285 9163306666 6 d-25     /
      data ak0cs( 31) / +.2091086893 1423843002 9632853333 3 d-26     /
      data ak0cs( 32) / -.4121713723 6462038274 1026133333 3 d-27     /
      data ak0cs( 33) / +.8193483971 1213076401 3568000000 0 d-28     /
      data ak0cs( 34) / -.1642000275 4592977267 8075733333 3 d-28     /
      data ak0cs( 35) / +.3316143281 4802271958 9034666666 6 d-29     /
      data ak0cs( 36) / -.6746863644 1452959410 8586666666 6 d-30     /
      data ak0cs( 37) / +.1382429146 3184246776 3541333333 3 d-30     /
      data ak0cs( 38) / -.2851874167 3598325708 1173333333 3 d-31     /
      data ak02cs(  1) / -.1201869826 3075922398 3934621245 2 d-1      /
      data ak02cs(  2) / -.9174852691 0256953106 5256107571 3 d-2      /
      data ak02cs(  3) / +.1444550931 7750058210 4884387805 7 d-3      /
      data ak02cs(  4) / -.4013614175 4357097286 7102107787 9 d-5      /
      data ak02cs(  5) / +.1567831810 8523106725 9034899033 3 d-6      /
      data ak02cs(  6) / -.7770110438 5217377103 1579975446 0 d-8      /
      data ak02cs(  7) / +.4611182576 1797178825 3313052958 6 d-9      /
      data ak02cs(  8) / -.3158592997 8605657705 2666580330 9 d-10     /
      data ak02cs(  9) / +.2435018039 3650411278 3588781432 9 d-11     /
      data ak02cs( 10) / -.2074331387 3983478977 0985337350 6 d-12     /
      data ak02cs( 11) / +.1925787280 5899170847 4273650469 3 d-13     /
      data ak02cs( 12) / -.1927554805 8389561036 0034718221 8 d-14     /
      data ak02cs( 13) / +.2062198029 1978182782 8523786964 4 d-15     /
      data ak02cs( 14) / -.2341685117 5792424026 0364019507 1 d-16     /
      data ak02cs( 15) / +.2805902810 6430422468 1517882845 8 d-17     /
      data ak02cs( 16) / -.3530507631 1618079458 1548246357 3 d-18     /
      data ak02cs( 17) / +.4645295422 9351082674 2421633706 6 d-19     /
      data ak02cs( 18) / -.6368625941 3442664739 2205346133 3 d-20     /
      data ak02cs( 19) / +.9069521310 9865155676 2234880000 0 d-21     /
      data ak02cs( 20) / -.1337974785 4236907398 4500531199 9 d-21     /
      data ak02cs( 21) / +.2039836021 8599523155 2208896000 0 d-22     /
      data ak02cs( 22) / -.3207027481 3678405000 6086997333 3 d-23     /
      data ak02cs( 23) / +.5189744413 6623099636 2635946666 6 d-24     /
      data ak02cs( 24) / -.8629501497 5405721929 6460799999 9 d-25     /
      data ak02cs( 25) / +.1472161183 1025598552 0803840000 0 d-25     /
      data ak02cs( 26) / -.2573069023 8670112838 1235199999 9 d-26     /
      data ak02cs( 27) / +.4601774086 6435165873 7664000000 0 d-27     /
      data ak02cs( 28) / -.8411555324 2010937371 3066666666 6 d-28     /
      data ak02cs( 29) / +.1569806306 6353689393 0154666666 6 d-28     /
      data ak02cs( 30) / -.2988226453 0057577889 7919999999 9 d-29     /
      data ak02cs( 31) / +.5796831375 2168365206 1866666666 6 d-30     /
      data ak02cs( 32) / -.1145035994 3476813321 5573333333 3 d-30     /
      data ak02cs( 33) / +.2301266594 2496828020 0533333333 3 d-31     /
      data first /.true./
c***first executable statement  dbsk0e
      if (first) then
         eta = 0.1*real(d1mach(3))
         ntk0 = initds (bk0cs, 16, eta)
         ntak0 = initds (ak0cs, 38, eta)
         ntak02 = initds (ak02cs, 33, eta)
         xsml = sqrt(4.0d0*d1mach(3))
      endif
      first = .false.
c
      if (x .le. 0.d0) call xermsg ('slatec', 'dbsk0e',
     +   'x is zero or negative', 2, 2)
      if (x.gt.2.0d0) go to 20
c
      y = 0.d0
      if (x.gt.xsml) y = x*x
      dbsk0e = exp(x)*(-log(0.5d0*x)*dbesi0(x) - 0.25d0 +
     1  dcsevl (.5d0*y-1.d0, bk0cs, ntk0))
      return
c
 20   if (x.le.8.d0) dbsk0e = (1.25d0 + dcsevl ((16.d0/x-5.d0)/3.d0,
     1  ak0cs, ntak0))/sqrt(x)
      if (x.gt.8.d0) dbsk0e = (1.25d0 +
     1  dcsevl (16.d0/x-1.d0, ak02cs, ntak02))/sqrt(x)
c
      return
      end
