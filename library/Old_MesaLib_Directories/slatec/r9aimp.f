*deck r9aimp
      subroutine r9aimp (x, ampl, theta)
c***begin prologue  r9aimp
c***subsidiary
c***purpose  evaluate the airy modulus and phase.
c***library   slatec (fnlib)
c***category  c10d
c***type      single precision (r9aimp-s, d9aimp-d)
c***keywords  airy function, fnlib, modulus, phase, special functions
c***author  fullerton, w., (lanl)
c***description
c
c evaluate the airy modulus and phase for x .le. -1.0
c
c series for am21       on the interval -1.25000d-01 to  0.
c                                        with weighted error   2.89e-17
c                                         log weighted error  16.54
c                               significant figures required  14.15
c                                    decimal places required  17.34
c
c series for ath1       on the interval -1.25000d-01 to  0.
c                                        with weighted error   2.53e-17
c                                         log weighted error  16.60
c                               significant figures required  15.15
c                                    decimal places required  17.38
c
c series for am22       on the interval -1.00000d+00 to -1.25000d-01
c                                        with weighted error   2.99e-17
c                                         log weighted error  16.52
c                               significant figures required  14.57
c                                    decimal places required  17.28
c
c series for ath2       on the interval -1.00000d+00 to -1.25000d-01
c                                        with weighted error   2.57e-17
c                                         log weighted error  16.59
c                               significant figures required  15.07
c                                    decimal places required  17.34
c
c***references  (none)
c***routines called  csevl, inits, r1mach, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890206  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900720  routine changed from user-callable to subsidiary.  (wrb)
c***end prologue  r9aimp
      dimension am21cs(40), ath1cs(36), am22cs(33), ath2cs(32)
      logical first
      save am21cs, ath1cs, am22cs, ath2cs, pi4, nam21,
     1 nath1, nam22, nath2, xsml, first
      data am21cs( 1) /    .0065809191 761485e0 /
      data am21cs( 2) /    .0023675984 685722e0 /
      data am21cs( 3) /    .0001324741 670371e0 /
      data am21cs( 4) /    .0000157600 904043e0 /
      data am21cs( 5) /    .0000027529 702663e0 /
      data am21cs( 6) /    .0000006102 679017e0 /
      data am21cs( 7) /    .0000001595 088468e0 /
      data am21cs( 8) /    .0000000471 033947e0 /
      data am21cs( 9) /    .0000000152 933871e0 /
      data am21cs(10) /    .0000000053 590722e0 /
      data am21cs(11) /    .0000000020 000910e0 /
      data am21cs(12) /    .0000000007 872292e0 /
      data am21cs(13) /    .0000000003 243103e0 /
      data am21cs(14) /    .0000000001 390106e0 /
      data am21cs(15) /    .0000000000 617011e0 /
      data am21cs(16) /    .0000000000 282491e0 /
      data am21cs(17) /    .0000000000 132979e0 /
      data am21cs(18) /    .0000000000 064188e0 /
      data am21cs(19) /    .0000000000 031697e0 /
      data am21cs(20) /    .0000000000 015981e0 /
      data am21cs(21) /    .0000000000 008213e0 /
      data am21cs(22) /    .0000000000 004296e0 /
      data am21cs(23) /    .0000000000 002284e0 /
      data am21cs(24) /    .0000000000 001232e0 /
      data am21cs(25) /    .0000000000 000675e0 /
      data am21cs(26) /    .0000000000 000374e0 /
      data am21cs(27) /    .0000000000 000210e0 /
      data am21cs(28) /    .0000000000 000119e0 /
      data am21cs(29) /    .0000000000 000068e0 /
      data am21cs(30) /    .0000000000 000039e0 /
      data am21cs(31) /    .0000000000 000023e0 /
      data am21cs(32) /    .0000000000 000013e0 /
      data am21cs(33) /    .0000000000 000008e0 /
      data am21cs(34) /    .0000000000 000005e0 /
      data am21cs(35) /    .0000000000 000003e0 /
      data am21cs(36) /    .0000000000 000001e0 /
      data am21cs(37) /    .0000000000 000001e0 /
      data am21cs(38) /    .0000000000 000000e0 /
      data am21cs(39) /    .0000000000 000000e0 /
      data am21cs(40) /    .0000000000 000000e0 /
      data ath1cs( 1) /   -.0712583781 5669365e0 /
      data ath1cs( 2) /   -.0059047197 9831451e0 /
      data ath1cs( 3) /   -.0001211454 4069499e0 /
      data ath1cs( 4) /   -.0000098860 8542270e0 /
      data ath1cs( 5) /   -.0000013808 4097352e0 /
      data ath1cs( 6) /   -.0000002614 2640172e0 /
      data ath1cs( 7) /   -.0000000605 0432589e0 /
      data ath1cs( 8) /   -.0000000161 8436223e0 /
      data ath1cs( 9) /   -.0000000048 3464911e0 /
      data ath1cs(10) /   -.0000000015 7655272e0 /
      data ath1cs(11) /   -.0000000005 5231518e0 /
      data ath1cs(12) /   -.0000000002 0545441e0 /
      data ath1cs(13) /   -.0000000000 8043412e0 /
      data ath1cs(14) /   -.0000000000 3291252e0 /
      data ath1cs(15) /   -.0000000000 1399875e0 /
      data ath1cs(16) /   -.0000000000 0616151e0 /
      data ath1cs(17) /   -.0000000000 0279614e0 /
      data ath1cs(18) /   -.0000000000 0130428e0 /
      data ath1cs(19) /   -.0000000000 0062373e0 /
      data ath1cs(20) /   -.0000000000 0030512e0 /
      data ath1cs(21) /   -.0000000000 0015239e0 /
      data ath1cs(22) /   -.0000000000 0007758e0 /
      data ath1cs(23) /   -.0000000000 0004020e0 /
      data ath1cs(24) /   -.0000000000 0002117e0 /
      data ath1cs(25) /   -.0000000000 0001132e0 /
      data ath1cs(26) /   -.0000000000 0000614e0 /
      data ath1cs(27) /   -.0000000000 0000337e0 /
      data ath1cs(28) /   -.0000000000 0000188e0 /
      data ath1cs(29) /   -.0000000000 0000105e0 /
      data ath1cs(30) /   -.0000000000 0000060e0 /
      data ath1cs(31) /   -.0000000000 0000034e0 /
      data ath1cs(32) /   -.0000000000 0000020e0 /
      data ath1cs(33) /   -.0000000000 0000011e0 /
      data ath1cs(34) /   -.0000000000 0000007e0 /
      data ath1cs(35) /   -.0000000000 0000004e0 /
      data ath1cs(36) /   -.0000000000 0000002e0 /
      data am22cs( 1) /   -.0156284448 0625341e0 /
      data am22cs( 2) /    .0077833644 5239681e0 /
      data am22cs( 3) /    .0008670577 7047718e0 /
      data am22cs( 4) /    .0001569662 7315611e0 /
      data am22cs( 5) /    .0000356396 2571432e0 /
      data am22cs( 6) /    .0000092459 8335425e0 /
      data am22cs( 7) /    .0000026211 0161850e0 /
      data am22cs( 8) /    .0000007918 8221651e0 /
      data am22cs( 9) /    .0000002510 4152792e0 /
      data am22cs(10) /    .0000000826 5223206e0 /
      data am22cs(11) /    .0000000280 5711662e0 /
      data am22cs(12) /    .0000000097 6821090e0 /
      data am22cs(13) /    .0000000034 7407923e0 /
      data am22cs(14) /    .0000000012 5828132e0 /
      data am22cs(15) /    .0000000004 6298826e0 /
      data am22cs(16) /    .0000000001 7272825e0 /
      data am22cs(17) /    .0000000000 6523192e0 /
      data am22cs(18) /    .0000000000 2490471e0 /
      data am22cs(19) /    .0000000000 0960156e0 /
      data am22cs(20) /    .0000000000 0373448e0 /
      data am22cs(21) /    .0000000000 0146417e0 /
      data am22cs(22) /    .0000000000 0057826e0 /
      data am22cs(23) /    .0000000000 0022991e0 /
      data am22cs(24) /    .0000000000 0009197e0 /
      data am22cs(25) /    .0000000000 0003700e0 /
      data am22cs(26) /    .0000000000 0001496e0 /
      data am22cs(27) /    .0000000000 0000608e0 /
      data am22cs(28) /    .0000000000 0000248e0 /
      data am22cs(29) /    .0000000000 0000101e0 /
      data am22cs(30) /    .0000000000 0000041e0 /
      data am22cs(31) /    .0000000000 0000017e0 /
      data am22cs(32) /    .0000000000 0000007e0 /
      data am22cs(33) /    .0000000000 0000002e0 /
      data ath2cs( 1) /    .0044052734 5871877e0 /
      data ath2cs( 2) /   -.0304291945 2318455e0 /
      data ath2cs( 3) /   -.0013856532 8377179e0 /
      data ath2cs( 4) /   -.0001804443 9089549e0 /
      data ath2cs( 5) /   -.0000338084 7108327e0 /
      data ath2cs( 6) /   -.0000076781 8353522e0 /
      data ath2cs( 7) /   -.0000019678 3944371e0 /
      data ath2cs( 8) /   -.0000005483 7271158e0 /
      data ath2cs( 9) /   -.0000001625 4615505e0 /
      data ath2cs(10) /   -.0000000505 3049981e0 /
      data ath2cs(11) /   -.0000000163 1580701e0 /
      data ath2cs(12) /   -.0000000054 3420411e0 /
      data ath2cs(13) /   -.0000000018 5739855e0 /
      data ath2cs(14) /   -.0000000006 4895120e0 /
      data ath2cs(15) /   -.0000000002 3105948e0 /
      data ath2cs(16) /   -.0000000000 8363282e0 /
      data ath2cs(17) /   -.0000000000 3071196e0 /
      data ath2cs(18) /   -.0000000000 1142367e0 /
      data ath2cs(19) /   -.0000000000 0429811e0 /
      data ath2cs(20) /   -.0000000000 0163389e0 /
      data ath2cs(21) /   -.0000000000 0062693e0 /
      data ath2cs(22) /   -.0000000000 0024260e0 /
      data ath2cs(23) /   -.0000000000 0009461e0 /
      data ath2cs(24) /   -.0000000000 0003716e0 /
      data ath2cs(25) /   -.0000000000 0001469e0 /
      data ath2cs(26) /   -.0000000000 0000584e0 /
      data ath2cs(27) /   -.0000000000 0000233e0 /
      data ath2cs(28) /   -.0000000000 0000093e0 /
      data ath2cs(29) /   -.0000000000 0000037e0 /
      data ath2cs(30) /   -.0000000000 0000015e0 /
      data ath2cs(31) /   -.0000000000 0000006e0 /
      data ath2cs(32) /   -.0000000000 0000002e0 /
      data pi4 / 0.7853981633 9744831 e0 /
      data first /.true./
c***first executable statement  r9aimp
      if (first) then
         eta = 0.1*r1mach(3)
         nam21 = inits (am21cs, 40, eta)
         nath1 = inits (ath1cs, 36, eta)
         nam22 = inits (am22cs, 33, eta)
         nath2 = inits (ath2cs, 32, eta)
c
         xsml = -1.0/r1mach(3)**0.3333
      endif
      first = .false.
c
      if (x.ge.(-2.0)) go to 20
      z = 1.0
      if (x.gt.xsml) z = 16.0/x**3 + 1.0
      ampl = 0.3125 + csevl(z, am21cs, nam21)
      theta = -0.625 + csevl (z, ath1cs, nath1)
      go to 30
c
 20   if (x .gt. (-1.0)) call xermsg ('slatec', 'r9aimp',
     +   'x must be le -1.0', 1, 2)
c
      z = (16.0/x**3 + 9.0)/7.0
      ampl = 0.3125 + csevl (z, am22cs, nam22)
      theta = -0.625 + csevl (z, ath2cs, nath2)
c
 30   sqrtx = sqrt(-x)
      ampl = sqrt (ampl/sqrtx)
      theta = pi4 - x*sqrtx * theta
c
      return
      end
