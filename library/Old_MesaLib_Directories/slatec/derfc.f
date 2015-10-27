*deck derfc
      double precision function derfc (x)
c***begin prologue  derfc
c***purpose  compute the complementary error function.
c***library   slatec (fnlib)
c***category  c8a, l5a1e
c***type      double precision (erfc-s, derfc-d)
c***keywords  complementary error function, erfc, fnlib,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c derfc(x) calculates the double precision complementary error function
c for double precision argument x.
c
c series for erf        on the interval  0.          to  1.00000e+00
c                                        with weighted error   1.28e-32
c                                         log weighted error  31.89
c                               significant figures required  31.05
c                                    decimal places required  32.55
c
c series for erc2       on the interval  2.50000e-01 to  1.00000e+00
c                                        with weighted error   2.67e-32
c                                         log weighted error  31.57
c                               significant figures required  30.31
c                                    decimal places required  32.42
c
c series for erfc       on the interval  0.          to  2.50000e-01
c                                        with weighted error   1.53e-31
c                                         log weighted error  30.82
c                               significant figures required  29.47
c                                    decimal places required  31.70
c
c***references  (none)
c***routines called  d1mach, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920618  removed space from variable names.  (rwc, wrb)
c***end prologue  derfc
      double precision x, erfcs(21), erfccs(59), erc2cs(49), sqeps,
     1  sqrtpi, xmax, txmax, xsml, y, d1mach, dcsevl
      logical first
      save erfcs, erc2cs, erfccs, sqrtpi, nterf,
     1 nterfc, nterc2, xsml, xmax, sqeps, first
      data erfcs(  1) / -.4904612123 4691808039 9845440333 76 d-1     /
      data erfcs(  2) / -.1422612051 0371364237 8247418996 31 d+0     /
      data erfcs(  3) / +.1003558218 7599795575 7546767129 33 d-1     /
      data erfcs(  4) / -.5768764699 7674847650 8270255091 67 d-3     /
      data erfcs(  5) / +.2741993125 2196061034 4221607914 71 d-4     /
      data erfcs(  6) / -.1104317550 7344507604 1353812959 05 d-5     /
      data erfcs(  7) / +.3848875542 0345036949 9613114981 74 d-7     /
      data erfcs(  8) / -.1180858253 3875466969 6317518015 81 d-8     /
      data erfcs(  9) / +.3233421582 6050909646 4029309533 54 d-10    /
      data erfcs( 10) / -.7991015947 0045487581 6073747085 95 d-12    /
      data erfcs( 11) / +.1799072511 3961455611 9672454866 34 d-13    /
      data erfcs( 12) / -.3718635487 8186926382 3168282094 93 d-15    /
      data erfcs( 13) / +.7103599003 7142529711 6899083946 66 d-17    /
      data erfcs( 14) / -.1261245511 9155225832 4954248533 33 d-18    /
      data erfcs( 15) / +.2091640694 1769294369 1705002666 66 d-20    /
      data erfcs( 16) / -.3253973102 9314072982 3641600000 00 d-22    /
      data erfcs( 17) / +.4766867209 7976748332 3733333333 33 d-24    /
      data erfcs( 18) / -.6598012078 2851343155 1999999999 99 d-26    /
      data erfcs( 19) / +.8655011469 9637626197 3333333333 33 d-28    /
      data erfcs( 20) / -.1078892517 7498064213 3333333333 33 d-29    /
      data erfcs( 21) / +.1281188399 3017002666 6666666666 66 d-31    /
      data erc2cs(  1) / -.6960134660 2309501127 3915082619 7 d-1      /
      data erc2cs(  2) / -.4110133936 2620893489 8221208466 6 d-1      /
      data erc2cs(  3) / +.3914495866 6896268815 6114370524 4 d-2      /
      data erc2cs(  4) / -.4906395650 5489791612 8093545077 4 d-3      /
      data erc2cs(  5) / +.7157479001 3770363807 6089414182 5 d-4      /
      data erc2cs(  6) / -.1153071634 1312328338 0823284791 2 d-4      /
      data erc2cs(  7) / +.1994670590 2019976350 5231486770 9 d-5      /
      data erc2cs(  8) / -.3642666471 5992228739 3611843071 1 d-6      /
      data erc2cs(  9) / +.6944372610 0050125899 3127721463 3 d-7      /
      data erc2cs( 10) / -.1371220902 1043660195 3460514121 0 d-7      /
      data erc2cs( 11) / +.2788389661 0071371319 6386034808 7 d-8      /
      data erc2cs( 12) / -.5814164724 3311615518 6479105031 6 d-9      /
      data erc2cs( 13) / +.1238920491 7527531811 8016881795 0 d-9      /
      data erc2cs( 14) / -.2690639145 3067434323 9042493788 9 d-10     /
      data erc2cs( 15) / +.5942614350 8479109824 4470968384 0 d-11     /
      data erc2cs( 16) / -.1332386735 7581195792 8775442057 0 d-11     /
      data erc2cs( 17) / +.3028046806 1771320171 7369724330 4 d-12     /
      data erc2cs( 18) / -.6966648814 9410325887 9586758895 4 d-13     /
      data erc2cs( 19) / +.1620854541 0539229698 1289322762 8 d-13     /
      data erc2cs( 20) / -.3809934465 2504919998 7691305772 9 d-14     /
      data erc2cs( 21) / +.9040487815 9788311493 6897101297 5 d-15     /
      data erc2cs( 22) / -.2164006195 0896073478 0981204700 3 d-15     /
      data erc2cs( 23) / +.5222102233 9958549846 0798024417 2 d-16     /
      data erc2cs( 24) / -.1269729602 3645553363 7241552778 0 d-16     /
      data erc2cs( 25) / +.3109145504 2761975838 3622741295 1 d-17     /
      data erc2cs( 26) / -.7663762920 3203855240 0956671481 1 d-18     /
      data erc2cs( 27) / +.1900819251 3627452025 3692973329 0 d-18     /
      data erc2cs( 28) / -.4742207279 0690395452 2565599996 5 d-19     /
      data erc2cs( 29) / +.1189649200 0765283828 8068307845 1 d-19     /
      data erc2cs( 30) / -.3000035590 3257802568 4527131306 6 d-20     /
      data erc2cs( 31) / +.7602993453 0432461730 1938527709 8 d-21     /
      data erc2cs( 32) / -.1935909447 6068728815 6981104913 0 d-21     /
      data erc2cs( 33) / +.4951399124 7733378810 0004238677 3 d-22     /
      data erc2cs( 34) / -.1271807481 3363718796 0862198988 8 d-22     /
      data erc2cs( 35) / +.3280049600 4695130433 1584165205 3 d-23     /
      data erc2cs( 36) / -.8492320176 8228965689 2479242239 9 d-24     /
      data erc2cs( 37) / +.2206917892 8075602235 1987998719 9 d-24     /
      data erc2cs( 38) / -.5755617245 6965284983 1281950719 9 d-25     /
      data erc2cs( 39) / +.1506191533 6392342503 5414405119 9 d-25     /
      data erc2cs( 40) / -.3954502959 0187969531 0428569599 9 d-26     /
      data erc2cs( 41) / +.1041529704 1515009799 8464505173 3 d-26     /
      data erc2cs( 42) / -.2751487795 2787650794 5017890133 3 d-27     /
      data erc2cs( 43) / +.7290058205 4975574089 9770368000 0 d-28     /
      data erc2cs( 44) / -.1936939645 9159478040 7750109866 6 d-28     /
      data erc2cs( 45) / +.5160357112 0514872983 7005482666 6 d-29     /
      data erc2cs( 46) / -.1378419322 1930940993 8964480000 0 d-29     /
      data erc2cs( 47) / +.3691326793 1070690422 5109333333 3 d-30     /
      data erc2cs( 48) / -.9909389590 6243654206 5322666666 6 d-31     /
      data erc2cs( 49) / +.2666491705 1953884133 2394666666 6 d-31     /
      data erfccs(  1) / +.7151793102 0292477450 3697709496 d-1        /
      data erfccs(  2) / -.2653243433 7606715755 8893386681 d-1        /
      data erfccs(  3) / +.1711153977 9208558833 2699194606 d-2        /
      data erfccs(  4) / -.1637516634 5851788416 3746404749 d-3        /
      data erfccs(  5) / +.1987129350 0552036499 5974806758 d-4        /
      data erfccs(  6) / -.2843712412 7665550875 0175183152 d-5        /
      data erfccs(  7) / +.4606161308 9631303696 9379968464 d-6        /
      data erfccs(  8) / -.8227753025 8792084205 7766536366 d-7        /
      data erfccs(  9) / +.1592141872 7709011298 9358340826 d-7        /
      data erfccs( 10) / -.3295071362 2528432148 6631665072 d-8        /
      data erfccs( 11) / +.7223439760 4005554658 1261153890 d-9        /
      data erfccs( 12) / -.1664855813 3987295934 4695966886 d-9        /
      data erfccs( 13) / +.4010392588 2376648207 7671768814 d-10       /
      data erfccs( 14) / -.1004816214 4257311327 2170176283 d-10       /
      data erfccs( 15) / +.2608275913 3003338085 9341009439 d-11       /
      data erfccs( 16) / -.6991110560 4040248655 7697812476 d-12       /
      data erfccs( 17) / +.1929492333 2617070862 4205749803 d-12       /
      data erfccs( 18) / -.5470131188 7543310649 0125085271 d-13       /
      data erfccs( 19) / +.1589663309 7626974483 9084032762 d-13       /
      data erfccs( 20) / -.4726893980 1975548392 0369584290 d-14       /
      data erfccs( 21) / +.1435873376 7849847867 2873997840 d-14       /
      data erfccs( 22) / -.4449510561 8173583941 7250062829 d-15       /
      data erfccs( 23) / +.1404810884 7682334373 7305537466 d-15       /
      data erfccs( 24) / -.4513818387 7642108962 5963281623 d-16       /
      data erfccs( 25) / +.1474521541 0451330778 7018713262 d-16       /
      data erfccs( 26) / -.4892621406 9457761543 6841552532 d-17       /
      data erfccs( 27) / +.1647612141 4106467389 5301522827 d-17       /
      data erfccs( 28) / -.5626817176 3294080929 9928521323 d-18       /
      data erfccs( 29) / +.1947443382 2320785142 9197867821 d-18       /
      data erfccs( 30) / -.6826305642 9484207295 6664144723 d-19       /
      data erfccs( 31) / +.2421988887 2986492401 8301125438 d-19       /
      data erfccs( 32) / -.8693414133 5030704256 3800861857 d-20       /
      data erfccs( 33) / +.3155180346 2280855712 2363401262 d-20       /
      data erfccs( 34) / -.1157372324 0496087426 1239486742 d-20       /
      data erfccs( 35) / +.4288947161 6056539462 3737097442 d-21       /
      data erfccs( 36) / -.1605030742 0576168500 5737770964 d-21       /
      data erfccs( 37) / +.6063298757 4538026449 5069923027 d-22       /
      data erfccs( 38) / -.2311404251 6979584909 8840801367 d-22       /
      data erfccs( 39) / +.8888778540 6618855255 4702955697 d-23       /
      data erfccs( 40) / -.3447260576 6513765223 0718495566 d-23       /
      data erfccs( 41) / +.1347865460 2069650682 7582774181 d-23       /
      data erfccs( 42) / -.5311794071 1250217364 5873201807 d-24       /
      data erfccs( 43) / +.2109341058 6197831682 8954734537 d-24       /
      data erfccs( 44) / -.8438365587 9237891159 8133256738 d-25       /
      data erfccs( 45) / +.3399982524 9452089062 7359576337 d-25       /
      data erfccs( 46) / -.1379452388 0732420900 2238377110 d-25       /
      data erfccs( 47) / +.5634490311 8332526151 3392634811 d-26       /
      data erfccs( 48) / -.2316490434 4770654482 3427752700 d-26       /
      data erfccs( 49) / +.9584462844 6018101526 3158381226 d-27       /
      data erfccs( 50) / -.3990722880 3301097262 4224850193 d-27       /
      data erfccs( 51) / +.1672129225 9444773601 7228709669 d-27       /
      data erfccs( 52) / -.7045991522 7660138563 8803782587 d-28       /
      data erfccs( 53) / +.2979768402 8642063541 2357989444 d-28       /
      data erfccs( 54) / -.1262522466 4606192972 2422632994 d-28       /
      data erfccs( 55) / +.5395438704 5424879398 5299653154 d-29       /
      data erfccs( 56) / -.2380992882 5314591867 5346190062 d-29       /
      data erfccs( 57) / +.1099052830 1027615735 9726683750 d-29       /
      data erfccs( 58) / -.4867713741 6449657273 2518677435 d-30       /
      data erfccs( 59) / +.1525877264 1103575676 3200828211 d-30       /
      data sqrtpi / 1.772453850 9055160272 9816748334 115d0 /
      data first /.true./
c***first executable statement  derfc
      if (first) then
         eta = 0.1*real(d1mach(3))
         nterf = initds (erfcs, 21, eta)
         nterfc = initds (erfccs, 59, eta)
         nterc2 = initds (erc2cs, 49, eta)
c
         xsml = -sqrt(-log(sqrtpi*d1mach(3)))
         txmax = sqrt(-log(sqrtpi*d1mach(1)))
         xmax = txmax - 0.5d0*log(txmax)/txmax - 0.01d0
         sqeps = sqrt(2.0d0*d1mach(3))
      endif
      first = .false.
c
      if (x.gt.xsml) go to 20
c
c erfc(x) = 1.0 - erf(x)  for  x .lt. xsml
c
      derfc = 2.0d0
      return
c
 20   if (x.gt.xmax) go to 40
      y = abs(x)
      if (y.gt.1.0d0) go to 30
c
c erfc(x) = 1.0 - erf(x)  for abs(x) .le. 1.0
c
      if (y.lt.sqeps) derfc = 1.0d0 - 2.0d0*x/sqrtpi
      if (y.ge.sqeps) derfc = 1.0d0 - x*(1.0d0 + dcsevl (2.d0*x*x-1.d0,
     1  erfcs, nterf))
      return
c
c erfc(x) = 1.0 - erf(x)  for  1.0 .lt. abs(x) .le. xmax
c
 30   y = y*y
      if (y.le.4.d0) derfc = exp(-y)/abs(x) * (0.5d0 + dcsevl (
     1  (8.d0/y-5.d0)/3.d0, erc2cs, nterc2) )
      if (y.gt.4.d0) derfc = exp(-y)/abs(x) * (0.5d0 + dcsevl (
     1  8.d0/y-1.d0, erfccs, nterfc) )
      if (x.lt.0.d0) derfc = 2.0d0 - derfc
      return
c
 40   call xermsg ('slatec', 'derfc', 'x so big erfc underflows', 1, 1)
      derfc = 0.d0
      return
c
      end
