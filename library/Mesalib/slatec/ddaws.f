*deck ddaws
      double precision function ddaws (x)
c***begin prologue  ddaws
c***purpose  compute dawson's function.
c***library   slatec (fnlib)
c***category  c8c
c***type      double precision (daws-s, ddaws-d)
c***keywords  dawson's function, fnlib, special functions
c***author  fullerton, w., (lanl)
c***description
c
c ddaws(x) calculates the double precision dawson's integral
c for double precision argument x.
c
c series for daw        on the interval  0.          to  1.00000e+00
c                                        with weighted error   8.95e-32
c                                         log weighted error  31.05
c                               significant figures required  30.41
c                                    decimal places required  31.71
c
c series for daw2       on the interval  0.          to  1.60000e+01
c                                        with weighted error   1.61e-32
c                                         log weighted error  31.79
c                               significant figures required  31.40
c                                    decimal places required  32.62
c
c series for dawa       on the interval  0.          to  6.25000e-02
c                                        with weighted error   1.97e-32
c                                         log weighted error  31.71
c                               significant figures required  29.79
c                                    decimal places required  32.64
c
c***references  (none)
c***routines called  d1mach, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   780401  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   920618  removed space from variable names.  (rwc, wrb)
c***end prologue  ddaws
      double precision x, dawcs(21), daw2cs(45), dawacs(75), xbig,
     1  xmax, xsml, y, dcsevl, d1mach
      logical first
      save dawcs, daw2cs, dawacs, ntdaw, ntdaw2, ntdawa,
     1 xsml, xbig, xmax, first
      data dawcs(  1) / -.6351734375 1459492010 6512773629 3 d-2      /
      data dawcs(  2) / -.2294071479 6773869398 9982412586 6 d+0      /
      data dawcs(  3) / +.2213050093 9084764416 8397916178 6 d-1      /
      data dawcs(  4) / -.1549265453 8929850467 4305775337 5 d-2      /
      data dawcs(  5) / +.8497327715 6849174567 7754294806 6 d-4      /
      data dawcs(  6) / -.3828266270 9720149249 9409952130 9 d-5      /
      data dawcs(  7) / +.1462854806 2501631977 5714894953 9 d-6      /
      data dawcs(  8) / -.4851982381 8259917988 4671542511 4 d-8      /
      data dawcs(  9) / +.1421463577 7591397903 4756818330 4 d-9      /
      data dawcs( 10) / -.3728836087 9205965253 3549305408 8 d-11     /
      data dawcs( 11) / +.8854942961 7782033701 9456523136 9 d-13     /
      data dawcs( 12) / -.1920757131 3502063554 2164841749 3 d-14     /
      data dawcs( 13) / +.3834325867 2463275882 4107443925 3 d-16     /
      data dawcs( 14) / -.7089154168 1758816335 8409932799 9 d-18     /
      data dawcs( 15) / +.1220552135 8894576744 1690112000 0 d-19     /
      data dawcs( 16) / -.1966204826 6053487602 9945173333 3 d-21     /
      data dawcs( 17) / +.2975845541 3765971891 1317333333 3 d-23     /
      data dawcs( 18) / -.4247069514 8005969510 3999999999 9 d-25     /
      data dawcs( 19) / +.5734270767 3917427985 0666666666 6 d-27     /
      data dawcs( 20) / -.7345836823 1784502613 3333333333 3 d-29     /
      data dawcs( 21) / +.8951937667 5165525333 3333333333 3 d-31     /
      data daw2cs(  1) / -.5688654410 5215527114 1605337336 74 d-1     /
      data daw2cs(  2) / -.3181134699 6168131279 3228780488 22 d+0     /
      data daw2cs(  3) / +.2087384541 3642236789 7415801988 58 d+0     /
      data daw2cs(  4) / -.1247540991 3779131214 0734983147 84 d+0     /
      data daw2cs(  5) / +.6786930518 6676777092 8475164236 76 d-1     /
      data daw2cs(  6) / -.3365914489 5270939503 0682309665 87 d-1     /
      data daw2cs(  7) / +.1526078127 1987971743 6824603816 40 d-1     /
      data daw2cs(  8) / -.6348370962 5962148230 5860947885 35 d-2     /
      data daw2cs(  9) / +.2432674092 0748520596 8659661093 43 d-2     /
      data daw2cs( 10) / -.8621954149 1065032038 5269835496 37 d-3     /
      data daw2cs( 11) / +.2837657333 6321625302 8576365382 95 d-3     /
      data daw2cs( 12) / -.8705754987 4170423699 3965814643 35 d-4     /
      data daw2cs( 13) / +.2498684998 5481658331 8000441372 76 d-4     /
      data daw2cs( 14) / -.6731928676 4160294344 6030503395 20 d-5     /
      data daw2cs( 15) / +.1707857878 5573543710 5045240478 44 d-5     /
      data daw2cs( 16) / -.4091755122 6475381271 8965924900 38 d-6     /
      data daw2cs( 17) / +.9282829221 6755773260 7517853122 73 d-7     /
      data daw2cs( 18) / -.1999140361 0147617829 8450963321 98 d-7     /
      data daw2cs( 19) / +.4096349064 4082195241 2104878689 17 d-8     /
      data daw2cs( 20) / -.8003240954 0993168075 7067817535 61 d-9     /
      data daw2cs( 21) / +.1493850312 8761465059 1432255501 10 d-9     /
      data daw2cs( 22) / -.2668799988 5622329284 9246510633 39 d-10    /
      data daw2cs( 23) / +.4571221698 5159458151 4056177241 03 d-11    /
      data daw2cs( 24) / -.7518730522 2043565872 2437273267 71 d-12    /
      data daw2cs( 25) / +.1189310005 2629681879 0298289873 02 d-12    /
      data daw2cs( 26) / -.1811690793 3852346973 4903182630 84 d-13    /
      data daw2cs( 27) / +.2661173368 4358969193 0016121996 26 d-14    /
      data daw2cs( 28) / -.3773886305 2129419795 4441099059 30 d-15    /
      data daw2cs( 29) / +.5172795378 9087172679 6800822293 29 d-16    /
      data daw2cs( 30) / -.6860368408 4077500979 4195646701 02 d-17    /
      data daw2cs( 31) / +.8812375135 4161071806 4693373217 45 d-18    /
      data daw2cs( 32) / -.1097424824 9996606292 1062996246 52 d-18    /
      data daw2cs( 33) / +.1326119932 6367178513 5955458916 35 d-19    /
      data daw2cs( 34) / -.1556273276 8137380785 4887765715 62 d-20    /
      data daw2cs( 35) / +.1775142558 3655720607 8334155707 73 d-21    /
      data daw2cs( 36) / -.1969500696 7006578384 9536087654 39 d-22    /
      data daw2cs( 37) / +.2127007489 6998699661 9240101205 33 d-23    /
      data daw2cs( 38) / -.2237539812 4627973794 1821139626 66 d-24    /
      data daw2cs( 39) / +.2294276857 8582348946 9713831253 33 d-25    /
      data daw2cs( 40) / -.2294378884 6552928693 3295923199 99 d-26    /
      data daw2cs( 41) / +.2239170210 0592453618 3422976000 00 d-27    /
      data daw2cs( 42) / -.2133823061 6608897703 6782250666 66 d-28    /
      data daw2cs( 43) / +.1986619658 5123531518 0284586666 66 d-29    /
      data daw2cs( 44) / -.1807929586 6694391771 9551999999 99 d-30    /
      data daw2cs( 45) / +.1609068601 5283030305 4506666666 66 d-31    /
      data dawacs(  1) / +.1690485637 7657037554 2263743884 9 d-1      /
      data dawacs(  2) / +.8683252278 4069579905 3610785076 8 d-2      /
      data dawacs(  3) / +.2424864042 4177154532 7770345988 9 d-3      /
      data dawacs(  4) / +.1261182399 5726900016 5194924037 7 d-4      /
      data dawacs(  5) / +.1066453314 6361769557 0569112590 6 d-5      /
      data dawacs(  6) / +.1358159794 7907276113 4842450572 8 d-6      /
      data dawacs(  7) / +.2171042356 5772983989 0431274474 3 d-7      /
      data dawacs(  8) / +.2867010501 8052952703 4367680481 3 d-8      /
      data dawacs(  9) / -.1901336393 0358201122 8249237802 4 d-9      /
      data dawacs( 10) / -.3097780484 3952011255 3206577426 8 d-9      /
      data dawacs( 11) / -.1029414876 0575092473 9813228641 3 d-9      /
      data dawacs( 12) / -.6260356459 4595761504 1758728312 1 d-11     /
      data dawacs( 13) / +.8563132497 4464512162 6230316627 6 d-11     /
      data dawacs( 14) / +.3033045148 0756592929 7626627625 7 d-11     /
      data dawacs( 15) / -.2523618306 8092913726 3088693882 6 d-12     /
      data dawacs( 16) / -.4210604795 4406645131 7546193451 0 d-12     /
      data dawacs( 17) / -.4431140826 6462383121 4342945203 6 d-13     /
      data dawacs( 18) / +.4911210272 8412052059 4003706511 7 d-13     /
      data dawacs( 19) / +.1235856242 2839034070 7647795473 9 d-13     /
      data dawacs( 20) / -.5788733199 0165692469 5576507106 9 d-14     /
      data dawacs( 21) / -.2282723294 8073586209 7818395703 0 d-14     /
      data dawacs( 22) / +.7637149411 0141264763 1236291759 0 d-15     /
      data dawacs( 23) / +.3851546883 5668117287 7759400209 5 d-15     /
      data dawacs( 24) / -.1199932056 9282905928 0323728304 5 d-15     /
      data dawacs( 25) / -.6313439150 0945723473 3427028525 0 d-16     /
      data dawacs( 26) / +.2239559965 9729753752 5491279023 7 d-16     /
      data dawacs( 27) / +.9987925830 0764959951 3289120074 9 d-17     /
      data dawacs( 28) / -.4681068274 3224953345 3624650725 2 d-17     /
      data dawacs( 29) / -.1436303644 3497213372 4162875153 4 d-17     /
      data dawacs( 30) / +.1020822731 4105411129 7790803213 0 d-17     /
      data dawacs( 31) / +.1538908873 1360920728 3738982237 2 d-18     /
      data dawacs( 32) / -.2189157877 6457938888 9479092605 6 d-18     /
      data dawacs( 33) / +.2156879197 9386517503 9235915251 7 d-20     /
      data dawacs( 34) / +.4370219827 4424498511 3479255739 5 d-19     /
      data dawacs( 35) / -.8234581460 9772072410 9892790517 7 d-20     /
      data dawacs( 36) / -.7498648721 2564662229 0320283542 0 d-20     /
      data dawacs( 37) / +.3282536720 7356716109 5761293003 9 d-20     /
      data dawacs( 38) / +.8858064309 5039211160 7656151515 1 d-21     /
      data dawacs( 39) / -.9185087111 7270029880 9446053148 5 d-21     /
      data dawacs( 40) / +.2978962223 7887489883 1416604579 1 d-22     /
      data dawacs( 41) / +.1972132136 6184718831 5950546804 1 d-21     /
      data dawacs( 42) / -.5974775596 3629066380 8958499511 7 d-22     /
      data dawacs( 43) / -.2834410031 5038509654 4382518244 1 d-22     /
      data dawacs( 44) / +.2209560791 1315545147 7715048901 2 d-22     /
      data dawacs( 45) / -.5439955741 8971443000 7948030771 1 d-25     /
      data dawacs( 46) / -.5213549243 2948486680 1713669647 0 d-23     /
      data dawacs( 47) / +.1702350556 8131141990 6567149907 6 d-23     /
      data dawacs( 48) / +.6917400860 8361483430 2218566019 7 d-24     /
      data dawacs( 49) / -.6540941793 0027525122 3944512580 2 d-24     /
      data dawacs( 50) / +.6093576580 4393289603 7182465463 6 d-25     /
      data dawacs( 51) / +.1408070432 9051874615 0194508027 2 d-24     /
      data dawacs( 52) / -.6785886121 0548463311 6767494375 5 d-25     /
      data dawacs( 53) / -.9799732036 2142957117 4158310222 5 d-26     /
      data dawacs( 54) / +.2121244903 0990413325 9896093916 0 d-25     /
      data dawacs( 55) / -.5954455022 5487909382 3880215448 7 d-26     /
      data dawacs( 56) / -.3093088861 8754701778 3884723204 9 d-26     /
      data dawacs( 57) / +.2854389216 3445246824 0069198610 4 d-26     /
      data dawacs( 58) / -.3951289447 3793055660 2347727181 1 d-27     /
      data dawacs( 59) / -.5906000648 6076284781 1684089445 3 d-27     /
      data dawacs( 60) / +.3670236964 6686870036 4788998060 9 d-27     /
      data dawacs( 61) / -.4839958238 0422762565 9830303894 1 d-29     /
      data dawacs( 62) / -.9799265984 2104438695 9740401702 2 d-28     /
      data dawacs( 63) / +.4684773732 6121306061 5890880430 0 d-28     /
      data dawacs( 64) / +.5030877696 9934610516 4766760315 5 d-29     /
      data dawacs( 65) / -.1547395051 7060282392 4755206829 5 d-28     /
      data dawacs( 66) / +.6112180185 0864192439 7600566271 4 d-29     /
      data dawacs( 67) / +.1357913399 1248116503 4360273615 8 d-29     /
      data dawacs( 68) / -.2417687752 7686730883 8530429904 4 d-29     /
      data dawacs( 69) / +.8369074582 0742989452 9288758729 1 d-30     /
      data dawacs( 70) / +.2665413042 7889791658 3831940156 6 d-30     /
      data dawacs( 71) / -.3811653692 3548903369 3569100371 2 d-30     /
      data dawacs( 72) / +.1230054721 8849514643 7170687258 5 d-30     /
      data dawacs( 73) / +.4622506399 0414935088 0553692998 3 d-31     /
      data dawacs( 74) / -.6120087296 8816777229 1143559300 1 d-31     /
      data dawacs( 75) / +.1966024640 1931646869 5623021789 6 d-31     /
      data first /.true./
c***first executable statement  ddaws
      if (first) then
         eps = d1mach(3)
         ntdaw = initds (dawcs, 21, 0.1*eps)
         ntdaw2 = initds (daw2cs, 45, 0.1*eps)
         ntdawa = initds (dawacs, 75, 0.1*eps)
c
         xsml = sqrt(1.5*eps)
         xbig = sqrt (0.5/eps)
         xmax = exp (min (-log(2.d0*d1mach(1)), log(d1mach(2)))
     1     - 0.001d0)
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.1.0d0) go to 20
c
      ddaws = x
      if (y.le.xsml) return
c
      ddaws = x * (.75d0 + dcsevl (2.d0*y*y-1.d0, dawcs, ntdaw))
      return
c
 20   if (y.gt.4.d0) go to 30
      ddaws = x * (.25d0 + dcsevl (.125d0*y*y-1.d0, daw2cs, ntdaw2))
      return
c
 30   if (y.gt.xmax) go to 40
      ddaws = 0.5d0/x
      if (y.gt.xbig) return
c
      ddaws = (0.5d0 + dcsevl (32.d0/y**2-1.d0, dawacs, ntdawa)) / x
      return
c
 40   call xermsg ('slatec', 'ddaws', 'abs(x) so large daws underflows',
     +   1, 1)
      ddaws = 0.0d0
      return
c
      end
