*deck dbie
      double precision function dbie (x)
c***begin prologue  dbie
c***purpose  calculate the bairy function for a negative argument and an
c            exponentially scaled bairy function for a non-negative
c            argument.
c***library   slatec (fnlib)
c***category  c10d
c***type      double precision (bie-s, dbie-d)
c***keywords  bairy function, exponentially scaled, fnlib,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c dbie(x) calculates the double precision airy function of the
c second kind or the double precision exponentially scaled airy
c function of the second kind, depending on the value of the
c double precision argument x.
c
c evaluate bi(x) for x .le. 0.0  and  bi(x)*exp(-zeta)  where
c zeta = 2/3 * x**(3/2)  for x .ge. 0.0
c
c
c series for bif        on the interval -1.00000e+00 to  1.00000e+00
c                                        with weighted error   1.45e-32
c                                         log weighted error  31.84
c                               significant figures required  30.85
c                                    decimal places required  32.40
c
c
c series for big        on the interval -1.00000e+00 to  1.00000e+00
c                                        with weighted error   1.29e-33
c                                         log weighted error  32.89
c                               significant figures required  31.48
c                                    decimal places required  33.45
c
c
c series for bif2       on the interval  1.00000e+00 to  8.00000e+00
c                                        with weighted error   6.08e-32
c                                         log weighted error  31.22
c                        approx significant figures required  30.8
c                                    decimal places required  31.80
c
c
c series for big2       on the interval  1.00000e+00 to  8.00000e+00
c                                        with weighted error   4.91e-33
c                                         log weighted error  32.31
c                        approx significant figures required  31.6
c                                    decimal places required  32.90
c
c
c series for bip1       on the interval  1.25000e-01 to  3.53553e-01
c                                        with weighted error   1.06e-32
c                                         log weighted error  31.98
c                               significant figures required  30.61
c                                    decimal places required  32.81
c
c
c series for bip2       on the interval  0.          to  1.25000e-01
c                                        with weighted error   4.04e-33
c                                         log weighted error  32.39
c                               significant figures required  31.15
c                                    decimal places required  33.37
c
c***references  (none)
c***routines called  d1mach, d9aimp, dcsevl, initds
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  dbie
      double precision x, bifcs(13), bigcs(13), bif2cs(15), big2cs(15),
     1  bip1cs(47), bip2cs(88), atr, btr, sqrtx, theta, xbig, xm, x3sml,
     2  x32sml, z, d1mach, dcsevl
      logical first
      save bifcs, bigcs, bif2cs, big2cs, bip1cs, bip2cs, atr, btr,
     1 nbif, nbig, nbif2, nbig2, nbip1, nbip2, x3sml, x32sml, xbig,
     2 first
      data bifcs(  1) / -.1673021647 1986649483 5374239281 76 d-1     /
      data bifcs(  2) / +.1025233583 4249445611 4263627777 57 d+0     /
      data bifcs(  3) / +.1708309250 7381516539 4296502420 13 d-2     /
      data bifcs(  4) / +.1186254546 7744681179 2164592100 40 d-4     /
      data bifcs(  5) / +.4493290701 7792133694 5318879272 42 d-7     /
      data bifcs(  6) / +.1069820714 3387889067 5677676636 28 d-9     /
      data bifcs(  7) / +.1748064339 9771824706 0105176285 73 d-12    /
      data bifcs(  8) / +.2081023107 1761711025 8818918343 99 d-15    /
      data bifcs(  9) / +.1884981469 5665416509 9279717333 33 d-18    /
      data bifcs( 10) / +.1342577917 3097804625 8826666666 66 d-21    /
      data bifcs( 11) / +.7715959342 9658887893 3333333333 33 d-25    /
      data bifcs( 12) / +.3653387961 7478566399 9999999999 99 d-28    /
      data bifcs( 13) / +.1449756592 7953066666 6666666666 66 d-31    /
      data bigcs(  1) / +.2246622324 8574522283 4682201390 24 d-1     /
      data bigcs(  2) / +.3736477545 3019545441 7275616667 52 d-1     /
      data bigcs(  3) / +.4447621895 7212285696 2152943266 39 d-3     /
      data bigcs(  4) / +.2470807563 6329384245 4945919488 82 d-5     /
      data bigcs(  5) / +.7919135339 5149635134 8624262855 96 d-8     /
      data bigcs(  6) / +.1649807985 1827779880 8878724027 06 d-10    /
      data bigcs(  7) / +.2411990666 4835455909 2475011228 41 d-13    /
      data bigcs(  8) / +.2610373623 6091436985 1847812693 33 d-16    /
      data bigcs(  9) / +.2175308297 7160323853 1237920000 00 d-19    /
      data bigcs( 10) / +.1438694640 0390433219 4837333333 33 d-22    /
      data bigcs( 11) / +.7734912561 2083468629 3333333333 33 d-26    /
      data bigcs( 12) / +.3446929203 3849002666 6666666666 66 d-29    /
      data bigcs( 13) / +.1293891927 3216000000 0000000000 00 d-32    /
      data bif2cs(  1) / +.0998457269 3816041044 6828425799 3 d+0      /
      data bif2cs(  2) / +.4786249778 6300553772 2114673182 31 d+0     /
      data bif2cs(  3) / +.2515521196 0433011771 3244154366 75 d-1     /
      data bif2cs(  4) / +.5820693885 2326456396 5156978722 16 d-3     /
      data bif2cs(  5) / +.7499765964 4377865943 8614573782 17 d-5     /
      data bif2cs(  6) / +.6134602870 3493836681 4030103564 74 d-7     /
      data bif2cs(  7) / +.3462753885 1480632900 4342687333 59 d-9     /
      data bif2cs(  8) / +.1428891008 0270254287 7708467489 31 d-11    /
      data bif2cs(  9) / +.4496270429 8334641895 0564721792 00 d-14    /
      data bif2cs( 10) / +.1114232306 5833011708 4283001066 66 d-16    /
      data bif2cs( 11) / +.2230479106 6175002081 5178666666 66 d-19    /
      data bif2cs( 12) / +.3681577873 6393142842 9226666666 66 d-22    /
      data bif2cs( 13) / +.5096086844 9338261333 3333333333 33 d-25    /
      data bif2cs( 14) / +.6000338692 6288554666 6666666666 66 d-28    /
      data bif2cs( 15) / +.6082749744 6570666666 6666666666 66 d-31    /
      data big2cs(  1) / +.0333056621 4551434046 5176188111 647 d+0    /
      data big2cs(  2) / +.1613092151 2319706761 3287532084 943 d+0    /
      data big2cs(  3) / +.6319007309 6134286912 1615634921 173 d-2    /
      data big2cs(  4) / +.1187904568 1625173638 9780192304 567 d-3    /
      data big2cs(  5) / +.1304534588 6200265614 7116485012 843 d-5    /
      data big2cs(  6) / +.9374125995 5352172954 6809615508 936 d-8    /
      data big2cs(  7) / +.4745801886 7472515378 8510169834 595 d-10   /
      data big2cs(  8) / +.1783107265 0948139980 0065667560 946 d-12   /
      data big2cs(  9) / +.5167591927 8495818037 4276356640 000 d-15   /
      data big2cs( 10) / +.1190045083 8682712512 9496251733 333 d-17   /
      data big2cs( 11) / +.2229828806 6640351727 7063466666 666 d-20   /
      data big2cs( 12) / +.3465519230 2768941972 2666666666 666 d-23   /
      data big2cs( 13) / +.4539263363 2050451413 3333333333 333 d-26   /
      data big2cs( 14) / +.5078849965 1352234666 6666666666 666 d-29   /
      data big2cs( 15) / +.4910206746 9653333333 3333333333 333 d-32   /
      data bip1cs(  1) / -.8322047477 9434474687 4718647079 73 d-1     /
      data bip1cs(  2) / +.1146118927 3711742889 9202261280 31 d-1     /
      data bip1cs(  3) / +.4289644071 8911509494 1344725666 35 d-3     /
      data bip1cs(  4) / -.1490663937 9950514017 8476777329 54 d-3     /
      data bip1cs(  5) / -.1307659726 7876290663 1363409988 81 d-4     /
      data bip1cs(  6) / +.6327598396 1030344754 5357160324 94 d-5     /
      data bip1cs(  7) / -.4222669698 2681924884 7785158894 33 d-6     /
      data bip1cs(  8) / -.1914718629 8654689632 8354941812 77 d-6     /
      data bip1cs(  9) / +.6453106284 5583173611 0381578809 34 d-7     /
      data bip1cs( 10) / -.7844854677 1397719289 7483104486 28 d-8     /
      data bip1cs( 11) / -.9607721662 3785085879 1985335654 32 d-9     /
      data bip1cs( 12) / +.7000471331 6443966339 0060744020 68 d-9     /
      data bip1cs( 13) / -.1773178913 2814932022 0831280566 98 d-9     /
      data bip1cs( 14) / +.2272089478 3465236347 2821263893 11 d-10    /
      data bip1cs( 15) / +.1654045631 3972049847 0328606818 91 d-11    /
      data bip1cs( 16) / -.1851712555 9292316390 7553698966 93 d-11    /
      data bip1cs( 17) / +.5957631247 7117290165 6807155342 77 d-12    /
      data bip1cs( 18) / -.1219434814 7346564781 0557694989 86 d-12    /
      data bip1cs( 19) / +.1334786925 3513048815 3863478135 97 d-13    /
      data bip1cs( 20) / +.1727831152 4339746664 3847928897 31 d-14    /
      data bip1cs( 21) / -.1459073201 3016720735 2688717131 66 d-14    /
      data bip1cs( 22) / +.4901031992 7115819978 9949895201 04 d-15    /
      data bip1cs( 23) / -.1155654551 9261548129 2629727625 21 d-15    /
      data bip1cs( 24) / +.1909880736 7072411430 6717324415 24 d-16    /
      data bip1cs( 25) / -.1176896685 4492179886 9139959578 62 d-17    /
      data bip1cs( 26) / -.6327192514 9530064474 5374596770 47 d-18    /
      data bip1cs( 27) / +.3386183888 0715361614 1301913223 16 d-18    /
      data bip1cs( 28) / -.1072582532 1758625254 9921622196 22 d-18    /
      data bip1cs( 29) / +.2599570960 5617169284 7869331155 62 d-19    /
      data bip1cs( 30) / -.4847758357 1081193660 9623094941 01 d-20    /
      data bip1cs( 31) / +.5529891398 2121625361 5055131989 33 d-21    /
      data bip1cs( 32) / +.4942166082 6069471371 7481974442 66 d-22    /
      data bip1cs( 33) / -.5516212192 4145707458 0697208149 33 d-22    /
      data bip1cs( 34) / +.2143756041 7632550086 6318844996 26 d-22    /
      data bip1cs( 35) / -.6191031338 7655605798 7850611370 66 d-23    /
      data bip1cs( 36) / +.1462936270 7391245659 8309673369 59 d-23    /
      data bip1cs( 37) / -.2791848447 1059005576 1778660693 33 d-24    /
      data bip1cs( 38) / +.3645570316 8570246150 9067953493 33 d-25    /
      data bip1cs( 39) / +.5851182190 6188711839 3824597333 33 d-27    /
      data bip1cs( 40) / -.2494695048 7566510969 7450475519 99 d-26    /
      data bip1cs( 41) / +.1097932398 0338380977 9195794773 33 d-26    /
      data bip1cs( 42) / -.3474338834 5961115015 0340881066 66 d-27    /
      data bip1cs( 43) / +.9137340263 5349697363 1710822400 00 d-28    /
      data bip1cs( 44) / -.2051035272 8210629186 2477209599 99 d-28    /
      data bip1cs( 45) / +.3797698569 8546461748 6516223999 99 d-29    /
      data bip1cs( 46) / -.4847945849 7755565887 8484480000 00 d-30    /
      data bip1cs( 47) / -.1055830694 1230714314 2058666666 66 d-31    /
      data bip2cs(  1) / -.1135967375 8598867913 7973108955 27 d+0     /
      data bip2cs(  2) / +.4138147394 7881595760 0520811714 44 d-2     /
      data bip2cs(  3) / +.1353470622 1193329857 6969217275 08 d-3     /
      data bip2cs(  4) / +.1042731665 3015353405 8871834567 80 d-4     /
      data bip2cs(  5) / +.1347495476 7849907889 5899119589 25 d-5     /
      data bip2cs(  6) / +.1696537405 4383983356 0625111637 56 d-6     /
      data bip2cs(  7) / -.1009650086 5641624301 3662283963 73 d-7     /
      data bip2cs(  8) / -.1672911949 3778475127 8369730959 43 d-7     /
      data bip2cs(  9) / -.4581536448 5068383217 1527956133 91 d-8     /
      data bip2cs( 10) / +.3736681366 5655477274 0647493842 84 d-9     /
      data bip2cs( 11) / +.5766930320 1452448119 5846435021 11 d-9     /
      data bip2cs( 12) / +.6218126508 7850324095 3934087923 71 d-10    /
      data bip2cs( 13) / -.6329412028 2743068241 5891772813 54 d-10    /
      data bip2cs( 14) / -.1491504790 8598767633 9990919894 87 d-10    /
      data bip2cs( 15) / +.7889621394 2486771938 1723942948 91 d-11    /
      data bip2cs( 16) / +.2496051372 1857797984 8880640001 27 d-11    /
      data bip2cs( 17) / -.1213007528 7291659477 7466647348 14 d-11    /
      data bip2cs( 18) / -.3740493910 8727277887 3434604027 16 d-12    /
      data bip2cs( 19) / +.2237727814 0321476798 7834469310 91 d-12    /
      data bip2cs( 20) / +.4749029631 2192466341 9860774725 14 d-13    /
      data bip2cs( 21) / -.4526160799 1821224810 6056558312 94 d-13    /
      data bip2cs( 22) / -.3017227184 1986072645 1122458760 20 d-14    /
      data bip2cs( 23) / +.9105860355 8754058327 5926834789 08 d-14    /
      data bip2cs( 24) / -.9814923803 3807062926 6438642077 09 d-15    /
      data bip2cs( 25) / -.1642940064 7889465253 6012452515 89 d-14    /
      data bip2cs( 26) / +.5533483421 4274215451 1821146351 64 d-15    /
      data bip2cs( 27) / +.2175047986 4482655984 3743819981 56 d-15    /
      data bip2cs( 28) / -.1737923620 0220656971 2870295580 87 d-15    /
      data bip2cs( 29) / -.1047002347 1443714959 2839093136 04 d-17    /
      data bip2cs( 30) / +.3921914598 6056386925 4414033114 62 d-16    /
      data bip2cs( 31) / -.1162129368 6345196925 8240056659 10 d-16    /
      data bip2cs( 32) / -.5402747449 1754245533 7354113077 73 d-17    /
      data bip2cs( 33) / +.4544158212 3884610882 6754285533 04 d-17    /
      data bip2cs( 34) / -.2877559962 5221075729 4275854800 86 d-18    /
      data bip2cs( 35) / -.1001734092 7225341243 5961629604 40 d-17    /
      data bip2cs( 36) / +.4482393121 5068369856 3325619063 13 d-18    /
      data bip2cs( 37) / +.7613596865 4908942328 9489823667 75 d-19    /
      data bip2cs( 38) / -.1444832409 4881347238 9560601454 22 d-18    /
      data bip2cs( 39) / +.4046085944 9205362251 6248473921 12 d-19    /
      data bip2cs( 40) / +.2032108570 0338446891 3251907072 77 d-19    /
      data bip2cs( 41) / -.1960279547 1446798718 2727580419 62 d-19    /
      data bip2cs( 42) / +.3427303844 3944824263 5189582117 38 d-20    /
      data bip2cs( 43) / +.3702370585 3905135480 0246515931 54 d-20    /
      data bip2cs( 44) / -.2687959517 2041591131 4003329667 12 d-20    /
      data bip2cs( 45) / +.2812167846 3531712209 7144546833 64 d-21    /
      data bip2cs( 46) / +.6093396363 6177797173 2711196803 29 d-21    /
      data bip2cs( 47) / -.3866662189 7150844994 1729778934 13 d-21    /
      data bip2cs( 48) / +.2598933125 3566943450 8956519272 28 d-22    /
      data bip2cs( 49) / +.9719439362 2938503767 2811752160 84 d-22    /
      data bip2cs( 50) / -.5939281783 4375098415 6304782045 91 d-22    /
      data bip2cs( 51) / +.3886494997 7113015409 5919604394 44 d-23    /
      data bip2cs( 52) / +.1533430739 3617272869 7215128687 69 d-22    /
      data bip2cs( 53) / -.9751355520 9762624036 3365214097 24 d-23    /
      data bip2cs( 54) / +.9634064444 0489471424 7413393837 26 d-24    /
      data bip2cs( 55) / +.2384199940 0208880109 9467487924 54 d-23    /
      data bip2cs( 56) / -.1689698631 5019706184 8480442052 07 d-23    /
      data bip2cs( 57) / +.2735271588 8928361222 5784448014 78 d-24    /
      data bip2cs( 58) / +.3566001618 5409578960 1116850257 30 d-24    /
      data bip2cs( 59) / -.3023402660 8258827249 5342806669 54 d-24    /
      data bip2cs( 60) / +.7500204160 5973930653 1442048232 32 d-25    /
      data bip2cs( 61) / +.4840328757 5851388827 4553198387 48 d-25    /
      data bip2cs( 62) / -.5436413765 4447888432 6980102977 66 d-25    /
      data bip2cs( 63) / +.1928121447 0820962653 3459788097 56 d-25    /
      data bip2cs( 64) / +.5011635502 0532656659 6118141721 72 d-26    /
      data bip2cs( 65) / -.9504074458 2693253786 0346208699 72 d-26    /
      data bip2cs( 66) / +.4637264615 7101975948 6963322456 11 d-26    /
      data bip2cs( 67) / +.2117717070 4466954163 7681705770 46 d-28    /
      data bip2cs( 68) / -.1540485026 8168594303 6922045487 26 d-26    /
      data bip2cs( 69) / +.1038794429 3201213662 0478891944 41 d-26    /
      data bip2cs( 70) / -.1989007815 6915416751 3167282351 53 d-27    /
      data bip2cs( 71) / -.2102217387 8658495471 1770445225 32 d-27    /
      data bip2cs( 72) / +.2135309972 4525793150 6333566704 91 d-27    /
      data bip2cs( 73) / -.7904081074 7961342319 0235376326 27 d-28    /
      data bip2cs( 74) / -.1657535996 0435585049 9737417635 92 d-28    /
      data bip2cs( 75) / +.3886834285 0124112587 6255864965 37 d-28    /
      data bip2cs( 76) / -.2230923733 0896866182 6215624247 17 d-28    /
      data bip2cs( 77) / +.2777724442 0176260265 6259774043 82 d-29    /
      data bip2cs( 78) / +.5707854347 2657725368 7124337827 72 d-29    /
      data bip2cs( 79) / -.5174308444 5303852800 1733715552 80 d-29    /
      data bip2cs( 80) / +.1841328075 1095837198 4509270715 69 d-29    /
      data bip2cs( 81) / +.4442256239 0957094598 5440710686 47 d-30    /
      data bip2cs( 82) / -.9850414263 9629801547 4649582269 43 d-30    /
      data bip2cs( 83) / +.5885720135 3585104884 7541988819 95 d-30    /
      data bip2cs( 84) / -.9763607544 0429787961 4023126285 95 d-31    /
      data bip2cs( 85) / -.1358101199 6074695047 0635978841 22 d-30    /
      data bip2cs( 86) / +.1399974351 8492413270 5680483803 45 d-30    /
      data bip2cs( 87) / -.5975490454 5248477620 8845629811 18 d-31    /
      data bip2cs( 88) / -.4039165387 5428313641 0453275298 56 d-32    /
      data atr / 8.750690570 8484345088 0771988210 148 d0 /
      data btr / -2.093836321 3560543136 0096498526 268 d0 /
      data first /.true./
c***first executable statement  dbie
      if (first) then
         eta = 0.1*real(d1mach(3))
         nbif = initds (bifcs, 13, eta)
         nbig = initds (bigcs, 13, eta)
         nbif2 = initds (bif2cs, 15, eta)
         nbig2 = initds (big2cs, 15, eta)
         nbip1 = initds (bip1cs, 47, eta)
         nbip2 = initds (bip2cs, 88, eta)
c
         x3sml = eta**0.3333
         x32sml = 1.3104d0*x3sml**2
         xbig = d1mach(2)**0.6666d0
      endif
      first = .false.
c
      if (x.ge.(-1.0d0)) go to 20
      call d9aimp (x, xm, theta)
      dbie = xm * sin(theta)
      return
c
 20   if (x.gt.1.0d0) go to 30
      z = 0.d0
      if (abs(x).gt.x3sml) z = x**3
      dbie = 0.625d0 + dcsevl (z, bifcs, nbif) + x*(0.4375d0 +
     1  dcsevl (z, bigcs, nbig) )
      if (x.gt.x32sml) dbie = dbie * exp(-2.0d0*x*sqrt(x)/3.0d0)
      return
c
 30   if (x.gt.2.0d0) go to 40
      z = (2.0d0*x**3 - 9.0d0)/7.0d0
      dbie = exp(-2.0d0*x*sqrt(x)/3.0d0) * (1.125d0 +
     1  dcsevl (z, bif2cs, nbif2) + x*(0.625d0 + dcsevl (z, big2cs,
     2  nbig2)) )
      return
c
 40   if (x.gt.4.0d0) go to 50
      sqrtx = sqrt(x)
      z = atr/(x*sqrtx) + btr
      dbie = (0.625d0 + dcsevl (z, bip1cs, nbip1))/sqrt(sqrtx)
      return
c
 50   sqrtx = sqrt(x)
      z = -1.0d0
      if (x.lt.xbig) z = 16.d0/(x*sqrtx) - 1.0d0
      dbie = (0.625d0 + dcsevl (z, bip2cs, nbip2))/sqrt(sqrtx)
      return
c
      end
