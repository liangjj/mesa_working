*deck d9b1mp
      subroutine d9b1mp (x, ampl, theta)
c***begin prologue  d9b1mp
c***subsidiary
c***purpose  evaluate the modulus and phase for the j1 and y1 bessel
c            functions.
c***library   slatec (fnlib)
c***category  c10a1
c***type      double precision (d9b1mp-d)
c***keywords  bessel function, fnlib, modulus, phase, special functions
c***author  fullerton, w., (lanl)
c***description
c
c evaluate the modulus and phase for the bessel j1 and y1 functions.
c
c series for bm1        on the interval  1.56250e-02 to  6.25000e-02
c                                        with weighted error   4.91e-32
c                                         log weighted error  31.31
c                               significant figures required  30.04
c                                    decimal places required  32.09
c
c series for bt12       on the interval  1.56250e-02 to  6.25000e-02
c                                        with weighted error   3.33e-32
c                                         log weighted error  31.48
c                               significant figures required  31.05
c                                    decimal places required  32.27
c
c series for bm12       on the interval  0.          to  1.56250e-02
c                                        with weighted error   5.01e-32
c                                         log weighted error  31.30
c                               significant figures required  29.99
c                                    decimal places required  32.10
c
c series for bth1       on the interval  0.          to  1.56250e-02
c                                        with weighted error   2.82e-32
c                                         log weighted error  31.55
c                               significant figures required  31.12
c                                    decimal places required  32.37
c
c***see also  dbesj1, dbesy1
c***references  (none)
c***routines called  d1mach, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900720  routine changed from user-callable to subsidiary.  (wrb)
c   920618  removed space from variable name and code restructured to
c           use if-then-else.  (rwc, wrb)
c***end prologue  d9b1mp
      real*8 x, ampl, theta, bm1cs(37), bt12cs(39),
     1  bm12cs(40), bth1cs(44), xmax, pi4, z, r1mach, csevl
      logical first
      save bm1cs, bt12cs, bth1cs, bm12cs, pi4, nbm1, nbt12,
     1 nbm12, nbth1, xmax, first
      data bm1cs(  1) / +.1069845452 6180630149 6998530853 8 d+0      /
      data bm1cs(  2) / +.3274915039 7159649007 2905514344 5 d-2      /
      data bm1cs(  3) / -.2987783266 8316985920 3044577793 8 d-4      /
      data bm1cs(  4) / +.8331237177 9919745313 9322266902 3 d-6      /
      data bm1cs(  5) / -.4112665690 3020073048 9638172549 8 d-7      /
      data bm1cs(  6) / +.2855344228 7892152207 1975766316 1 d-8      /
      data bm1cs(  7) / -.2485408305 4156238780 6002659605 5 d-9      /
      data bm1cs(  8) / +.2543393338 0725824427 4248439717 4 d-10     /
      data bm1cs(  9) / -.2941045772 8229675234 8975082790 9 d-11     /
      data bm1cs( 10) / +.3743392025 4939033092 6505615362 6 d-12     /
      data bm1cs( 11) / -.5149118293 8211672187 2054824352 7 d-13     /
      data bm1cs( 12) / +.7552535949 8651439080 3404076419 9 d-14     /
      data bm1cs( 13) / -.1169409706 8288464441 6629062246 4 d-14     /
      data bm1cs( 14) / +.1896562449 4347915717 2182460506 0 d-15     /
      data bm1cs( 15) / -.3201955368 6932864206 6477531639 4 d-16     /
      data bm1cs( 16) / +.5599548399 3162041144 8416990549 3 d-17     /
      data bm1cs( 17) / -.1010215894 7304324431 1939044454 4 d-17     /
      data bm1cs( 18) / +.1873844985 7275629833 0204271957 3 d-18     /
      data bm1cs( 19) / -.3563537470 3285802192 7430143999 9 d-19     /
      data bm1cs( 20) / +.6931283819 9712383304 2276351999 9 d-20     /
      data bm1cs( 21) / -.1376059453 4065001522 5140893013 3 d-20     /
      data bm1cs( 22) / +.2783430784 1070802205 9977932799 9 d-21     /
      data bm1cs( 23) / -.5727595364 3205616893 4866943999 9 d-22     /
      data bm1cs( 24) / +.1197361445 9188926725 3575679999 9 d-22     /
      data bm1cs( 25) / -.2539928509 8918719766 4144042666 6 d-23     /
      data bm1cs( 26) / +.5461378289 6572959730 6961919999 9 d-24     /
      data bm1cs( 27) / -.1189211341 7733202889 8628949333 3 d-24     /
      data bm1cs( 28) / +.2620150977 3400815949 5782400000 0 d-25     /
      data bm1cs( 29) / -.5836810774 2556859019 2093866666 6 d-26     /
      data bm1cs( 30) / +.1313743500 0805957734 2361599999 9 d-26     /
      data bm1cs( 31) / -.2985814622 5103803553 3277866666 6 d-27     /
      data bm1cs( 32) / +.6848390471 3346049376 2559999999 9 d-28     /
      data bm1cs( 33) / -.1584401568 2224767211 9296000000 0 d-28     /
      data bm1cs( 34) / +.3695641006 5709380543 0101333333 3 d-29     /
      data bm1cs( 35) / -.8687115921 1446682430 1226666666 6 d-30     /
      data bm1cs( 36) / +.2057080846 1587634629 2906666666 6 d-30     /
      data bm1cs( 37) / -.4905225761 1162255185 2373333333 3 d-31     /
      data bt12cs(  1) / +.7382386012 8742974662 6208397927 64 d+0     /
      data bt12cs(  2) / -.3336111317 4483906384 4701476811 89 d-2     /
      data bt12cs(  3) / +.6146345488 8046964698 5148994201 86 d-4     /
      data bt12cs(  4) / -.2402458516 1602374264 9776354695 68 d-5     /
      data bt12cs(  5) / +.1466355557 7509746153 2105919972 04 d-6     /
      data bt12cs(  6) / -.1184191730 5589180567 0051475049 83 d-7     /
      data bt12cs(  7) / +.1157419896 3919197052 1254663030 55 d-8     /
      data bt12cs(  8) / -.1300116112 9439187449 3660077945 71 d-9     /
      data bt12cs(  9) / +.1624539114 1361731937 7421662736 67 d-10    /
      data bt12cs( 10) / -.2208963682 1403188752 1554417701 28 d-11    /
      data bt12cs( 11) / +.3218030425 8553177090 4743586537 78 d-12    /
      data bt12cs( 12) / -.4965314793 2768480785 5520211353 81 d-13    /
      data bt12cs( 13) / +.8043890043 2847825985 5588826393 17 d-14    /
      data bt12cs( 14) / -.1358912131 0161291384 6947126822 82 d-14    /
      data bt12cs( 15) / +.2381050439 7147214869 6765296059 73 d-15    /
      data bt12cs( 16) / -.4308146636 3849106724 4712414207 99 d-16    /
      data bt12cs( 17) / +.8020254403 2771002434 9935125504 00 d-17    /
      data bt12cs( 18) / -.1531631064 2462311864 2300274687 99 d-17    /
      data bt12cs( 19) / +.2992860635 2715568924 0730405546 66 d-18    /
      data bt12cs( 20) / -.5970996465 8085443393 8156366506 66 d-19    /
      data bt12cs( 21) / +.1214028966 9415185024 1608526506 66 d-19    /
      data bt12cs( 22) / -.2511511469 6612948901 0069777066 66 d-20    /
      data bt12cs( 23) / +.5279056717 0328744850 7383807999 99 d-21    /
      data bt12cs( 24) / -.1126050922 7550498324 3611613866 66 d-21    /
      data bt12cs( 25) / +.2434827735 9576326659 6634624000 00 d-22    /
      data bt12cs( 26) / -.5331726123 6931800130 0384426666 66 d-23    /
      data bt12cs( 27) / +.1181361505 9707121039 2059903999 99 d-23    /
      data bt12cs( 28) / -.2646536828 3353523514 8567893333 33 d-24    /
      data bt12cs( 29) / +.5990339404 1361503945 5778133333 33 d-25    /
      data bt12cs( 30) / -.1369085463 0829503109 1363839999 99 d-25    /
      data bt12cs( 31) / +.3157679015 4380228326 4136533333 33 d-26    /
      data bt12cs( 32) / -.7345791508 2084356491 4005333333 33 d-27    /
      data bt12cs( 33) / +.1722808148 0722747930 7059200000 00 d-27    /
      data bt12cs( 34) / -.4071690796 1286507941 0688000000 00 d-28    /
      data bt12cs( 35) / +.9693474513 6779622700 3733333333 33 d-29    /
      data bt12cs( 36) / -.2323763633 7765716765 3546666666 66 d-29    /
      data bt12cs( 37) / +.5607451067 3522029406 8906666666 66 d-30    /
      data bt12cs( 38) / -.1361646539 1539005860 5226666666 66 d-30    /
      data bt12cs( 39) / +.3326310923 3894654388 9066666666 66 d-31    /
      data bm12cs(  1) / +.9807979156 2330500272 7209354693 7 d-1      /
      data bm12cs(  2) / +.1150961189 5046853061 7548348460 2 d-2      /
      data bm12cs(  3) / -.4312482164 3382054098 8935809773 2 d-5      /
      data bm12cs(  4) / +.5951839610 0888163078 1302980183 2 d-7      /
      data bm12cs(  5) / -.1704844019 8269098574 0070158647 8 d-8      /
      data bm12cs(  6) / +.7798265413 6111095086 5817382740 1 d-10     /
      data bm12cs(  7) / -.4958986126 7664158094 9175495186 5 d-11     /
      data bm12cs(  8) / +.4038432416 4211415168 3820226514 4 d-12     /
      data bm12cs(  9) / -.3993046163 7251754457 6548384664 5 d-13     /
      data bm12cs( 10) / +.4619886183 1189664943 1334243277 5 d-14     /
      data bm12cs( 11) / -.6089208019 0953833013 4547261933 3 d-15     /
      data bm12cs( 12) / +.8960930916 4338764821 5704804124 9 d-16     /
      data bm12cs( 13) / -.1449629423 9420231229 1651891892 5 d-16     /
      data bm12cs( 14) / +.2546463158 5377760561 6514964806 8 d-17     /
      data bm12cs( 15) / -.4809472874 6478364442 5926371862 0 d-18     /
      data bm12cs( 16) / +.9687684668 2925990490 8727583912 4 d-19     /
      data bm12cs( 17) / -.2067213372 2779660232 4503811755 1 d-19     /
      data bm12cs( 18) / +.4646651559 1503847318 0276780959 0 d-20     /
      data bm12cs( 19) / -.1094966128 8483341382 4135132833 9 d-20     /
      data bm12cs( 20) / +.2693892797 2886828609 0570761278 5 d-21     /
      data bm12cs( 21) / -.6894992910 9303744778 1897002685 7 d-22     /
      data bm12cs( 22) / +.1830268262 7520629098 9066855474 0 d-22     /
      data bm12cs( 23) / -.5025064246 3519164281 5611355322 4 d-23     /
      data bm12cs( 24) / +.1423545194 4548060396 3169363419 4 d-23     /
      data bm12cs( 25) / -.4152191203 6164503880 6888676980 1 d-24     /
      data bm12cs( 26) / +.1244609201 5039793258 8233007654 7 d-24     /
      data bm12cs( 27) / -.3827336370 5693042994 3191866128 6 d-25     /
      data bm12cs( 28) / +.1205591357 8156175353 7472398183 5 d-25     /
      data bm12cs( 29) / -.3884536246 3764880764 3185936112 4 d-26     /
      data bm12cs( 30) / +.1278689528 7204097219 0489528346 1 d-26     /
      data bm12cs( 31) / -.4295146689 4479462720 6193691591 2 d-27     /
      data bm12cs( 32) / +.1470689117 8290708864 5680270798 3 d-27     /
      data bm12cs( 33) / -.5128315665 1060731281 8037401779 6 d-28     /
      data bm12cs( 34) / +.1819509585 4711693854 8143737328 6 d-28     /
      data bm12cs( 35) / -.6563031314 8419808676 1863505037 3 d-29     /
      data bm12cs( 36) / +.2404898976 9199606531 9891487583 4 d-29     /
      data bm12cs( 37) / -.8945966744 6906124732 3495824297 9 d-30     /
      data bm12cs( 38) / +.3376085160 6572310266 3714897824 0 d-30     /
      data bm12cs( 39) / -.1291791454 6206563609 1309991696 6 d-30     /
      data bm12cs( 40) / +.5008634462 9588105206 8495150125 4 d-31     /
      data bth1cs(  1) / +.7474995720 3587276055 4434839696 95 d+0     /
      data bth1cs(  2) / -.1240077714 4651711252 5457775413 84 d-2     /
      data bth1cs(  3) / +.9925244240 4424527376 6414976895 92 d-5     /
      data bth1cs(  4) / -.2030369073 7159711052 4193753756 08 d-6     /
      data bth1cs(  5) / +.7535961770 5690885712 1840175836 29 d-8     /
      data bth1cs(  6) / -.4166161271 5343550107 6300238562 28 d-9     /
      data bth1cs(  7) / +.3070161807 0834890481 2451020912 16 d-10    /
      data bth1cs(  8) / -.2817849963 7605213992 3240088839 24 d-11    /
      data bth1cs(  9) / +.3079069673 9040295476 0281468216 47 d-12    /
      data bth1cs( 10) / -.3880330026 2803434112 7873475547 81 d-13    /
      data bth1cs( 11) / +.5509603960 8630904934 5617262085 62 d-14    /
      data bth1cs( 12) / -.8659006076 8383779940 1033989539 94 d-15    /
      data bth1cs( 13) / +.1485604914 1536749003 4236890606 83 d-15    /
      data bth1cs( 14) / -.2751952981 5904085805 3712121250 09 d-16    /
      data bth1cs( 15) / +.5455079609 0481089625 0362236409 23 d-17    /
      data bth1cs( 16) / -.1148653450 1983642749 5436310271 77 d-17    /
      data bth1cs( 17) / +.2553521337 7973900223 1990525335 22 d-18    /
      data bth1cs( 18) / -.5962149019 7413450395 7682879078 49 d-19    /
      data bth1cs( 19) / +.1455662290 2372718620 2883020058 33 d-19    /
      data bth1cs( 20) / -.3702218542 2450538201 5797760195 93 d-20    /
      data bth1cs( 21) / +.9776307412 5345357664 1684345179 24 d-21    /
      data bth1cs( 22) / -.2672682163 9668488468 7237753930 52 d-21    /
      data bth1cs( 23) / +.7545330038 4983271794 0381906557 64 d-22    /
      data bth1cs( 24) / -.2194789991 9802744897 8923833716 47 d-22    /
      data bth1cs( 25) / +.6564839462 3955262178 9069998174 93 d-23    /
      data bth1cs( 26) / -.2015560429 8370207570 7840768695 19 d-23    /
      data bth1cs( 27) / +.6341776855 6776143492 1446671856 70 d-24    /
      data bth1cs( 28) / -.2041927788 5337895634 8137699555 91 d-24    /
      data bth1cs( 29) / +.6719146422 0720567486 6589800185 51 d-25    /
      data bth1cs( 30) / -.2256907911 0207573595 7090036873 36 d-25    /
      data bth1cs( 31) / +.7729771989 2989706370 9269598719 29 d-26    /
      data bth1cs( 32) / -.2696744451 2294640913 2114240809 20 d-26    /
      data bth1cs( 33) / +.9574934451 8502698072 2955219336 27 d-27    /
      data bth1cs( 34) / -.3456916844 8890113000 1756808276 27 d-27    /
      data bth1cs( 35) / +.1268123481 7398436504 2119862383 74 d-27    /
      data bth1cs( 36) / -.4723253663 0722639860 4649937134 45 d-28    /
      data bth1cs( 37) / +.1785000847 8186376177 8586197964 17 d-28    /
      data bth1cs( 38) / -.6840436100 4510395406 2152235667 46 d-29    /
      data bth1cs( 39) / +.2656602867 1720419358 2934226722 12 d-29    /
      data bth1cs( 40) / -.1045040252 7914452917 7141614846 70 d-29    /
      data bth1cs( 41) / +.4161829082 5377144306 8619171970 64 d-30    /
      data bth1cs( 42) / -.1677163920 3643714856 5013478828 87 d-30    /
      data bth1cs( 43) / +.6836199777 6664389173 5359280285 28 d-31    /
      data bth1cs( 44) / -.2817224786 1233641166 7395746228 10 d-31    /
      data pi4 / 0.7853981633 9744830961 5660845819 876 d0 /
      data first /.true./
c***first executable statement  d9b1mp
      if (first) then
         eta = 0.1*real(r1mach(3))
         nbm1 = initds (bm1cs, 37, eta)
         nbt12 = initds (bt12cs, 39, eta)
         nbm12 = initds (bm12cs, 40, eta)
         nbth1 = initds (bth1cs, 44, eta)
c
         xmax = 1.0d0/r1mach(4)
      endif
      first = .false.
c
      if (x .lt. 4.0d0) then
         call lnkerr ('d9b1mp:x must be .ge. 4')
         ampl = 0.0d0
         theta = 0.0d0
      else if (x .le. 8.0d0) then
         z = (128.0d0/(x*x) - 5.0d0)/3.0d0
         ampl = (0.75d0 + csevl (z, bm1cs, nbm1))/sqrt(x)
         theta = x - 3.0d0*pi4 + csevl (z, bt12cs, nbt12)/x
      else
         if (x .gt. xmax) call lnkerr ('d9b1mp:no precision x too big')
c
         z = 128.0d0/(x*x) - 1.0d0
         ampl = (0.75d0 + csevl (z, bm12cs, nbm12))/sqrt(x)
         theta = x - 3.0d0*pi4 + csevl (z, bth1cs, nbth1)/x
      endif
      return
      end
