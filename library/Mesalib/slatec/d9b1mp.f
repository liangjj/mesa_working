*DECK D9B1MP
      SUBROUTINE D9B1MP (X, AMPL, THETA)
C***BEGIN PROLOGUE  D9B1MP
C***SUBSIDIARY
C***PURPOSE  Evaluate the modulus and phase for the J1 and Y1 Bessel
C            functions.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C10A1
C***TYPE      DOUBLE PRECISION (D9B1MP-D)
C***KEYWORDS  BESSEL FUNCTION, FNLIB, MODULUS, PHASE, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Evaluate the modulus and phase for the Bessel J1 and Y1 functions.
C
C Series for BM1        on the interval  1.56250E-02 to  6.25000E-02
C                                        with weighted error   4.91E-32
C                                         log weighted error  31.31
C                               significant figures required  30.04
C                                    decimal places required  32.09
C
C Series for BT12       on the interval  1.56250E-02 to  6.25000E-02
C                                        with weighted error   3.33E-32
C                                         log weighted error  31.48
C                               significant figures required  31.05
C                                    decimal places required  32.27
C
C Series for BM12       on the interval  0.          to  1.56250E-02
C                                        with weighted error   5.01E-32
C                                         log weighted error  31.30
C                               significant figures required  29.99
C                                    decimal places required  32.10
C
C Series for BTH1       on the interval  0.          to  1.56250E-02
C                                        with weighted error   2.82E-32
C                                         log weighted error  31.55
C                               significant figures required  31.12
C                                    decimal places required  32.37
C
C***SEE ALSO  DBESJ1, DBESY1
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900720  Routine changed from user-callable to subsidiary.  (WRB)
C   920618  Removed space from variable name and code restructured to
C           use IF-THEN-ELSE.  (RWC, WRB)
C***END PROLOGUE  D9B1MP
      DOUBLE PRECISION X, AMPL, THETA, BM1CS(37), BT12CS(39),
     1  BM12CS(40), BTH1CS(44), XMAX, PI4, Z, D1MACH, DCSEVL
      LOGICAL FIRST
      SAVE BM1CS, BT12CS, BTH1CS, BM12CS, PI4, NBM1, NBT12,
     1 NBM12, NBTH1, XMAX, FIRST
      DATA BM1CS(  1) / +.1069845452 6180630149 6998530853 8 D+0      /
      DATA BM1CS(  2) / +.3274915039 7159649007 2905514344 5 D-2      /
      DATA BM1CS(  3) / -.2987783266 8316985920 3044577793 8 D-4      /
      DATA BM1CS(  4) / +.8331237177 9919745313 9322266902 3 D-6      /
      DATA BM1CS(  5) / -.4112665690 3020073048 9638172549 8 D-7      /
      DATA BM1CS(  6) / +.2855344228 7892152207 1975766316 1 D-8      /
      DATA BM1CS(  7) / -.2485408305 4156238780 6002659605 5 D-9      /
      DATA BM1CS(  8) / +.2543393338 0725824427 4248439717 4 D-10     /
      DATA BM1CS(  9) / -.2941045772 8229675234 8975082790 9 D-11     /
      DATA BM1CS( 10) / +.3743392025 4939033092 6505615362 6 D-12     /
      DATA BM1CS( 11) / -.5149118293 8211672187 2054824352 7 D-13     /
      DATA BM1CS( 12) / +.7552535949 8651439080 3404076419 9 D-14     /
      DATA BM1CS( 13) / -.1169409706 8288464441 6629062246 4 D-14     /
      DATA BM1CS( 14) / +.1896562449 4347915717 2182460506 0 D-15     /
      DATA BM1CS( 15) / -.3201955368 6932864206 6477531639 4 D-16     /
      DATA BM1CS( 16) / +.5599548399 3162041144 8416990549 3 D-17     /
      DATA BM1CS( 17) / -.1010215894 7304324431 1939044454 4 D-17     /
      DATA BM1CS( 18) / +.1873844985 7275629833 0204271957 3 D-18     /
      DATA BM1CS( 19) / -.3563537470 3285802192 7430143999 9 D-19     /
      DATA BM1CS( 20) / +.6931283819 9712383304 2276351999 9 D-20     /
      DATA BM1CS( 21) / -.1376059453 4065001522 5140893013 3 D-20     /
      DATA BM1CS( 22) / +.2783430784 1070802205 9977932799 9 D-21     /
      DATA BM1CS( 23) / -.5727595364 3205616893 4866943999 9 D-22     /
      DATA BM1CS( 24) / +.1197361445 9188926725 3575679999 9 D-22     /
      DATA BM1CS( 25) / -.2539928509 8918719766 4144042666 6 D-23     /
      DATA BM1CS( 26) / +.5461378289 6572959730 6961919999 9 D-24     /
      DATA BM1CS( 27) / -.1189211341 7733202889 8628949333 3 D-24     /
      DATA BM1CS( 28) / +.2620150977 3400815949 5782400000 0 D-25     /
      DATA BM1CS( 29) / -.5836810774 2556859019 2093866666 6 D-26     /
      DATA BM1CS( 30) / +.1313743500 0805957734 2361599999 9 D-26     /
      DATA BM1CS( 31) / -.2985814622 5103803553 3277866666 6 D-27     /
      DATA BM1CS( 32) / +.6848390471 3346049376 2559999999 9 D-28     /
      DATA BM1CS( 33) / -.1584401568 2224767211 9296000000 0 D-28     /
      DATA BM1CS( 34) / +.3695641006 5709380543 0101333333 3 D-29     /
      DATA BM1CS( 35) / -.8687115921 1446682430 1226666666 6 D-30     /
      DATA BM1CS( 36) / +.2057080846 1587634629 2906666666 6 D-30     /
      DATA BM1CS( 37) / -.4905225761 1162255185 2373333333 3 D-31     /
      DATA BT12CS(  1) / +.7382386012 8742974662 6208397927 64 D+0     /
      DATA BT12CS(  2) / -.3336111317 4483906384 4701476811 89 D-2     /
      DATA BT12CS(  3) / +.6146345488 8046964698 5148994201 86 D-4     /
      DATA BT12CS(  4) / -.2402458516 1602374264 9776354695 68 D-5     /
      DATA BT12CS(  5) / +.1466355557 7509746153 2105919972 04 D-6     /
      DATA BT12CS(  6) / -.1184191730 5589180567 0051475049 83 D-7     /
      DATA BT12CS(  7) / +.1157419896 3919197052 1254663030 55 D-8     /
      DATA BT12CS(  8) / -.1300116112 9439187449 3660077945 71 D-9     /
      DATA BT12CS(  9) / +.1624539114 1361731937 7421662736 67 D-10    /
      DATA BT12CS( 10) / -.2208963682 1403188752 1554417701 28 D-11    /
      DATA BT12CS( 11) / +.3218030425 8553177090 4743586537 78 D-12    /
      DATA BT12CS( 12) / -.4965314793 2768480785 5520211353 81 D-13    /
      DATA BT12CS( 13) / +.8043890043 2847825985 5588826393 17 D-14    /
      DATA BT12CS( 14) / -.1358912131 0161291384 6947126822 82 D-14    /
      DATA BT12CS( 15) / +.2381050439 7147214869 6765296059 73 D-15    /
      DATA BT12CS( 16) / -.4308146636 3849106724 4712414207 99 D-16    /
      DATA BT12CS( 17) / +.8020254403 2771002434 9935125504 00 D-17    /
      DATA BT12CS( 18) / -.1531631064 2462311864 2300274687 99 D-17    /
      DATA BT12CS( 19) / +.2992860635 2715568924 0730405546 66 D-18    /
      DATA BT12CS( 20) / -.5970996465 8085443393 8156366506 66 D-19    /
      DATA BT12CS( 21) / +.1214028966 9415185024 1608526506 66 D-19    /
      DATA BT12CS( 22) / -.2511511469 6612948901 0069777066 66 D-20    /
      DATA BT12CS( 23) / +.5279056717 0328744850 7383807999 99 D-21    /
      DATA BT12CS( 24) / -.1126050922 7550498324 3611613866 66 D-21    /
      DATA BT12CS( 25) / +.2434827735 9576326659 6634624000 00 D-22    /
      DATA BT12CS( 26) / -.5331726123 6931800130 0384426666 66 D-23    /
      DATA BT12CS( 27) / +.1181361505 9707121039 2059903999 99 D-23    /
      DATA BT12CS( 28) / -.2646536828 3353523514 8567893333 33 D-24    /
      DATA BT12CS( 29) / +.5990339404 1361503945 5778133333 33 D-25    /
      DATA BT12CS( 30) / -.1369085463 0829503109 1363839999 99 D-25    /
      DATA BT12CS( 31) / +.3157679015 4380228326 4136533333 33 D-26    /
      DATA BT12CS( 32) / -.7345791508 2084356491 4005333333 33 D-27    /
      DATA BT12CS( 33) / +.1722808148 0722747930 7059200000 00 D-27    /
      DATA BT12CS( 34) / -.4071690796 1286507941 0688000000 00 D-28    /
      DATA BT12CS( 35) / +.9693474513 6779622700 3733333333 33 D-29    /
      DATA BT12CS( 36) / -.2323763633 7765716765 3546666666 66 D-29    /
      DATA BT12CS( 37) / +.5607451067 3522029406 8906666666 66 D-30    /
      DATA BT12CS( 38) / -.1361646539 1539005860 5226666666 66 D-30    /
      DATA BT12CS( 39) / +.3326310923 3894654388 9066666666 66 D-31    /
      DATA BM12CS(  1) / +.9807979156 2330500272 7209354693 7 D-1      /
      DATA BM12CS(  2) / +.1150961189 5046853061 7548348460 2 D-2      /
      DATA BM12CS(  3) / -.4312482164 3382054098 8935809773 2 D-5      /
      DATA BM12CS(  4) / +.5951839610 0888163078 1302980183 2 D-7      /
      DATA BM12CS(  5) / -.1704844019 8269098574 0070158647 8 D-8      /
      DATA BM12CS(  6) / +.7798265413 6111095086 5817382740 1 D-10     /
      DATA BM12CS(  7) / -.4958986126 7664158094 9175495186 5 D-11     /
      DATA BM12CS(  8) / +.4038432416 4211415168 3820226514 4 D-12     /
      DATA BM12CS(  9) / -.3993046163 7251754457 6548384664 5 D-13     /
      DATA BM12CS( 10) / +.4619886183 1189664943 1334243277 5 D-14     /
      DATA BM12CS( 11) / -.6089208019 0953833013 4547261933 3 D-15     /
      DATA BM12CS( 12) / +.8960930916 4338764821 5704804124 9 D-16     /
      DATA BM12CS( 13) / -.1449629423 9420231229 1651891892 5 D-16     /
      DATA BM12CS( 14) / +.2546463158 5377760561 6514964806 8 D-17     /
      DATA BM12CS( 15) / -.4809472874 6478364442 5926371862 0 D-18     /
      DATA BM12CS( 16) / +.9687684668 2925990490 8727583912 4 D-19     /
      DATA BM12CS( 17) / -.2067213372 2779660232 4503811755 1 D-19     /
      DATA BM12CS( 18) / +.4646651559 1503847318 0276780959 0 D-20     /
      DATA BM12CS( 19) / -.1094966128 8483341382 4135132833 9 D-20     /
      DATA BM12CS( 20) / +.2693892797 2886828609 0570761278 5 D-21     /
      DATA BM12CS( 21) / -.6894992910 9303744778 1897002685 7 D-22     /
      DATA BM12CS( 22) / +.1830268262 7520629098 9066855474 0 D-22     /
      DATA BM12CS( 23) / -.5025064246 3519164281 5611355322 4 D-23     /
      DATA BM12CS( 24) / +.1423545194 4548060396 3169363419 4 D-23     /
      DATA BM12CS( 25) / -.4152191203 6164503880 6888676980 1 D-24     /
      DATA BM12CS( 26) / +.1244609201 5039793258 8233007654 7 D-24     /
      DATA BM12CS( 27) / -.3827336370 5693042994 3191866128 6 D-25     /
      DATA BM12CS( 28) / +.1205591357 8156175353 7472398183 5 D-25     /
      DATA BM12CS( 29) / -.3884536246 3764880764 3185936112 4 D-26     /
      DATA BM12CS( 30) / +.1278689528 7204097219 0489528346 1 D-26     /
      DATA BM12CS( 31) / -.4295146689 4479462720 6193691591 2 D-27     /
      DATA BM12CS( 32) / +.1470689117 8290708864 5680270798 3 D-27     /
      DATA BM12CS( 33) / -.5128315665 1060731281 8037401779 6 D-28     /
      DATA BM12CS( 34) / +.1819509585 4711693854 8143737328 6 D-28     /
      DATA BM12CS( 35) / -.6563031314 8419808676 1863505037 3 D-29     /
      DATA BM12CS( 36) / +.2404898976 9199606531 9891487583 4 D-29     /
      DATA BM12CS( 37) / -.8945966744 6906124732 3495824297 9 D-30     /
      DATA BM12CS( 38) / +.3376085160 6572310266 3714897824 0 D-30     /
      DATA BM12CS( 39) / -.1291791454 6206563609 1309991696 6 D-30     /
      DATA BM12CS( 40) / +.5008634462 9588105206 8495150125 4 D-31     /
      DATA BTH1CS(  1) / +.7474995720 3587276055 4434839696 95 D+0     /
      DATA BTH1CS(  2) / -.1240077714 4651711252 5457775413 84 D-2     /
      DATA BTH1CS(  3) / +.9925244240 4424527376 6414976895 92 D-5     /
      DATA BTH1CS(  4) / -.2030369073 7159711052 4193753756 08 D-6     /
      DATA BTH1CS(  5) / +.7535961770 5690885712 1840175836 29 D-8     /
      DATA BTH1CS(  6) / -.4166161271 5343550107 6300238562 28 D-9     /
      DATA BTH1CS(  7) / +.3070161807 0834890481 2451020912 16 D-10    /
      DATA BTH1CS(  8) / -.2817849963 7605213992 3240088839 24 D-11    /
      DATA BTH1CS(  9) / +.3079069673 9040295476 0281468216 47 D-12    /
      DATA BTH1CS( 10) / -.3880330026 2803434112 7873475547 81 D-13    /
      DATA BTH1CS( 11) / +.5509603960 8630904934 5617262085 62 D-14    /
      DATA BTH1CS( 12) / -.8659006076 8383779940 1033989539 94 D-15    /
      DATA BTH1CS( 13) / +.1485604914 1536749003 4236890606 83 D-15    /
      DATA BTH1CS( 14) / -.2751952981 5904085805 3712121250 09 D-16    /
      DATA BTH1CS( 15) / +.5455079609 0481089625 0362236409 23 D-17    /
      DATA BTH1CS( 16) / -.1148653450 1983642749 5436310271 77 D-17    /
      DATA BTH1CS( 17) / +.2553521337 7973900223 1990525335 22 D-18    /
      DATA BTH1CS( 18) / -.5962149019 7413450395 7682879078 49 D-19    /
      DATA BTH1CS( 19) / +.1455662290 2372718620 2883020058 33 D-19    /
      DATA BTH1CS( 20) / -.3702218542 2450538201 5797760195 93 D-20    /
      DATA BTH1CS( 21) / +.9776307412 5345357664 1684345179 24 D-21    /
      DATA BTH1CS( 22) / -.2672682163 9668488468 7237753930 52 D-21    /
      DATA BTH1CS( 23) / +.7545330038 4983271794 0381906557 64 D-22    /
      DATA BTH1CS( 24) / -.2194789991 9802744897 8923833716 47 D-22    /
      DATA BTH1CS( 25) / +.6564839462 3955262178 9069998174 93 D-23    /
      DATA BTH1CS( 26) / -.2015560429 8370207570 7840768695 19 D-23    /
      DATA BTH1CS( 27) / +.6341776855 6776143492 1446671856 70 D-24    /
      DATA BTH1CS( 28) / -.2041927788 5337895634 8137699555 91 D-24    /
      DATA BTH1CS( 29) / +.6719146422 0720567486 6589800185 51 D-25    /
      DATA BTH1CS( 30) / -.2256907911 0207573595 7090036873 36 D-25    /
      DATA BTH1CS( 31) / +.7729771989 2989706370 9269598719 29 D-26    /
      DATA BTH1CS( 32) / -.2696744451 2294640913 2114240809 20 D-26    /
      DATA BTH1CS( 33) / +.9574934451 8502698072 2955219336 27 D-27    /
      DATA BTH1CS( 34) / -.3456916844 8890113000 1756808276 27 D-27    /
      DATA BTH1CS( 35) / +.1268123481 7398436504 2119862383 74 D-27    /
      DATA BTH1CS( 36) / -.4723253663 0722639860 4649937134 45 D-28    /
      DATA BTH1CS( 37) / +.1785000847 8186376177 8586197964 17 D-28    /
      DATA BTH1CS( 38) / -.6840436100 4510395406 2152235667 46 D-29    /
      DATA BTH1CS( 39) / +.2656602867 1720419358 2934226722 12 D-29    /
      DATA BTH1CS( 40) / -.1045040252 7914452917 7141614846 70 D-29    /
      DATA BTH1CS( 41) / +.4161829082 5377144306 8619171970 64 D-30    /
      DATA BTH1CS( 42) / -.1677163920 3643714856 5013478828 87 D-30    /
      DATA BTH1CS( 43) / +.6836199777 6664389173 5359280285 28 D-31    /
      DATA BTH1CS( 44) / -.2817224786 1233641166 7395746228 10 D-31    /
      DATA PI4 / 0.7853981633 9744830961 5660845819 876 D0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  D9B1MP
      IF (FIRST) THEN
         ETA = 0.1*REAL(D1MACH(3))
         NBM1 = INITDS (BM1CS, 37, ETA)
         NBT12 = INITDS (BT12CS, 39, ETA)
         NBM12 = INITDS (BM12CS, 40, ETA)
         NBTH1 = INITDS (BTH1CS, 44, ETA)
C
         XMAX = 1.0D0/D1MACH(4)
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LT. 4.0D0) THEN
         CALL XERMSG ('SLATEC', 'D9B1MP', 'X must be .GE. 4', 1, 2)
         AMPL = 0.0D0
         THETA = 0.0D0
      ELSE IF (X .LE. 8.0D0) THEN
         Z = (128.0D0/(X*X) - 5.0D0)/3.0D0
         AMPL = (0.75D0 + DCSEVL (Z, BM1CS, NBM1))/SQRT(X)
         THETA = X - 3.0D0*PI4 + DCSEVL (Z, BT12CS, NBT12)/X
      ELSE
         IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'D9B1MP',
     +      'No precision because X is too big', 2, 2)
C
         Z = 128.0D0/(X*X) - 1.0D0
         AMPL = (0.75D0 + DCSEVL (Z, BM12CS, NBM12))/SQRT(X)
         THETA = X - 3.0D0*PI4 + DCSEVL (Z, BTH1CS, NBTH1)/X
      ENDIF
      RETURN
      END
