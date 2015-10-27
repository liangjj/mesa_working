*deck daie
      double precision function daie (x)
c***begin prologue  daie
c***purpose  calculate the airy function for a negative argument and an
c            exponentially scaled airy function for a non-negative
c            argument.
c***library   slatec (fnlib)
c***category  c10d
c***type      double precision (aie-s, daie-d)
c***keywords  exponentially scaled airy function, fnlib,
c             special functions
c***author  fullerton, w., (lanl)
c***description
c
c daie(x) calculates the airy function or the exponentially scaled
c airy function depending on the value of the argument.  the function
c and argument are both double precision.
c
c evaluate ai(x) for x .le. 0.0 and ai(x)*exp(zeta) where
c zeta = 2/3 * x**(3/2)  for x .ge. 0.0
c
c series for aif        on the interval -1.00000e+00 to  1.00000e+00
c                                        with weighted error   8.37e-33
c                                         log weighted error  32.08
c                               significant figures required  30.87
c                                    decimal places required  32.63
c
c series for aig        on the interval -1.00000e+00 to  1.00000e+00
c                                        with weighted error   7.47e-34
c                                         log weighted error  33.13
c                               significant figures required  31.50
c                                    decimal places required  33.68
c
c series for aip1       on the interval  1.25000e-01 to  1.00000e+00
c                                        with weighted error   3.69e-32
c                                         log weighted error  31.43
c                               significant figures required  29.55
c                                    decimal places required  32.31
c
c series for aip2       on the interval  0.          to  1.25000e-01
c                                        with weighted error   3.48e-32
c                                         log weighted error  31.46
c                               significant figures required  28.74
c                                    decimal places required  32.24
c
c***references  (none)
c***routines called  d1mach, d9aimp, dcsevl, initds
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920618  removed space from variable names.  (rwc, wrb)
c***end prologue  daie
      double precision x, aifcs(13), aigcs(13), aip1cs(57), aip2cs(37),
     1  sqrtx, theta, xbig, xm, x3sml, x32sml, z, d1mach, dcsevl
      logical first
      save aifcs, aigcs, aip1cs, aip2cs, naif, naig, naip1,
     1 naip2, x3sml, x32sml, xbig, first
      data aifcs(  1) / -.3797135849 6669997496 1970894694 14 d-1     /
      data aifcs(  2) / +.5919188853 7263638574 3197280137 77 d-1     /
      data aifcs(  3) / +.9862928057 7279975365 6038910440 60 d-3     /
      data aifcs(  4) / +.6848843819 0765667554 8548301824 12 d-5     /
      data aifcs(  5) / +.2594202596 2194713019 4892790814 03 d-7     /
      data aifcs(  6) / +.6176612774 0813750329 4457496972 36 d-10    /
      data aifcs(  7) / +.1009245417 2466117901 4295562246 01 d-12    /
      data aifcs(  8) / +.1201479251 1179938141 2880332253 33 d-15    /
      data aifcs(  9) / +.1088294558 8716991878 5252954666 66 d-18    /
      data aifcs( 10) / +.7751377219 6684887039 2384000000 00 d-22    /
      data aifcs( 11) / +.4454811203 7175638391 4666666666 66 d-25    /
      data aifcs( 12) / +.2109284523 1692343466 6666666666 66 d-28    /
      data aifcs( 13) / +.8370173591 0741333333 3333333333 33 d-32    /
      data aigcs(  1) / +.1815236558 1161273011 5562099578 64 d-1     /
      data aigcs(  2) / +.2157256316 6010755534 0306388199 68 d-1     /
      data aigcs(  3) / +.2567835698 7483249659 0524280901 33 d-3     /
      data aigcs(  4) / +.1426521411 9792403898 8294969217 21 d-5     /
      data aigcs(  5) / +.4572114920 0180426070 4340975581 91 d-8     /
      data aigcs(  6) / +.9525170843 5647098607 3922788405 92 d-11    /
      data aigcs(  7) / +.1392563460 5771399051 1504206861 90 d-13    /
      data aigcs(  8) / +.1507099914 2762379592 3069911386 66 d-16    /
      data aigcs(  9) / +.1255914831 2567778822 7032053333 33 d-19    /
      data aigcs( 10) / +.8306307377 0821340343 8293333333 33 d-23    /
      data aigcs( 11) / +.4465753849 3718567445 3333333333 33 d-26    /
      data aigcs( 12) / +.1990085503 4518869333 3333333333 33 d-29    /
      data aigcs( 13) / +.7470288525 6533333333 3333333333 33 d-33    /
      data aip1cs(  1) / -.2146951858 9105384554 6086346777 8 d-1      /
      data aip1cs(  2) / -.7535382535 0433011662 1972086556 5 d-2      /
      data aip1cs(  3) / +.5971527949 0263808520 3538888199 4 d-3      /
      data aip1cs(  4) / -.7283251254 2076106485 0236829154 8 d-4      /
      data aip1cs(  5) / +.1110297130 7392996665 1738182114 0 d-4      /
      data aip1cs(  6) / -.1950386152 2844057103 4693031403 3 d-5      /
      data aip1cs(  7) / +.3786973885 1595151938 8531967005 7 d-6      /
      data aip1cs(  8) / -.7929675297 3509782790 3907287915 4 d-7      /
      data aip1cs(  9) / +.1762247638 6742560755 6842012220 2 d-7      /
      data aip1cs( 10) / -.4110767539 6671950450 2989659389 3 d-8      /
      data aip1cs( 11) / +.9984770057 8578922471 8341410754 4 d-9      /
      data aip1cs( 12) / -.2510093251 3871222113 4986773003 4 d-9      /
      data aip1cs( 13) / +.6500501929 8606954092 7203860172 5 d-10     /
      data aip1cs( 14) / -.1727818405 3936165154 7887710736 6 d-10     /
      data aip1cs( 15) / +.4699378842 8245125783 6229287230 7 d-11     /
      data aip1cs( 16) / -.1304675656 2977439144 9124124627 2 d-11     /
      data aip1cs( 17) / +.3689698478 4626788104 7394838228 2 d-12     /
      data aip1cs( 18) / -.1061087206 6468061736 5035967903 5 d-12     /
      data aip1cs( 19) / +.3098414384 8781874386 6021007011 0 d-13     /
      data aip1cs( 20) / -.9174908079 8241393078 3342354785 1 d-14     /
      data aip1cs( 21) / +.2752049140 3472108956 9357906227 1 d-14     /
      data aip1cs( 22) / -.8353750115 9220465580 9139330188 0 d-15     /
      data aip1cs( 23) / +.2563931129 3579349475 6863616861 2 d-15     /
      data aip1cs( 24) / -.7950633762 5988549832 7374728982 2 d-16     /
      data aip1cs( 25) / +.2489283634 6030699774 3728117564 4 d-16     /
      data aip1cs( 26) / -.7864326933 9287355696 6462622129 6 d-17     /
      data aip1cs( 27) / +.2505687311 4399756723 2447064501 9 d-17     /
      data aip1cs( 28) / -.8047420364 1639095245 3795868224 1 d-18     /
      data aip1cs( 29) / +.2604097118 9520539644 4340110439 2 d-18     /
      data aip1cs( 30) / -.8486954164 0564122594 8248883418 4 d-19     /
      data aip1cs( 31) / +.2784706882 1423378433 5942918602 7 d-19     /
      data aip1cs( 32) / -.9195858953 4986129136 8722415135 4 d-20     /
      data aip1cs( 33) / +.3055304318 3742387422 4766822558 3 d-20     /
      data aip1cs( 34) / -.1021035455 4794778759 0217704843 9 d-20     /
      data aip1cs( 35) / +.3431118190 7437578440 0055568083 6 d-21     /
      data aip1cs( 36) / -.1159129341 7977495133 7692246310 9 d-21     /
      data aip1cs( 37) / +.3935772844 2002556108 3626822915 4 d-22     /
      data aip1cs( 38) / -.1342880980 2967176119 5671898903 8 d-22     /
      data aip1cs( 39) / +.4603287883 5200027416 5919030531 4 d-23     /
      data aip1cs( 40) / -.1585043927 0040642278 1077249938 7 d-23     /
      data aip1cs( 41) / +.5481275667 7296759089 2552375500 8 d-24     /
      data aip1cs( 42) / -.1903349371 8550472590 6401794894 5 d-24     /
      data aip1cs( 43) / +.6635682302 3740087167 7761211596 8 d-25     /
      data aip1cs( 44) / -.2322311650 0263143079 7520098645 3 d-25     /
      data aip1cs( 45) / +.8157640113 4291793131 4274369535 9 d-26     /
      data aip1cs( 46) / -.2875824240 6329004900 5748992955 7 d-26     /
      data aip1cs( 47) / +.1017329450 9429014350 7971431901 8 d-26     /
      data aip1cs( 48) / -.3610879108 7422164465 7570349055 9 d-27     /
      data aip1cs( 49) / +.1285788540 3639934212 5664034269 8 d-27     /
      data aip1cs( 50) / -.4592901037 3785474251 6069302271 9 d-28     /
      data aip1cs( 51) / +.1645597033 8207137258 1210248533 3 d-28     /
      data aip1cs( 52) / -.5913421299 8435018420 8792027136 0 d-29     /
      data aip1cs( 53) / +.2131057006 6049933034 7936950954 6 d-29     /
      data aip1cs( 54) / -.7701158157 7875982169 8276174506 6 d-30     /
      data aip1cs( 55) / +.2790533307 9689304175 8178377728 0 d-30     /
      data aip1cs( 56) / -.1013807715 1112840064 5224136703 9 d-30     /
      data aip1cs( 57) / +.3692580158 7196240936 5828621653 3 d-31     /
      data aip2cs(  1) / -.1743144969 2937551339 0355844011 d-2        /
      data aip2cs(  2) / -.1678938543 2554167163 2190613480 d-2        /
      data aip2cs(  3) / +.3596534033 5216603588 5983858114 d-4        /
      data aip2cs(  4) / -.1380818602 7392283545 7399383100 d-5        /
      data aip2cs(  5) / +.7411228077 3150529884 8699095233 d-7        /
      data aip2cs(  6) / -.5002382039 0013301313 0422866325 d-8        /
      data aip2cs(  7) / +.4006939174 1718424067 5446866355 d-9        /
      data aip2cs(  8) / -.3673312427 9590504419 9318496207 d-10       /
      data aip2cs(  9) / +.3760344395 9237385243 9592002918 d-11       /
      data aip2cs( 10) / -.4223213327 1874753802 6564938968 d-12       /
      data aip2cs( 11) / +.5135094540 3365707091 9618754120 d-13       /
      data aip2cs( 12) / -.6690958503 9047759565 1681356676 d-14       /
      data aip2cs( 13) / +.9266675456 4129064823 9550724382 d-15       /
      data aip2cs( 14) / -.1355143824 1607057633 3397356591 d-15       /
      data aip2cs( 15) / +.2081154963 1283099529 9006549335 d-16       /
      data aip2cs( 16) / -.3341164991 5917685687 1277570256 d-17       /
      data aip2cs( 17) / +.5585785845 8592431686 8032946585 d-18       /
      data aip2cs( 18) / -.9692190401 5236524751 8658209109 d-19       /
      data aip2cs( 19) / +.1740457001 2889320646 5696557738 d-19       /
      data aip2cs( 20) / -.3226409797 3113040024 7846333098 d-20       /
      data aip2cs( 21) / +.6160744711 0662525853 3259618986 d-21       /
      data aip2cs( 22) / -.1209363479 8249005907 6420676266 d-21       /
      data aip2cs( 23) / +.2436327633 1013810826 1570095786 d-22       /
      data aip2cs( 24) / -.5029142214 9745746894 3403144533 d-23       /
      data aip2cs( 25) / +.1062241755 4363568949 5470626133 d-23       /
      data aip2cs( 26) / -.2292842848 9598924150 9856324266 d-24       /
      data aip2cs( 27) / +.5051817339 2950374498 6884778666 d-25       /
      data aip2cs( 28) / -.1134981237 1441240497 9793920000 d-25       /
      data aip2cs( 29) / +.2597655659 8560698069 8374144000 d-26       /
      data aip2cs( 30) / -.6051246215 4293950617 2231679999 d-27       /
      data aip2cs( 31) / +.1433597779 6677280072 0295253333 d-27       /
      data aip2cs( 32) / -.3451477570 6089998628 0721066666 d-28       /
      data aip2cs( 33) / +.8438751902 1364674042 7025066666 d-29       /
      data aip2cs( 34) / -.2093961422 9818816943 4453333333 d-29       /
      data aip2cs( 35) / +.5270088734 7894550318 2848000000 d-30       /
      data aip2cs( 36) / -.1344574330 1455338578 9030399999 d-30       /
      data aip2cs( 37) / +.3475709645 2660114734 0117333333 d-31       /
      data first /.true./
c***first executable statement  daie
      if (first) then
         eta = 0.1*real(d1mach(3))
         naif = initds (aifcs, 13, eta)
         naig = initds (aigcs, 13, eta)
         naip1 = initds (aip1cs, 57, eta)
         naip2 = initds (aip2cs, 37, eta)
c
         x3sml = eta**0.3333e0
         x32sml = 1.3104d0*x3sml**2
         xbig = d1mach(2)**0.6666d0
      endif
      first = .false.
c
      if (x.ge.(-1.0d0)) go to 20
      call d9aimp (x, xm, theta)
      daie = xm * cos(theta)
      return
c
 20   if (x.gt.1.0d0) go to 30
      z = 0.0d0
      if (abs(x).gt.x3sml) z = x**3
      daie = 0.375d0 + (dcsevl (z, aifcs, naif) - x*(0.25d0 +
     1  dcsevl (z, aigcs, naig)) )
      if (x.gt.x32sml) daie = daie * exp (2.0d0*x*sqrt(x)/3.0d0)
      return
c
 30   if (x.gt.4.0d0) go to 40
      sqrtx = sqrt(x)
      z = (16.d0/(x*sqrtx) - 9.d0)/7.d0
      daie = (0.28125d0 + dcsevl (z, aip1cs, naip1))/sqrt(sqrtx)
      return
c
 40   sqrtx = sqrt(x)
      z = -1.0d0
      if (x.lt.xbig) z = 16.0d0/(x*sqrtx) - 1.0d0
      daie = (0.28125d0 + dcsevl (z, aip2cs, naip2))/sqrt(sqrtx)
      return
c
      end
