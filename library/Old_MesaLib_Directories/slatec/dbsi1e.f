*deck dbsi1e
      double precision function dbsi1e (x)
c***begin prologue  dbsi1e
c***purpose  compute the exponentially scaled modified (hyperbolic)
c            bessel function of the first kind of order one.
c***library   slatec (fnlib)
c***category  c10b1
c***type      double precision (besi1e-s, dbsi1e-d)
c***keywords  exponentially scaled, first kind, fnlib,
c             hyperbolic bessel function, modified bessel function,
c             order one, special functions
c***author  fullerton, w., (lanl)
c***description
c
c dbsi1e(x) calculates the double precision exponentially scaled
c modified (hyperbolic) bessel function of the first kind of order
c one for double precision argument x.  the result is i1(x)
c multiplied by exp(-abs(x)).
c
c series for bi1        on the interval  0.          to  9.00000e+00
c                                        with weighted error   1.44e-32
c                                         log weighted error  31.84
c                               significant figures required  31.45
c                                    decimal places required  32.46
c
c series for ai1        on the interval  1.25000e-01 to  3.33333e-01
c                                        with weighted error   2.81e-32
c                                         log weighted error  31.55
c                               significant figures required  29.93
c                                    decimal places required  32.38
c
c series for ai12       on the interval  0.          to  1.25000e-01
c                                        with weighted error   1.83e-32
c                                         log weighted error  31.74
c                               significant figures required  29.97
c                                    decimal places required  32.66
c
c***references  (none)
c***routines called  d1mach, dcsevl, initds, xermsg
c***revision history  (yymmdd)
c   770701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c***end prologue  dbsi1e
      double precision x, bi1cs(17), ai1cs(46), ai12cs(69), xmin,
     1  xsml, y, d1mach, dcsevl
      logical first
      save bi1cs, ai1cs, ai12cs, nti1, ntai1, ntai12, xmin, xsml,
     1  first
      data bi1cs(  1) / -.1971713261 0998597316 1385032181 49 d-2     /
      data bi1cs(  2) / +.4073488766 7546480608 1553936520 14 d+0     /
      data bi1cs(  3) / +.3483899429 9959455866 2450377837 87 d-1     /
      data bi1cs(  4) / +.1545394556 3001236038 5984010584 89 d-2     /
      data bi1cs(  5) / +.4188852109 8377784129 4588320041 20 d-4     /
      data bi1cs(  6) / +.7649026764 8362114741 9597039660 69 d-6     /
      data bi1cs(  7) / +.1004249392 4741178689 1798080372 38 d-7     /
      data bi1cs(  8) / +.9932207791 9238106481 3712980548 63 d-10    /
      data bi1cs(  9) / +.7663801791 8447637275 2001716813 49 d-12    /
      data bi1cs( 10) / +.4741418923 8167394980 3880919481 60 d-14    /
      data bi1cs( 11) / +.2404114404 0745181799 8631720320 00 d-16    /
      data bi1cs( 12) / +.1017150500 7093713649 1211007999 99 d-18    /
      data bi1cs( 13) / +.3645093565 7866949458 4917333333 33 d-21    /
      data bi1cs( 14) / +.1120574950 2562039344 8106666666 66 d-23    /
      data bi1cs( 15) / +.2987544193 4468088832 0000000000 00 d-26    /
      data bi1cs( 16) / +.6973231093 9194709333 3333333333 33 d-29    /
      data bi1cs( 17) / +.1436794822 0620800000 0000000000 00 d-31    /
      data ai1cs(  1) / -.2846744181 8814786741 0037246830 7 d-1      /
      data ai1cs(  2) / -.1922953231 4432206510 4444877497 9 d-1      /
      data ai1cs(  3) / -.6115185857 9437889822 5624991778 5 d-3      /
      data ai1cs(  4) / -.2069971253 3502277088 8282377797 9 d-4      /
      data ai1cs(  5) / +.8585619145 8107255655 3694467313 8 d-5      /
      data ai1cs(  6) / +.1049498246 7115908625 1745399786 0 d-5      /
      data ai1cs(  7) / -.2918338918 4479022020 9343232669 7 d-6      /
      data ai1cs(  8) / -.1559378146 6317390001 6068096907 7 d-7      /
      data ai1cs(  9) / +.1318012367 1449447055 2530287390 9 d-7      /
      data ai1cs( 10) / -.1448423418 1830783176 3913446781 5 d-8      /
      data ai1cs( 11) / -.2908512243 9931420948 2504099301 0 d-9      /
      data ai1cs( 12) / +.1266388917 8753823873 1115969040 3 d-9      /
      data ai1cs( 13) / -.1664947772 9192206706 2417839858 0 d-10     /
      data ai1cs( 14) / -.1666653644 6094329760 9593715499 9 d-11     /
      data ai1cs( 15) / +.1242602414 2907682652 3216847201 7 d-11     /
      data ai1cs( 16) / -.2731549379 6724323972 5146142863 3 d-12     /
      data ai1cs( 17) / +.2023947881 6458037807 0026268898 1 d-13     /
      data ai1cs( 18) / +.7307950018 1168836361 9869812612 3 d-14     /
      data ai1cs( 19) / -.3332905634 4046749438 1377861713 3 d-14     /
      data ai1cs( 20) / +.7175346558 5129537435 4225466567 0 d-15     /
      data ai1cs( 21) / -.6982530324 7962563558 5062922365 6 d-16     /
      data ai1cs( 22) / -.1299944201 5627607600 6044608058 7 d-16     /
      data ai1cs( 23) / +.8120942864 2427988920 5467834286 0 d-17     /
      data ai1cs( 24) / -.2194016207 4107368981 5626664378 3 d-17     /
      data ai1cs( 25) / +.3630516170 0296548482 7986093233 4 d-18     /
      data ai1cs( 26) / -.1695139772 4391041663 0686679039 9 d-19     /
      data ai1cs( 27) / -.1288184829 8979078071 1688253822 2 d-19     /
      data ai1cs( 28) / +.5694428604 9670527801 0999107310 9 d-20     /
      data ai1cs( 29) / -.1459597009 0904800565 4550990028 7 d-20     /
      data ai1cs( 30) / +.2514546010 6757173140 8469133448 5 d-21     /
      data ai1cs( 31) / -.1844758883 1391248181 6040002901 3 d-22     /
      data ai1cs( 32) / -.6339760596 2279486419 2860979199 9 d-23     /
      data ai1cs( 33) / +.3461441102 0310111111 0814662656 0 d-23     /
      data ai1cs( 34) / -.1017062335 3713935475 9654102357 3 d-23     /
      data ai1cs( 35) / +.2149877147 0904314459 6250077866 6 d-24     /
      data ai1cs( 36) / -.3045252425 2386764017 4620617386 6 d-25     /
      data ai1cs( 37) / +.5238082144 7212859821 7763498666 6 d-27     /
      data ai1cs( 38) / +.1443583107 0893824464 1678950399 9 d-26     /
      data ai1cs( 39) / -.6121302074 8900427332 0067071999 9 d-27     /
      data ai1cs( 40) / +.1700011117 4678184183 4918980266 6 d-27     /
      data ai1cs( 41) / -.3596589107 9842441585 3521578666 6 d-28     /
      data ai1cs( 42) / +.5448178578 9484185766 5051306666 6 d-29     /
      data ai1cs( 43) / -.2731831789 6890849891 6256426666 6 d-30     /
      data ai1cs( 44) / -.1858905021 7086007157 7190399999 9 d-30     /
      data ai1cs( 45) / +.9212682974 5139334411 2776533333 3 d-31     /
      data ai1cs( 46) / -.2813835155 6535611063 7083306666 6 d-31     /
      data ai12cs(  1) / +.2857623501 8280120474 4984594846 9 d-1      /
      data ai12cs(  2) / -.9761097491 3614684077 6516445730 2 d-2      /
      data ai12cs(  3) / -.1105889387 6262371629 1256921277 5 d-3      /
      data ai12cs(  4) / -.3882564808 8776903934 5654477627 4 d-5      /
      data ai12cs(  5) / -.2512236237 8702089252 9452002212 1 d-6      /
      data ai12cs(  6) / -.2631468846 8895195068 3705236523 2 d-7      /
      data ai12cs(  7) / -.3835380385 9642370220 4500678796 8 d-8      /
      data ai12cs(  8) / -.5589743462 1965838068 6811252222 9 d-9      /
      data ai12cs(  9) / -.1897495812 3505412344 9892503323 8 d-10     /
      data ai12cs( 10) / +.3252603583 0154882385 5508067994 9 d-10     /
      data ai12cs( 11) / +.1412580743 6613781331 6336633284 6 d-10     /
      data ai12cs( 12) / +.2035628544 1470895072 2452613684 0 d-11     /
      data ai12cs( 13) / -.7198551776 2459085120 9258989044 6 d-12     /
      data ai12cs( 14) / -.4083551111 0921973182 2849963969 1 d-12     /
      data ai12cs( 15) / -.2101541842 7726643130 1984572746 2 d-13     /
      data ai12cs( 16) / +.4272440016 7119513542 9778833699 7 d-13     /
      data ai12cs( 17) / +.1042027698 4128802764 1741449994 8 d-13     /
      data ai12cs( 18) / -.3814403072 4370078047 6707253539 6 d-14     /
      data ai12cs( 19) / -.1880354775 5107824485 1273453396 3 d-14     /
      data ai12cs( 20) / +.3308202310 9209282827 3190335240 5 d-15     /
      data ai12cs( 21) / +.2962628997 6459501390 6854654205 2 d-15     /
      data ai12cs( 22) / -.3209525921 9934239587 7837353288 7 d-16     /
      data ai12cs( 23) / -.4650305368 4893583255 7128281897 9 d-16     /
      data ai12cs( 24) / +.4414348323 0717079499 4611375964 1 d-17     /
      data ai12cs( 25) / +.7517296310 8421048054 2545808029 5 d-17     /
      data ai12cs( 26) / -.9314178867 3268833756 8484784515 7 d-18     /
      data ai12cs( 27) / -.1242193275 1948909561 1678448869 7 d-17     /
      data ai12cs( 28) / +.2414276719 4548484690 0515390217 6 d-18     /
      data ai12cs( 29) / +.2026944384 0532851789 7192286069 2 d-18     /
      data ai12cs( 30) / -.6394267188 2690977870 4391988681 1 d-19     /
      data ai12cs( 31) / -.3049812452 3730958960 8488450357 1 d-19     /
      data ai12cs( 32) / +.1612841851 6514802251 3462230769 1 d-19     /
      data ai12cs( 33) / +.3560913964 3099250545 1027090462 0 d-20     /
      data ai12cs( 34) / -.3752017947 9364390796 6682800324 6 d-20     /
      data ai12cs( 35) / -.5787037427 0747993459 5198231074 1 d-22     /
      data ai12cs( 36) / +.7759997511 6481619619 8236963209 2 d-21     /
      data ai12cs( 37) / -.1452790897 2022333940 6445987408 5 d-21     /
      data ai12cs( 38) / -.1318225286 7390367021 2192275337 4 d-21     /
      data ai12cs( 39) / +.6116654862 9030707018 7999133171 7 d-22     /
      data ai12cs( 40) / +.1376279762 4271264277 3024338363 4 d-22     /
      data ai12cs( 41) / -.1690837689 9593478849 1983938230 6 d-22     /
      data ai12cs( 42) / +.1430596088 5954331539 8720108538 5 d-23     /
      data ai12cs( 43) / +.3409557828 0905940204 0536772990 2 d-23     /
      data ai12cs( 44) / -.1309457666 2707602278 4573872642 4 d-23     /
      data ai12cs( 45) / -.3940706411 2402574360 9352141755 7 d-24     /
      data ai12cs( 46) / +.4277137426 9808765808 0616679735 2 d-24     /
      data ai12cs( 47) / -.4424634830 9826068819 0028312302 9 d-25     /
      data ai12cs( 48) / -.8734113196 2307149721 1530978874 7 d-25     /
      data ai12cs( 49) / +.4045401335 6835333921 4340414242 8 d-25     /
      data ai12cs( 50) / +.7067100658 0946894656 5160771780 6 d-26     /
      data ai12cs( 51) / -.1249463344 5651052230 0286451860 5 d-25     /
      data ai12cs( 52) / +.2867392244 4034370329 7948339142 6 d-26     /
      data ai12cs( 53) / +.2044292892 5042926702 8177957421 0 d-26     /
      data ai12cs( 54) / -.1518636633 8204625683 7134680291 1 d-26     /
      data ai12cs( 55) / +.8110181098 1875758861 3227910703 7 d-28     /
      data ai12cs( 56) / +.3580379354 7735860911 2717370327 0 d-27     /
      data ai12cs( 57) / -.1692929018 9279025095 9305717544 8 d-27     /
      data ai12cs( 58) / -.2222902499 7024276390 6775852777 4 d-28     /
      data ai12cs( 59) / +.5424535127 1459696550 4860040112 8 d-28     /
      data ai12cs( 60) / -.1787068401 5780186887 6491299330 4 d-28     /
      data ai12cs( 61) / -.6565479068 7228149388 2392943788 0 d-29     /
      data ai12cs( 62) / +.7807013165 0611452809 2206770683 9 d-29     /
      data ai12cs( 63) / -.1816595260 6689797173 7933315222 1 d-29     /
      data ai12cs( 64) / -.1287704952 6600848203 7687559895 9 d-29     /
      data ai12cs( 65) / +.1114548172 9881645474 1370927369 4 d-29     /
      data ai12cs( 66) / -.1808343145 0393369391 5936887668 7 d-30     /
      data ai12cs( 67) / -.2231677718 2037719522 3244822893 9 d-30     /
      data ai12cs( 68) / +.1619029596 0803415106 1790980361 4 d-30     /
      data ai12cs( 69) / -.1834079908 8049414139 0130843921 0 d-31     /
      data first /.true./
c***first executable statement  dbsi1e
      if (first) then
         eta = 0.1*real(d1mach(3))
         nti1 = initds (bi1cs, 17, eta)
         ntai1 = initds (ai1cs, 46, eta)
         ntai12 = initds (ai12cs, 69, eta)
c
         xmin = 2.0d0*d1mach(1)
         xsml = sqrt(4.5d0*d1mach(3))
      endif
      first = .false.
c
      y = abs(x)
      if (y.gt.3.0d0) go to 20
c
      dbsi1e = 0.0d0
      if (y.eq.0.d0)  return
c
      if (y .le. xmin) call xermsg ('slatec', 'dbsi1e',
     +   'abs(x) so small i1 underflows', 1, 1)
      if (y.gt.xmin) dbsi1e = 0.5d0*x
      if (y.gt.xsml) dbsi1e = x*(0.875d0 + dcsevl (y*y/4.5d0-1.d0,
     1  bi1cs, nti1) )
      dbsi1e = exp(-y) * dbsi1e
      return
c
 20   if (y.le.8.d0) dbsi1e = (0.375d0 + dcsevl ((48.d0/y-11.d0)/5.d0,
     1  ai1cs, ntai1))/sqrt(y)
      if (y.gt.8.d0) dbsi1e = (0.375d0 + dcsevl (16.d0/y-1.d0, ai12cs,
     1  ntai12))/sqrt(y)
      dbsi1e = sign (dbsi1e, x)
c
      return
      end
