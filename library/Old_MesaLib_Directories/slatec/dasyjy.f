*deck dasyjy
      subroutine dasyjy (funjy, x, fnu, flgjy, in, y, wk, iflw)
c***begin prologue  dasyjy
c***subsidiary
c***purpose  subsidiary to dbesj and dbesy
c***library   slatec
c***type      double precision (asyjy-s, dasyjy-d)
c***author  amos, d. e., (snla)
c***description
c
c                 dasyjy computes bessel functions j and y
c               for arguments x.gt.0.0 and orders fnu .ge. 35.0
c               on flgjy = 1 and flgjy = -1 respectively
c
c                                  input
c
c      funjy - external subroutine jairy or yairy
c          x - argument, x.gt.0.0d0
c        fnu - order of the first bessel function
c      flgjy - selection flag
c              flgjy =  1.0d0 gives the j function
c              flgjy = -1.0d0 gives the y function
c         in - number of functions desired, in = 1 or 2
c
c                                  output
c
c         y  - a vector whose first in components contain the sequence
c       iflw - a flag indicating underflow or overflow
c                    return variables for besj only
c      wk(1) = 1 - (x/fnu)**2 = w**2
c      wk(2) = sqrt(abs(wk(1)))
c      wk(3) = abs(wk(2) - atan(wk(2)))  or
c              abs(ln((1 + wk(2))/(x/fnu)) - wk(2))
c            = abs((2/3)*zeta**(3/2))
c      wk(4) = fnu*wk(3)
c      wk(5) = (1.5*wk(3)*fnu)**(1/3) = sqrt(zeta)*fnu**(1/3)
c      wk(6) = sign(1.,w**2)*wk(5)**2 = sign(1.,w**2)*zeta*fnu**(2/3)
c      wk(7) = fnu**(1/3)
c
c     abstract   **** a double precision routine ****
c         dasyjy implements the uniform asymptotic expansion of
c         the j and y bessel functions for fnu.ge.35 and real
c         x.gt.0.0d0. the forms are identical except for a change
c         in sign of some of the terms. this change in sign is
c         accomplished by means of the flag flgjy = 1 or -1. on
c         flgjy = 1 the airy functions ai(x) and dai(x) are
c         supplied by the external function jairy, and on
c         flgjy = -1 the airy functions bi(x) and dbi(x) are
c         supplied by the external function yairy.
c
c***see also  dbesj, dbesy
c***routines called  d1mach, i1mach
c***revision history  (yymmdd)
c   750101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   891004  correction computation of elim.  (wrb)
c   891009  removed unreferenced variable.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910408  updated the author section.  (wrb)
c***end prologue  dasyjy
      integer i, iflw, in, j, jn,jr,ju,k, kb,klast,kmax,kp1, ks, ksp1,
     * kstemp, l, lr, lrp1, iseta, isetb
      integer i1mach
      double precision abw2, akm, alfa, alfa1, alfa2, ap, ar, asum, az,
     * beta, beta1, beta2, beta3, br, bsum, c, con1, con2,
     * con548,cr,crz32, dfi,elim, dr,fi, flgjy, fn, fnu,
     * fn2, gama, phi,  rcz, rden, relb, rfn2,  rtz, rzden,
     * sa, sb, suma, sumb, s1, ta, tau, tb, tfn, tol, tols, t2, upol,
     *  wk, x, xx, y, z, z32
      double precision d1mach
      dimension y(*), wk(*), c(65)
      dimension alfa(26,4), beta(26,5)
      dimension alfa1(26,2), alfa2(26,2)
      dimension beta1(26,2), beta2(26,2), beta3(26,1)
      dimension gama(26), kmax(5), ar(8), br(10), upol(10)
      dimension cr(10), dr(10)
      equivalence (alfa(1,1),alfa1(1,1))
      equivalence (alfa(1,3),alfa2(1,1))
      equivalence (beta(1,1),beta1(1,1))
      equivalence (beta(1,3),beta2(1,1))
      equivalence (beta(1,5),beta3(1,1))
      save tols, con1, con2, con548, ar, br, c,
     1 alfa1, alfa2, beta1, beta2, beta3, gama
      data tols            /-6.90775527898214d+00/
      data con1,con2,con548/
     1 6.66666666666667d-01, 3.33333333333333d-01, 1.04166666666667d-01/
      data  ar(1),  ar(2),  ar(3),  ar(4),  ar(5),  ar(6),  ar(7),
     a      ar(8)          / 8.35503472222222d-02, 1.28226574556327d-01,
     1 2.91849026464140d-01, 8.81627267443758d-01, 3.32140828186277d+00,
     2 1.49957629868626d+01, 7.89230130115865d+01, 4.74451538868264d+02/
      data  br(1), br(2), br(3), br(4), br(5), br(6), br(7), br(8),
     a      br(9), br(10)  /-1.45833333333333d-01,-9.87413194444444d-02,
     1-1.43312053915895d-01,-3.17227202678414d-01,-9.42429147957120d-01,
     2-3.51120304082635d+00,-1.57272636203680d+01,-8.22814390971859d+01,
     3-4.92355370523671d+02,-3.31621856854797d+03/
      data c(1), c(2), c(3), c(4), c(5), c(6), c(7), c(8), c(9), c(10),
     1     c(11), c(12), c(13), c(14), c(15), c(16), c(17), c(18),
     2     c(19), c(20), c(21), c(22), c(23), c(24)/
     3       -2.08333333333333d-01,        1.25000000000000d-01,
     4        3.34201388888889d-01,       -4.01041666666667d-01,
     5        7.03125000000000d-02,       -1.02581259645062d+00,
     6        1.84646267361111d+00,       -8.91210937500000d-01,
     7        7.32421875000000d-02,        4.66958442342625d+00,
     8       -1.12070026162230d+01,        8.78912353515625d+00,
     9       -2.36408691406250d+00,        1.12152099609375d-01,
     a       -2.82120725582002d+01,        8.46362176746007d+01,
     b       -9.18182415432400d+01,        4.25349987453885d+01,
     c       -7.36879435947963d+00,        2.27108001708984d-01,
     d        2.12570130039217d+02,       -7.65252468141182d+02,
     e        1.05999045252800d+03,       -6.99579627376133d+02/
      data c(25), c(26), c(27), c(28), c(29), c(30), c(31), c(32),
     1     c(33), c(34), c(35), c(36), c(37), c(38), c(39), c(40),
     2     c(41), c(42), c(43), c(44), c(45), c(46), c(47), c(48)/
     3        2.18190511744212d+02,       -2.64914304869516d+01,
     4        5.72501420974731d-01,       -1.91945766231841d+03,
     5        8.06172218173731d+03,       -1.35865500064341d+04,
     6        1.16553933368645d+04,       -5.30564697861340d+03,
     7        1.20090291321635d+03,       -1.08090919788395d+02,
     8        1.72772750258446d+00,        2.02042913309661d+04,
     9       -9.69805983886375d+04,        1.92547001232532d+05,
     a       -2.03400177280416d+05,        1.22200464983017d+05,
     b       -4.11926549688976d+04,        7.10951430248936d+03,
     c       -4.93915304773088d+02,        6.07404200127348d+00,
     d       -2.42919187900551d+05,        1.31176361466298d+06,
     e       -2.99801591853811d+06,        3.76327129765640d+06/
      data c(49), c(50), c(51), c(52), c(53), c(54), c(55), c(56),
     1     c(57), c(58), c(59), c(60), c(61), c(62), c(63), c(64),
     2     c(65)/
     3       -2.81356322658653d+06,        1.26836527332162d+06,
     4       -3.31645172484564d+05,        4.52187689813627d+04,
     5       -2.49983048181121d+03,        2.43805296995561d+01,
     6        3.28446985307204d+06,       -1.97068191184322d+07,
     7        5.09526024926646d+07,       -7.41051482115327d+07,
     8        6.63445122747290d+07,       -3.75671766607634d+07,
     9        1.32887671664218d+07,       -2.78561812808645d+06,
     a        3.08186404612662d+05,       -1.38860897537170d+04,
     b        1.10017140269247d+02/
      data alfa1(1,1), alfa1(2,1), alfa1(3,1), alfa1(4,1), alfa1(5,1),
     1     alfa1(6,1), alfa1(7,1), alfa1(8,1), alfa1(9,1), alfa1(10,1),
     2     alfa1(11,1),alfa1(12,1),alfa1(13,1),alfa1(14,1),alfa1(15,1),
     3     alfa1(16,1),alfa1(17,1),alfa1(18,1),alfa1(19,1),alfa1(20,1),
     4     alfa1(21,1),alfa1(22,1),alfa1(23,1),alfa1(24,1),alfa1(25,1),
     5     alfa1(26,1)     /-4.44444444444444d-03,-9.22077922077922d-04,
     6-8.84892884892885d-05, 1.65927687832450d-04, 2.46691372741793d-04,
     7 2.65995589346255d-04, 2.61824297061501d-04, 2.48730437344656d-04,
     8 2.32721040083232d-04, 2.16362485712365d-04, 2.00738858762752d-04,
     9 1.86267636637545d-04, 1.73060775917876d-04, 1.61091705929016d-04,
     1 1.50274774160908d-04, 1.40503497391270d-04, 1.31668816545923d-04,
     2 1.23667445598253d-04, 1.16405271474738d-04, 1.09798298372713d-04,
     3 1.03772410422993d-04, 9.82626078369363d-05, 9.32120517249503d-05,
     4 8.85710852478712d-05, 8.42963105715700d-05, 8.03497548407791d-05/
      data alfa1(1,2), alfa1(2,2), alfa1(3,2), alfa1(4,2), alfa1(5,2),
     1     alfa1(6,2), alfa1(7,2), alfa1(8,2), alfa1(9,2), alfa1(10,2),
     2     alfa1(11,2),alfa1(12,2),alfa1(13,2),alfa1(14,2),alfa1(15,2),
     3     alfa1(16,2),alfa1(17,2),alfa1(18,2),alfa1(19,2),alfa1(20,2),
     4     alfa1(21,2),alfa1(22,2),alfa1(23,2),alfa1(24,2),alfa1(25,2),
     5     alfa1(26,2)     / 6.93735541354589d-04, 2.32241745182922d-04,
     6-1.41986273556691d-05,-1.16444931672049d-04,-1.50803558053049d-04,
     7-1.55121924918096d-04,-1.46809756646466d-04,-1.33815503867491d-04,
     8-1.19744975684254d-04,-1.06184319207974d-04,-9.37699549891194d-05,
     9-8.26923045588193d-05,-7.29374348155221d-05,-6.44042357721016d-05,
     1-5.69611566009369d-05,-5.04731044303562d-05,-4.48134868008883d-05,
     2-3.98688727717599d-05,-3.55400532972042d-05,-3.17414256609022d-05,
     3-2.83996793904175d-05,-2.54522720634871d-05,-2.28459297164725d-05,
     4-2.05352753106481d-05,-1.84816217627666d-05,-1.66519330021394d-05/
      data alfa2(1,1), alfa2(2,1), alfa2(3,1), alfa2(4,1), alfa2(5,1),
     1     alfa2(6,1), alfa2(7,1), alfa2(8,1), alfa2(9,1), alfa2(10,1),
     2     alfa2(11,1),alfa2(12,1),alfa2(13,1),alfa2(14,1),alfa2(15,1),
     3     alfa2(16,1),alfa2(17,1),alfa2(18,1),alfa2(19,1),alfa2(20,1),
     4     alfa2(21,1),alfa2(22,1),alfa2(23,1),alfa2(24,1),alfa2(25,1),
     5     alfa2(26,1)     /-3.54211971457744d-04,-1.56161263945159d-04,
     6 3.04465503594936d-05, 1.30198655773243d-04, 1.67471106699712d-04,
     7 1.70222587683593d-04, 1.56501427608595d-04, 1.36339170977445d-04,
     8 1.14886692029825d-04, 9.45869093034688d-05, 7.64498419250898d-05,
     9 6.07570334965197d-05, 4.74394299290509d-05, 3.62757512005344d-05,
     1 2.69939714979225d-05, 1.93210938247939d-05, 1.30056674793963d-05,
     2 7.82620866744497d-06, 3.59257485819352d-06, 1.44040049814252d-07,
     3-2.65396769697939d-06,-4.91346867098486d-06,-6.72739296091248d-06,
     4-8.17269379678658d-06,-9.31304715093561d-06,-1.02011418798016d-05/
      data alfa2(1,2), alfa2(2,2), alfa2(3,2), alfa2(4,2), alfa2(5,2),
     1     alfa2(6,2), alfa2(7,2), alfa2(8,2), alfa2(9,2), alfa2(10,2),
     2     alfa2(11,2),alfa2(12,2),alfa2(13,2),alfa2(14,2),alfa2(15,2),
     3     alfa2(16,2),alfa2(17,2),alfa2(18,2),alfa2(19,2),alfa2(20,2),
     4     alfa2(21,2),alfa2(22,2),alfa2(23,2),alfa2(24,2),alfa2(25,2),
     5     alfa2(26,2)     / 3.78194199201773d-04, 2.02471952761816d-04,
     6-6.37938506318862d-05,-2.38598230603006d-04,-3.10916256027362d-04,
     7-3.13680115247576d-04,-2.78950273791323d-04,-2.28564082619141d-04,
     8-1.75245280340847d-04,-1.25544063060690d-04,-8.22982872820208d-05,
     9-4.62860730588116d-05,-1.72334302366962d-05, 5.60690482304602d-06,
     1 2.31395443148287d-05, 3.62642745856794d-05, 4.58006124490189d-05,
     2 5.24595294959114d-05, 5.68396208545815d-05, 5.94349820393104d-05,
     3 6.06478527578422d-05, 6.08023907788436d-05, 6.01577894539460d-05,
     4 5.89199657344698d-05, 5.72515823777593d-05, 5.52804375585853d-05/
      data beta1(1,1), beta1(2,1), beta1(3,1), beta1(4,1), beta1(5,1),
     1     beta1(6,1), beta1(7,1), beta1(8,1), beta1(9,1), beta1(10,1),
     2     beta1(11,1),beta1(12,1),beta1(13,1),beta1(14,1),beta1(15,1),
     3     beta1(16,1),beta1(17,1),beta1(18,1),beta1(19,1),beta1(20,1),
     4     beta1(21,1),beta1(22,1),beta1(23,1),beta1(24,1),beta1(25,1),
     5     beta1(26,1)     / 1.79988721413553d-02, 5.59964911064388d-03,
     6 2.88501402231133d-03, 1.80096606761054d-03, 1.24753110589199d-03,
     7 9.22878876572938d-04, 7.14430421727287d-04, 5.71787281789705d-04,
     8 4.69431007606482d-04, 3.93232835462917d-04, 3.34818889318298d-04,
     9 2.88952148495752d-04, 2.52211615549573d-04, 2.22280580798883d-04,
     1 1.97541838033063d-04, 1.76836855019718d-04, 1.59316899661821d-04,
     2 1.44347930197334d-04, 1.31448068119965d-04, 1.20245444949303d-04,
     3 1.10449144504599d-04, 1.01828770740567d-04, 9.41998224204238d-05,
     4 8.74130545753834d-05, 8.13466262162801d-05, 7.59002269646219d-05/
      data beta1(1,2), beta1(2,2), beta1(3,2), beta1(4,2), beta1(5,2),
     1     beta1(6,2), beta1(7,2), beta1(8,2), beta1(9,2), beta1(10,2),
     2     beta1(11,2),beta1(12,2),beta1(13,2),beta1(14,2),beta1(15,2),
     3     beta1(16,2),beta1(17,2),beta1(18,2),beta1(19,2),beta1(20,2),
     4     beta1(21,2),beta1(22,2),beta1(23,2),beta1(24,2),beta1(25,2),
     5     beta1(26,2)     /-1.49282953213429d-03,-8.78204709546389d-04,
     6-5.02916549572035d-04,-2.94822138512746d-04,-1.75463996970783d-04,
     7-1.04008550460816d-04,-5.96141953046458d-05,-3.12038929076098d-05,
     8-1.26089735980230d-05,-2.42892608575730d-07, 8.05996165414274d-06,
     9 1.36507009262147d-05, 1.73964125472926d-05, 1.98672978842134d-05,
     1 2.14463263790823d-05, 2.23954659232457d-05, 2.28967783814713d-05,
     2 2.30785389811178d-05, 2.30321976080909d-05, 2.28236073720349d-05,
     3 2.25005881105292d-05, 2.20981015361991d-05, 2.16418427448104d-05,
     4 2.11507649256221d-05, 2.06388749782171d-05, 2.01165241997082d-05/
      data beta2(1,1), beta2(2,1), beta2(3,1), beta2(4,1), beta2(5,1),
     1     beta2(6,1), beta2(7,1), beta2(8,1), beta2(9,1), beta2(10,1),
     2     beta2(11,1),beta2(12,1),beta2(13,1),beta2(14,1),beta2(15,1),
     3     beta2(16,1),beta2(17,1),beta2(18,1),beta2(19,1),beta2(20,1),
     4     beta2(21,1),beta2(22,1),beta2(23,1),beta2(24,1),beta2(25,1),
     5     beta2(26,1)     / 5.52213076721293d-04, 4.47932581552385d-04,
     6 2.79520653992021d-04, 1.52468156198447d-04, 6.93271105657044d-05,
     7 1.76258683069991d-05,-1.35744996343269d-05,-3.17972413350427d-05,
     8-4.18861861696693d-05,-4.69004889379141d-05,-4.87665447413787d-05,
     9-4.87010031186735d-05,-4.74755620890087d-05,-4.55813058138628d-05,
     1-4.33309644511266d-05,-4.09230193157750d-05,-3.84822638603221d-05,
     2-3.60857167535411d-05,-3.37793306123367d-05,-3.15888560772110d-05,
     3-2.95269561750807d-05,-2.75978914828336d-05,-2.58006174666884d-05,
     4-2.41308356761280d-05,-2.25823509518346d-05,-2.11479656768913d-05/
      data beta2(1,2), beta2(2,2), beta2(3,2), beta2(4,2), beta2(5,2),
     1     beta2(6,2), beta2(7,2), beta2(8,2), beta2(9,2), beta2(10,2),
     2     beta2(11,2),beta2(12,2),beta2(13,2),beta2(14,2),beta2(15,2),
     3     beta2(16,2),beta2(17,2),beta2(18,2),beta2(19,2),beta2(20,2),
     4     beta2(21,2),beta2(22,2),beta2(23,2),beta2(24,2),beta2(25,2),
     5     beta2(26,2)     /-4.74617796559960d-04,-4.77864567147321d-04,
     6-3.20390228067038d-04,-1.61105016119962d-04,-4.25778101285435d-05,
     7 3.44571294294968d-05, 7.97092684075675d-05, 1.03138236708272d-04,
     8 1.12466775262204d-04, 1.13103642108481d-04, 1.08651634848774d-04,
     9 1.01437951597662d-04, 9.29298396593364d-05, 8.40293133016090d-05,
     1 7.52727991349134d-05, 6.69632521975731d-05, 5.92564547323195d-05,
     2 5.22169308826976d-05, 4.58539485165361d-05, 4.01445513891487d-05,
     3 3.50481730031328d-05, 3.05157995034347d-05, 2.64956119950516d-05,
     4 2.29363633690998d-05, 1.97893056664022d-05, 1.70091984636413d-05/
      data beta3(1,1), beta3(2,1), beta3(3,1), beta3(4,1), beta3(5,1),
     1     beta3(6,1), beta3(7,1), beta3(8,1), beta3(9,1), beta3(10,1),
     2     beta3(11,1),beta3(12,1),beta3(13,1),beta3(14,1),beta3(15,1),
     3     beta3(16,1),beta3(17,1),beta3(18,1),beta3(19,1),beta3(20,1),
     4     beta3(21,1),beta3(22,1),beta3(23,1),beta3(24,1),beta3(25,1),
     5     beta3(26,1)     / 7.36465810572578d-04, 8.72790805146194d-04,
     6 6.22614862573135d-04, 2.85998154194304d-04, 3.84737672879366d-06,
     7-1.87906003636972d-04,-2.97603646594555d-04,-3.45998126832656d-04,
     8-3.53382470916038d-04,-3.35715635775049d-04,-3.04321124789040d-04,
     9-2.66722723047613d-04,-2.27654214122820d-04,-1.89922611854562d-04,
     1-1.55058918599094d-04,-1.23778240761874d-04,-9.62926147717644d-05,
     2-7.25178327714425d-05,-5.22070028895634d-05,-3.50347750511901d-05,
     3-2.06489761035552d-05,-8.70106096849767d-06, 1.13698686675100d-06,
     4 9.16426474122779d-06, 1.56477785428873d-05, 2.08223629482467d-05/
      data gama(1),   gama(2),   gama(3),   gama(4),   gama(5),
     1     gama(6),   gama(7),   gama(8),   gama(9),   gama(10),
     2     gama(11),  gama(12),  gama(13),  gama(14),  gama(15),
     3     gama(16),  gama(17),  gama(18),  gama(19),  gama(20),
     4     gama(21),  gama(22),  gama(23),  gama(24),  gama(25),
     5     gama(26)        / 6.29960524947437d-01, 2.51984209978975d-01,
     6 1.54790300415656d-01, 1.10713062416159d-01, 8.57309395527395d-02,
     7 6.97161316958684d-02, 5.86085671893714d-02, 5.04698873536311d-02,
     8 4.42600580689155d-02, 3.93720661543510d-02, 3.54283195924455d-02,
     9 3.21818857502098d-02, 2.94646240791158d-02, 2.71581677112934d-02,
     1 2.51768272973862d-02, 2.34570755306079d-02, 2.19508390134907d-02,
     2 2.06210828235646d-02, 1.94388240897881d-02, 1.83810633800683d-02,
     3 1.74293213231963d-02, 1.65685837786612d-02, 1.57865285987918d-02,
     4 1.50729501494096d-02, 1.44193250839955d-02, 1.38184805735342d-02/
c***first executable statement  dasyjy
      ta = d1mach(3)
      tol = max(ta,1.0d-15)
      tb = d1mach(5)
      ju = i1mach(15)
      if(flgjy.eq.1.0d0) go to 6
      jr = i1mach(14)
      elim = -2.303d0*tb*(ju+jr)
      go to 7
    6 continue
      elim = -2.303d0*(tb*ju+3.0d0)
    7 continue
      fn = fnu
      iflw = 0
      do 170 jn=1,in
        xx = x/fn
        wk(1) = 1.0d0 - xx*xx
        abw2 = abs(wk(1))
        wk(2) = sqrt(abw2)
        wk(7) = fn**con2
        if (abw2.gt.0.27750d0) go to 80
c
c     asymptotic expansion
c     cases near x=fn, abs(1.-(x/fn)**2).le.0.2775
c     coefficients of asymptotic expansion by series
c
c     zeta and truncation for a(zeta) and b(zeta) series
c
c     kmax is truncation index for a(zeta) and b(zeta) series=max(2,sa)
c
        sa = 0.0d0
        if (abw2.eq.0.0d0) go to 10
        sa = tols/log(abw2)
   10   sb = sa
        do 20 i=1,5
          akm = max(sa,2.0d0)
          kmax(i) = int(akm)
          sa = sa + sb
   20   continue
        kb = kmax(5)
        klast = kb - 1
        sa = gama(kb)
        do 30 k=1,klast
          kb = kb - 1
          sa = sa*wk(1) + gama(kb)
   30   continue
        z = wk(1)*sa
        az = abs(z)
        rtz = sqrt(az)
        wk(3) = con1*az*rtz
        wk(4) = wk(3)*fn
        wk(5) = rtz*wk(7)
        wk(6) = -wk(5)*wk(5)
        if(z.le.0.0d0) go to 35
        if(wk(4).gt.elim) go to 75
        wk(6) = -wk(6)
   35   continue
        phi = sqrt(sqrt(sa+sa+sa+sa))
c
c     b(zeta) for s=0
c
        kb = kmax(5)
        klast = kb - 1
        sb = beta(kb,1)
        do 40 k=1,klast
          kb = kb - 1
          sb = sb*wk(1) + beta(kb,1)
   40   continue
        ksp1 = 1
        fn2 = fn*fn
        rfn2 = 1.0d0/fn2
        rden = 1.0d0
        asum = 1.0d0
        relb = tol*abs(sb)
        bsum = sb
        do 60 ks=1,4
          ksp1 = ksp1 + 1
          rden = rden*rfn2
c
c     a(zeta) and b(zeta) for s=1,2,3,4
c
          kstemp = 5 - ks
          kb = kmax(kstemp)
          klast = kb - 1
          sa = alfa(kb,ks)
          sb = beta(kb,ksp1)
          do 50 k=1,klast
            kb = kb - 1
            sa = sa*wk(1) + alfa(kb,ks)
            sb = sb*wk(1) + beta(kb,ksp1)
   50     continue
          ta = sa*rden
          tb = sb*rden
          asum = asum + ta
          bsum = bsum + tb
          if (abs(ta).le.tol .and. abs(tb).le.relb) go to 70
   60   continue
   70   continue
        bsum = bsum/(fn*wk(7))
        go to 160
c
   75   continue
        iflw = 1
        return
c
   80   continue
        upol(1) = 1.0d0
        tau = 1.0d0/wk(2)
        t2 = 1.0d0/wk(1)
        if (wk(1).ge.0.0d0) go to 90
c
c     cases for (x/fn).gt.sqrt(1.2775)
c
        wk(3) = abs(wk(2)-atan(wk(2)))
        wk(4) = wk(3)*fn
        rcz = -con1/wk(4)
        z32 = 1.5d0*wk(3)
        rtz = z32**con2
        wk(5) = rtz*wk(7)
        wk(6) = -wk(5)*wk(5)
        go to 100
   90   continue
c
c     cases for (x/fn).lt.sqrt(0.7225)
c
        wk(3) = abs(log((1.0d0+wk(2))/xx)-wk(2))
        wk(4) = wk(3)*fn
        rcz = con1/wk(4)
        if(wk(4).gt.elim) go to 75
        z32 = 1.5d0*wk(3)
        rtz = z32**con2
        wk(7) = fn**con2
        wk(5) = rtz*wk(7)
        wk(6) = wk(5)*wk(5)
  100   continue
        phi = sqrt((rtz+rtz)*tau)
        tb = 1.0d0
        asum = 1.0d0
        tfn = tau/fn
        rden=1.0d0/fn
        rfn2=rden*rden
        rden=1.0d0
        upol(2) = (c(1)*t2+c(2))*tfn
        crz32 = con548*rcz
        bsum = upol(2) + crz32
        relb = tol*abs(bsum)
        ap = tfn
        ks = 0
        kp1 = 2
        rzden = rcz
        l = 2
        iseta=0
        isetb=0
        do 140 lr=2,8,2
c
c     compute two u polynomials for next a(zeta) and b(zeta)
c
          lrp1 = lr + 1
          do 120 k=lr,lrp1
            ks = ks + 1
            kp1 = kp1 + 1
            l = l + 1
            s1 = c(l)
            do 110 j=2,kp1
              l = l + 1
              s1 = s1*t2 + c(l)
  110       continue
            ap = ap*tfn
            upol(kp1) = ap*s1
            cr(ks) = br(ks)*rzden
            rzden = rzden*rcz
            dr(ks) = ar(ks)*rzden
  120     continue
          suma = upol(lrp1)
          sumb = upol(lr+2) + upol(lrp1)*crz32
          ju = lrp1
          do 130 jr=1,lr
            ju = ju - 1
            suma = suma + cr(jr)*upol(ju)
            sumb = sumb + dr(jr)*upol(ju)
  130     continue
          rden=rden*rfn2
          tb = -tb
          if (wk(1).gt.0.0d0) tb = abs(tb)
          if(rden.lt.tol) go to 131
          asum = asum + suma*tb
          bsum = bsum + sumb*tb
          go to 140
  131     if(iseta.eq.1) go to 132
          if(abs(suma).lt.tol) iseta=1
          asum=asum+suma*tb
  132     if(isetb.eq.1) go to 133
          if(abs(sumb).lt.relb) isetb=1
          bsum=bsum+sumb*tb
  133     if(iseta.eq.1 .and. isetb.eq.1) go to 150
  140   continue
  150   tb = wk(5)
        if (wk(1).gt.0.0d0) tb = -tb
        bsum = bsum/tb
c
  160   continue
        call funjy(wk(6), wk(5), wk(4), fi, dfi)
        ta=1.0d0/tol
        tb=d1mach(1)*ta*1.0d+3
        if(abs(fi).gt.tb) go to 165
        fi=fi*ta
        dfi=dfi*ta
        phi=phi*tol
  165   continue
        y(jn) = flgjy*phi*(fi*asum+dfi*bsum)/wk(7)
        fn = fn - flgjy
  170 continue
      return
      end
