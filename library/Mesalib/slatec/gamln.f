*deck gamln
      real function gamln (z, ierr)
c***begin prologue  gamln
c***subsidiary
c***purpose  compute the logarithm of the gamma function
c***library   slatec
c***category  c7a
c***type      single precision (gamln-s, dgamln-d)
c***keywords  logarithm of gamma function
c***author  amos, d. e., (snl)
c***description
c
c         gamln computes the natural log of the gamma function for
c         z.gt.0.  the asymptotic expansion is used to generate values
c         greater than zmin which are adjusted by the recursion
c         g(z+1)=z*g(z) for z.le.zmin.  the function was made as
c         portable as possible by computing zmin from the number of base
c         10 digits in a word, rln=max(-alog10(r1mach(4)),0.5e-18)
c         limited to 18 digits of (relative) accuracy.
c
c         since integer arguments are common, a table look up on 100
c         values is used for speed of execution.
c
c     description of arguments
c
c         input
c           z      - real argument, z.gt.0.0e0
c
c         output
c           gamln  - natural log of the gamma function at z
c           ierr   - error flag
c                    ierr=0, normal return, computation completed
c                    ierr=1, z.le.0.0e0,    no computation
c
c***references  computation of bessel functions of complex argument
c                 by d. e. amos, sand83-0083, may, 1983.
c***routines called  i1mach, r1mach
c***revision history  (yymmdd)
c   830501  date written
c   830501  revision date from version 3.2
c   910415  prologue converted to version 4.0 format.  (bab)
c   920128  category corrected.  (wrb)
c   921215  gamln defined for z negative.  (wrb)
c***end prologue  gamln
c
      integer i, i1m, k, mz, nz, ierr, i1mach
      real cf, con, fln, fz, gln, rln, s, tlg, trm, tst, t1, wdtol, z,
     * zdmy, zinc, zm, zmin, zp, zsq
      real r1mach
      dimension cf(22), gln(100)
c           lngamma(n), n=1,100
      data gln(1), gln(2), gln(3), gln(4), gln(5), gln(6), gln(7),
     1     gln(8), gln(9), gln(10), gln(11), gln(12), gln(13), gln(14),
     2     gln(15), gln(16), gln(17), gln(18), gln(19), gln(20),
     3     gln(21), gln(22)/
     4     0.00000000000000000e+00,     0.00000000000000000e+00,
     5     6.93147180559945309e-01,     1.79175946922805500e+00,
     6     3.17805383034794562e+00,     4.78749174278204599e+00,
     7     6.57925121201010100e+00,     8.52516136106541430e+00,
     8     1.06046029027452502e+01,     1.28018274800814696e+01,
     9     1.51044125730755153e+01,     1.75023078458738858e+01,
     a     1.99872144956618861e+01,     2.25521638531234229e+01,
     b     2.51912211827386815e+01,     2.78992713838408916e+01,
     c     3.06718601060806728e+01,     3.35050734501368889e+01,
     d     3.63954452080330536e+01,     3.93398841871994940e+01,
     e     4.23356164607534850e+01,     4.53801388984769080e+01/
      data gln(23), gln(24), gln(25), gln(26), gln(27), gln(28),
     1     gln(29), gln(30), gln(31), gln(32), gln(33), gln(34),
     2     gln(35), gln(36), gln(37), gln(38), gln(39), gln(40),
     3     gln(41), gln(42), gln(43), gln(44)/
     4     4.84711813518352239e+01,     5.16066755677643736e+01,
     5     5.47847293981123192e+01,     5.80036052229805199e+01,
     6     6.12617017610020020e+01,     6.45575386270063311e+01,
     7     6.78897431371815350e+01,     7.12570389671680090e+01,
     8     7.46582363488301644e+01,     7.80922235533153106e+01,
     9     8.15579594561150372e+01,     8.50544670175815174e+01,
     a     8.85808275421976788e+01,     9.21361756036870925e+01,
     b     9.57196945421432025e+01,     9.93306124547874269e+01,
     c     1.02968198614513813e+02,     1.06631760260643459e+02,
     d     1.10320639714757395e+02,     1.14034211781461703e+02,
     e     1.17771881399745072e+02,     1.21533081515438634e+02/
      data gln(45), gln(46), gln(47), gln(48), gln(49), gln(50),
     1     gln(51), gln(52), gln(53), gln(54), gln(55), gln(56),
     2     gln(57), gln(58), gln(59), gln(60), gln(61), gln(62),
     3     gln(63), gln(64), gln(65), gln(66)/
     4     1.25317271149356895e+02,     1.29123933639127215e+02,
     5     1.32952575035616310e+02,     1.36802722637326368e+02,
     6     1.40673923648234259e+02,     1.44565743946344886e+02,
     7     1.48477766951773032e+02,     1.52409592584497358e+02,
     8     1.56360836303078785e+02,     1.60331128216630907e+02,
     9     1.64320112263195181e+02,     1.68327445448427652e+02,
     a     1.72352797139162802e+02,     1.76395848406997352e+02,
     b     1.80456291417543771e+02,     1.84533828861449491e+02,
     c     1.88628173423671591e+02,     1.92739047287844902e+02,
     d     1.96866181672889994e+02,     2.01009316399281527e+02,
     e     2.05168199482641199e+02,     2.09342586752536836e+02/
      data gln(67), gln(68), gln(69), gln(70), gln(71), gln(72),
     1     gln(73), gln(74), gln(75), gln(76), gln(77), gln(78),
     2     gln(79), gln(80), gln(81), gln(82), gln(83), gln(84),
     3     gln(85), gln(86), gln(87), gln(88)/
     4     2.13532241494563261e+02,     2.17736934113954227e+02,
     5     2.21956441819130334e+02,     2.26190548323727593e+02,
     6     2.30439043565776952e+02,     2.34701723442818268e+02,
     7     2.38978389561834323e+02,     2.43268849002982714e+02,
     8     2.47572914096186884e+02,     2.51890402209723194e+02,
     9     2.56221135550009525e+02,     2.60564940971863209e+02,
     a     2.64921649798552801e+02,     2.69291097651019823e+02,
     b     2.73673124285693704e+02,     2.78067573440366143e+02,
     c     2.82474292687630396e+02,     2.86893133295426994e+02,
     d     2.91323950094270308e+02,     2.95766601350760624e+02,
     e     3.00220948647014132e+02,     3.04686856765668715e+02/
      data gln(89), gln(90), gln(91), gln(92), gln(93), gln(94),
     1     gln(95), gln(96), gln(97), gln(98), gln(99), gln(100)/
     2     3.09164193580146922e+02,     3.13652829949879062e+02,
     3     3.18152639620209327e+02,     3.22663499126726177e+02,
     4     3.27185287703775217e+02,     3.31717887196928473e+02,
     5     3.36261181979198477e+02,     3.40815058870799018e+02,
     6     3.45379407062266854e+02,     3.49954118040770237e+02,
     7     3.54539085519440809e+02,     3.59134205369575399e+02/
c             coefficients of asymptotic expansion
      data cf(1), cf(2), cf(3), cf(4), cf(5), cf(6), cf(7), cf(8),
     1     cf(9), cf(10), cf(11), cf(12), cf(13), cf(14), cf(15),
     2     cf(16), cf(17), cf(18), cf(19), cf(20), cf(21), cf(22)/
     3     8.33333333333333333e-02,    -2.77777777777777778e-03,
     4     7.93650793650793651e-04,    -5.95238095238095238e-04,
     5     8.41750841750841751e-04,    -1.91752691752691753e-03,
     6     6.41025641025641026e-03,    -2.95506535947712418e-02,
     7     1.79644372368830573e-01,    -1.39243221690590112e+00,
     8     1.34028640441683920e+01,    -1.56848284626002017e+02,
     9     2.19310333333333333e+03,    -3.61087712537249894e+04,
     a     6.91472268851313067e+05,    -1.52382215394074162e+07,
     b     3.82900751391414141e+08,    -1.08822660357843911e+10,
     c     3.47320283765002252e+11,    -1.23696021422692745e+13,
     d     4.88788064793079335e+14,    -2.13203339609193739e+16/
c
c             ln(2*pi)
      data con                    /     1.83787706640934548e+00/
c
c***first executable statement  gamln
      ierr=0
      if (z.le.0.0e0) go to 70
      if (z.gt.101.0e0) go to 10
      nz = z
      fz = z - nz
      if (fz.gt.0.0e0) go to 10
      if (nz.gt.100) go to 10
      gamln = gln(nz)
      return
   10 continue
      wdtol = r1mach(4)
      wdtol = max(wdtol,0.5e-18)
      i1m = i1mach(11)
      rln = r1mach(5)*i1m
      fln = min(rln,20.0e0)
      fln = max(fln,3.0e0)
      fln = fln - 3.0e0
      zm = 1.8000e0 + 0.3875e0*fln
      mz = zm + 1
      zmin = mz
      zdmy = z
      zinc = 0.0e0
      if (z.ge.zmin) go to 20
      zinc = zmin - nz
      zdmy = z + zinc
   20 continue
      zp = 1.0e0/zdmy
      t1 = cf(1)*zp
      s = t1
      if (zp.lt.wdtol) go to 40
      zsq = zp*zp
      tst = t1*wdtol
      do 30 k=2,22
        zp = zp*zsq
        trm = cf(k)*zp
        if (abs(trm).lt.tst) go to 40
        s = s + trm
   30 continue
   40 continue
      if (zinc.ne.0.0e0) go to 50
      tlg = alog(z)
      gamln = z*(tlg-1.0e0) + 0.5e0*(con-tlg) + s
      return
   50 continue
      zp = 1.0e0
      nz = zinc
      do 60 i=1,nz
        zp = zp*(z+(i-1))
   60 continue
      tlg = alog(zdmy)
      gamln = zdmy*(tlg-1.0e0) - alog(zp) + 0.5e0*(con-tlg) + s
      return
c
c
   70 continue
      gamln = r1mach(2)
      ierr=1
      return
      end
