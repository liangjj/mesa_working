*deck dbesi
      subroutine dbesi (x, alpha, kode, n, y, nz)
c***begin prologue  dbesi
c***purpose  compute an n member sequence of i bessel functions
c            i/sub(alpha+k-1)/(x), k=1,...,n or scaled bessel functions
c            exp(-x)*i/sub(alpha+k-1)/(x), k=1,...,n for nonnegative
c            alpha and x.
c***library   slatec
c***category  c10b3
c***type      double precision (besi-s, dbesi-d)
c***keywords  i bessel function, special functions
c***author  amos, d. e., (snla)
c           daniel, s. l., (snla)
c***description
c
c     abstract  **** a double precision routine ****
c         dbesi computes an n member sequence of i bessel functions
c         i/sub(alpha+k-1)/(x), k=1,...,n or scaled bessel functions
c         exp(-x)*i/sub(alpha+k-1)/(x), k=1,...,n for nonnegative alpha
c         and x.  a combination of the power series, the asymptotic
c         expansion for x to infinity, and the uniform asymptotic
c         expansion for nu to infinity are applied over subdivisions of
c         the (nu,x) plane.  for values not covered by one of these
c         formulae, the order is incremented by an integer so that one
c         of these formulae apply.  backward recursion is used to reduce
c         orders by integer values.  the asymptotic expansion for x to
c         infinity is used only when the entire sequence (specifically
c         the last member) lies within the region covered by the
c         expansion.  leading terms of these expansions are used to test
c         for over or underflow where appropriate.  if a sequence is
c         requested and the last member would underflow, the result is
c         set to zero and the next lower order tried, etc., until a
c         member comes on scale or all are set to zero.  an overflow
c         cannot occur with scaling.
c
c         the maximum number of significant digits obtainable
c         is the smaller of 14 and the number of digits carried in
c         double precision arithmetic.
c
c     description of arguments
c
c         input      x,alpha are double precision
c           x      - x .ge. 0.0d0
c           alpha  - order of first member of the sequence,
c                    alpha .ge. 0.0d0
c           kode   - a parameter to indicate the scaling option
c                    kode=1 returns
c                           y(k)=        i/sub(alpha+k-1)/(x),
c                                k=1,...,n
c                    kode=2 returns
c                           y(k)=exp(-x)*i/sub(alpha+k-1)/(x),
c                                k=1,...,n
c           n      - number of members in the sequence, n .ge. 1
c
c         output     y is double precision
c           y      - a vector whose first n components contain
c                    values for i/sub(alpha+k-1)/(x) or scaled
c                    values for exp(-x)*i/sub(alpha+k-1)/(x),
c                    k=1,...,n depending on kode
c           nz     - number of components of y set to zero due to
c                    underflow,
c                    nz=0   , normal return, computation completed
c                    nz .ne. 0, last nz components of y set to zero,
c                             y(k)=0.0d0, k=n-nz+1,...,n.
c
c     error conditions
c         improper input arguments - a fatal error
c         overflow with kode=1 - a fatal error
c         underflow - a non-fatal error(nz .ne. 0)
c
c***references  d. e. amos, s. l. daniel and m. k. weston, cdc 6600
c                 subroutines ibess and jbess for bessel functions
c                 i(nu,x) and j(nu,x), x .ge. 0, nu .ge. 0, acm
c                 transactions on mathematical software 3, (1977),
c                 pp. 76-92.
c               f. w. j. olver, tables of bessel functions of moderate
c                 or large orders, npl mathematical tables 6, her
c                 majesty's stationery office, london, 1962.
c***routines called  d1mach, dasyik, dlngam, i1mach, xermsg
c***revision history  (yymmdd)
c   750101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   890911  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dbesi
c
      integer i, ialp, in, inlim, is, i1, k, kk, km, kode, kt,
     1 n, nn, ns, nz
      integer i1mach
      double precision ain,ak,akm,alpha,ans,ap,arg,atol,tolln,dfn,
     1 dtm, dx, earg, elim, etx, flgik,fn, fnf, fni,fnp1,fnu,gln,ra,
     2 rttpi, s, sx, sxo2, s1, s2, t, ta, tb, temp, tfn, tm, tol,
     3 trx, t2, x, xo2, xo2l, y, z
      double precision d1mach, dlngam
      dimension y(*), temp(3)
      save rttpi, inlim
      data rttpi           / 3.98942280401433d-01/
      data inlim           /          80         /
c***first executable statement  dbesi
      nz = 0
      kt = 1
c     i1mach(15) replaces i1mach(12) in a double precision code
c     i1mach(14) replaces i1mach(11) in a double precision code
      ra = d1mach(3)
      tol = max(ra,1.0d-15)
      i1 = -i1mach(15)
      gln = d1mach(5)
      elim = 2.303d0*(i1*gln-3.0d0)
c     tolln = -ln(tol)
      i1 = i1mach(14)+1
      tolln = 2.303d0*gln*i1
      tolln = min(tolln,34.5388d0)
      if (n-1) 590, 10, 20
   10 kt = 2
   20 nn = n
      if (kode.lt.1 .or. kode.gt.2) go to 570
      if (x) 600, 30, 80
   30 if (alpha) 580, 40, 50
   40 y(1) = 1.0d0
      if (n.eq.1) return
      i1 = 2
      go to 60
   50 i1 = 1
   60 do 70 i=i1,n
        y(i) = 0.0d0
   70 continue
      return
   80 continue
      if (alpha.lt.0.0d0) go to 580
c
      ialp = int(alpha)
      fni = ialp + n - 1
      fnf = alpha - ialp
      dfn = fni + fnf
      fnu = dfn
      in = 0
      xo2 = x*0.5d0
      sxo2 = xo2*xo2
      etx = kode - 1
      sx = etx*x
c
c     decision tree for region where series, asymptotic expansion for x
c     to infinity and asymptotic expansion for nu to infinity are
c     applied.
c
      if (sxo2.le.(fnu+1.0d0)) go to 90
      if (x.le.12.0d0) go to 110
      fn = 0.55d0*fnu*fnu
      fn = max(17.0d0,fn)
      if (x.ge.fn) go to 430
      ans = max(36.0d0-fnu,0.0d0)
      ns = int(ans)
      fni = fni + ns
      dfn = fni + fnf
      fn = dfn
      is = kt
      km = n - 1 + ns
      if (km.gt.0) is = 3
      go to 120
   90 fn = fnu
      fnp1 = fn + 1.0d0
      xo2l = log(xo2)
      is = kt
      if (x.le.0.5d0) go to 230
      ns = 0
  100 fni = fni + ns
      dfn = fni + fnf
      fn = dfn
      fnp1 = fn + 1.0d0
      is = kt
      if (n-1+ns.gt.0) is = 3
      go to 230
  110 xo2l = log(xo2)
      ns = int(sxo2-fnu)
      go to 100
  120 continue
c
c     overflow test on uniform asymptotic expansion
c
      if (kode.eq.2) go to 130
      if (alpha.lt.1.0d0) go to 150
      z = x/alpha
      ra = sqrt(1.0d0+z*z)
      gln = log((1.0d0+ra)/z)
      t = ra*(1.0d0-etx) + etx/(z+ra)
      arg = alpha*(t-gln)
      if (arg.gt.elim) go to 610
      if (km.eq.0) go to 140
  130 continue
c
c     underflow test on uniform asymptotic expansion
c
      z = x/fn
      ra = sqrt(1.0d0+z*z)
      gln = log((1.0d0+ra)/z)
      t = ra*(1.0d0-etx) + etx/(z+ra)
      arg = fn*(t-gln)
  140 if (arg.lt.(-elim)) go to 280
      go to 190
  150 if (x.gt.elim) go to 610
      go to 130
c
c     uniform asymptotic expansion for nu to infinity
c
  160 if (km.ne.0) go to 170
      y(1) = temp(3)
      return
  170 temp(1) = temp(3)
      in = ns
      kt = 1
      i1 = 0
  180 continue
      is = 2
      fni = fni - 1.0d0
      dfn = fni + fnf
      fn = dfn
      if(i1.eq.2) go to 350
      z = x/fn
      ra = sqrt(1.0d0+z*z)
      gln = log((1.0d0+ra)/z)
      t = ra*(1.0d0-etx) + etx/(z+ra)
      arg = fn*(t-gln)
  190 continue
      i1 = abs(3-is)
      i1 = max(i1,1)
      flgik = 1.0d0
      call dasyik(x,fn,kode,flgik,ra,arg,i1,temp(is))
      go to (180, 350, 510), is
c
c     series for (x/2)**2.le.nu+1
c
  230 continue
      gln = dlngam(fnp1)
      arg = fn*xo2l - gln - sx
      if (arg.lt.(-elim)) go to 300
      earg = exp(arg)
  240 continue
      s = 1.0d0
      if (x.lt.tol) go to 260
      ak = 3.0d0
      t2 = 1.0d0
      t = 1.0d0
      s1 = fn
      do 250 k=1,17
        s2 = t2 + s1
        t = t*sxo2/s2
        s = s + t
        if (abs(t).lt.tol) go to 260
        t2 = t2 + ak
        ak = ak + 2.0d0
        s1 = s1 + fn
  250 continue
  260 continue
      temp(is) = s*earg
      go to (270, 350, 500), is
  270 earg = earg*fn/xo2
      fni = fni - 1.0d0
      dfn = fni + fnf
      fn = dfn
      is = 2
      go to 240
c
c     set underflow value and update parameters
c
  280 y(nn) = 0.0d0
      nn = nn - 1
      fni = fni - 1.0d0
      dfn = fni + fnf
      fn = dfn
      if (nn-1) 340, 290, 130
  290 kt = 2
      is = 2
      go to 130
  300 y(nn) = 0.0d0
      nn = nn - 1
      fnp1 = fn
      fni = fni - 1.0d0
      dfn = fni + fnf
      fn = dfn
      if (nn-1) 340, 310, 320
  310 kt = 2
      is = 2
  320 if (sxo2.le.fnp1) go to 330
      go to 130
  330 arg = arg - xo2l + log(fnp1)
      if (arg.lt.(-elim)) go to 300
      go to 230
  340 nz = n - nn
      return
c
c     backward recursion section
c
  350 continue
      nz = n - nn
  360 continue
      if(kt.eq.2) go to 420
      s1 = temp(1)
      s2 = temp(2)
      trx = 2.0d0/x
      dtm = fni
      tm = (dtm+fnf)*trx
      if (in.eq.0) go to 390
c     backward recur to index alpha+nn-1
      do 380 i=1,in
        s = s2
        s2 = tm*s2 + s1
        s1 = s
        dtm = dtm - 1.0d0
        tm = (dtm+fnf)*trx
  380 continue
      y(nn) = s1
      if (nn.eq.1) return
      y(nn-1) = s2
      if (nn.eq.2) return
      go to 400
  390 continue
c     backward recur from index alpha+nn-1 to alpha
      y(nn) = s1
      y(nn-1) = s2
      if (nn.eq.2) return
  400 k = nn + 1
      do 410 i=3,nn
        k = k - 1
        y(k-2) = tm*y(k-1) + y(k)
        dtm = dtm - 1.0d0
        tm = (dtm+fnf)*trx
  410 continue
      return
  420 y(1) = temp(2)
      return
c
c     asymptotic expansion for x to infinity
c
  430 continue
      earg = rttpi/sqrt(x)
      if (kode.eq.2) go to 440
      if (x.gt.elim) go to 610
      earg = earg*exp(x)
  440 etx = 8.0d0*x
      is = kt
      in = 0
      fn = fnu
  450 dx = fni + fni
      tm = 0.0d0
      if (fni.eq.0.0d0 .and. abs(fnf).lt.tol) go to 460
      tm = 4.0d0*fnf*(fni+fni+fnf)
  460 continue
      dtm = dx*dx
      s1 = etx
      trx = dtm - 1.0d0
      dx = -(trx+tm)/etx
      t = dx
      s = 1.0d0 + dx
      atol = tol*abs(s)
      s2 = 1.0d0
      ak = 8.0d0
      do 470 k=1,25
        s1 = s1 + etx
        s2 = s2 + ak
        dx = dtm - s2
        ap = dx + tm
        t = -t*ap/s1
        s = s + t
        if (abs(t).le.atol) go to 480
        ak = ak + 8.0d0
  470 continue
  480 temp(is) = s*earg
      if(is.eq.2) go to 360
      is = 2
      fni = fni - 1.0d0
      dfn = fni + fnf
      fn = dfn
      go to 450
c
c     backward recursion with normalization by
c     asymptotic expansion for nu to infinity or power series.
c
  500 continue
c     computation of last order for series normalization
      akm = max(3.0d0-fn,0.0d0)
      km = int(akm)
      tfn = fn + km
      ta = (gln+tfn-0.9189385332d0-0.0833333333d0/tfn)/(tfn+0.5d0)
      ta = xo2l - ta
      tb = -(1.0d0-1.0d0/tfn)/tfn
      ain = tolln/(-ta+sqrt(ta*ta-tolln*tb)) + 1.5d0
      in = int(ain)
      in = in + km
      go to 520
  510 continue
c     computation of last order for asymptotic expansion normalization
      t = 1.0d0/(fn*ra)
      ain = tolln/(gln+sqrt(gln*gln+t*tolln)) + 1.5d0
      in = int(ain)
      if (in.gt.inlim) go to 160
  520 continue
      trx = 2.0d0/x
      dtm = fni + in
      tm = (dtm+fnf)*trx
      ta = 0.0d0
      tb = tol
      kk = 1
  530 continue
c
c     backward recur unindexed
c
      do 540 i=1,in
        s = tb
        tb = tm*tb + ta
        ta = s
        dtm = dtm - 1.0d0
        tm = (dtm+fnf)*trx
  540 continue
c     normalization
      if (kk.ne.1) go to 550
      ta = (ta/tb)*temp(3)
      tb = temp(3)
      kk = 2
      in = ns
      if (ns.ne.0) go to 530
  550 y(nn) = tb
      nz = n - nn
      if (nn.eq.1) return
      tb = tm*tb + ta
      k = nn - 1
      y(k) = tb
      if (nn.eq.2) return
      dtm = dtm - 1.0d0
      tm = (dtm+fnf)*trx
      km = k - 1
c
c     backward recur indexed
c
      do 560 i=1,km
        y(k-1) = tm*y(k) + y(k+1)
        dtm = dtm - 1.0d0
        tm = (dtm+fnf)*trx
        k = k - 1
  560 continue
      return
c
c
c
  570 continue
      call xermsg ('slatec', 'dbesi',
     +   'scaling option, kode, not 1 or 2.', 2, 1)
      return
  580 continue
      call xermsg ('slatec', 'dbesi', 'order, alpha, less than zero.',
     +   2, 1)
      return
  590 continue
      call xermsg ('slatec', 'dbesi', 'n less than one.', 2, 1)
      return
  600 continue
      call xermsg ('slatec', 'dbesi', 'x less than zero.', 2, 1)
      return
  610 continue
      call xermsg ('slatec', 'dbesi',
     +   'overflow, x too large for kode = 1.', 6, 1)
      return
      end
