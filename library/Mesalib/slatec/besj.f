*deck besj
      subroutine besj (x, alpha, n, y, nz)
c***begin prologue  besj
c***purpose  compute an n member sequence of j bessel functions
c            j/sub(alpha+k-1)/(x), k=1,...,n for non-negative alpha
c            and x.
c***library   slatec
c***category  c10a3
c***type      single precision (besj-s, dbesj-d)
c***keywords  j bessel function, special functions
c***author  amos, d. e., (snla)
c           daniel, s. l., (snla)
c           weston, m. k., (snla)
c***description
c
c     abstract
c         besj computes an n member sequence of j bessel functions
c         j/sub(alpha+k-1)/(x), k=1,...,n for non-negative alpha and x.
c         a combination of the power series, the asymptotic expansion
c         for x to infinity and the uniform asymptotic expansion for
c         nu to infinity are applied over subdivisions of the (nu,x)
c         plane.  for values of (nu,x) not covered by one of these
c         formulae, the order is incremented or decremented by integer
c         values into a region where one of the formulae apply. backward
c         recursion is applied to reduce orders by integer values except
c         where the entire sequence lies in the oscillatory region.  in
c         this case forward recursion is stable and values from the
c         asymptotic expansion for x to infinity start the recursion
c         when it is efficient to do so.  leading terms of the series
c         and uniform expansion are tested for underflow.  if a sequence
c         is requested and the last member would underflow, the result
c         is set to zero and the next lower order tried, etc., until a
c         member comes on scale or all members are set to zero.
c         overflow cannot occur.
c
c     description of arguments
c
c         input
c           x      - x .ge. 0.0e0
c           alpha  - order of first member of the sequence,
c                    alpha .ge. 0.0e0
c           n      - number of members in the sequence, n .ge. 1
c
c         output
c           y      - a vector whose first  n components contain
c                    values for j/sub(alpha+k-1)/(x), k=1,...,n
c           nz     - number of components of y set to zero due to
c                    underflow,
c                    nz=0   , normal return, computation completed
c                    nz .ne. 0, last nz components of y set to zero,
c                             y(k)=0.0e0, k=n-nz+1,...,n.
c
c     error conditions
c         improper input arguments - a fatal error
c         underflow  - a non-fatal error (nz .ne. 0)
c
c***references  d. e. amos, s. l. daniel and m. k. weston, cdc 6600
c                 subroutines ibess and jbess for bessel functions
c                 i(nu,x) and j(nu,x), x .ge. 0, nu .ge. 0, acm
c                 transactions on mathematical software 3, (1977),
c                 pp. 76-92.
c               f. w. j. olver, tables of bessel functions of moderate
c                 or large orders, npl mathematical tables 6, her
c                 majesty's stationery office, london, 1962.
c***routines called  alngam, asyjy, i1mach, jairy, r1mach, xermsg
c***revision history  (yymmdd)
c   750101  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  besj
      external jairy
      integer i,ialp,idalp,iflw,in,inlim,is,i1,i2,k,kk,km,kt,n,nn,
     1        ns,nz
      integer i1mach
      real       ak,akm,alpha,ans,ap,arg,coef,dalpha,dfn,dtm,earg,
     1           elim1,etx,fidal,flgjy,fn,fnf,fni,fnp1,fnu,fnulim,
     2           gln,pdf,pidt,pp,rden,relb,rttp,rtwo,rtx,rzden,
     3           s,sa,sb,sxo2,s1,s2,t,ta,tau,tb,temp,tfn,tm,tol,
     4           tolln,trx,tx,t1,t2,wk,x,xo2,xo2l,y,rtol,slim
      save rtwo, pdf, rttp, pidt, pp, inlim, fnulim
      real r1mach, alngam
      dimension y(*), temp(3), fnulim(2), pp(4), wk(7)
      data rtwo,pdf,rttp,pidt                    / 1.34839972492648e+00,
     1 7.85398163397448e-01, 7.97884560802865e-01, 1.57079632679490e+00/
      data  pp(1),  pp(2),  pp(3),  pp(4)        / 8.72909153935547e+00,
     1 2.65693932265030e-01, 1.24578576865586e-01, 7.70133747430388e-04/
      data inlim           /      150            /
      data fnulim(1), fnulim(2) /      100.0e0,     60.0e0     /
c***first executable statement  besj
      nz = 0
      kt = 1
      ns=0
c     i1mach(14) replaces i1mach(11) in a double precision code
c     i1mach(15) replaces i1mach(12) in a double precision code
      ta = r1mach(3)
      tol = max(ta,1.0e-15)
      i1 = i1mach(11) + 1
      i2 = i1mach(12)
      tb = r1mach(5)
      elim1 = -2.303e0*(i2*tb+3.0e0)
      rtol=1.0e0/tol
      slim=r1mach(1)*1.0e+3*rtol
c     tolln = -ln(tol)
      tolln = 2.303e0*tb*i1
      tolln = min(tolln,34.5388e0)
      if (n-1) 720, 10, 20
   10 kt = 2
   20 nn = n
      if (x) 730, 30, 80
   30 if (alpha) 710, 40, 50
   40 y(1) = 1.0e0
      if (n.eq.1) return
      i1 = 2
      go to 60
   50 i1 = 1
   60 do 70 i=i1,n
        y(i) = 0.0e0
   70 continue
      return
   80 continue
      if (alpha.lt.0.0e0) go to 710
c
      ialp = int(alpha)
      fni = ialp + n - 1
      fnf = alpha - ialp
      dfn = fni + fnf
      fnu = dfn
      xo2 = x*0.5e0
      sxo2 = xo2*xo2
c
c     decision tree for region where series, asymptotic expansion for x
c     to infinity and asymptotic expansion for nu to infinity are
c     applied.
c
      if (sxo2.le.(fnu+1.0e0)) go to 90
      ta = max(20.0e0,fnu)
      if (x.gt.ta) go to 120
      if (x.gt.12.0e0) go to 110
      xo2l = log(xo2)
      ns = int(sxo2-fnu) + 1
      go to 100
   90 fn = fnu
      fnp1 = fn + 1.0e0
      xo2l = log(xo2)
      is = kt
      if (x.le.0.50e0) go to 330
      ns = 0
  100 fni = fni + ns
      dfn = fni + fnf
      fn = dfn
      fnp1 = fn + 1.0e0
      is = kt
      if (n-1+ns.gt.0) is = 3
      go to 330
  110 ans = max(36.0e0-fnu,0.0e0)
      ns = int(ans)
      fni = fni + ns
      dfn = fni + fnf
      fn = dfn
      is = kt
      if (n-1+ns.gt.0) is = 3
      go to 130
  120 continue
      rtx = sqrt(x)
      tau = rtwo*rtx
      ta = tau + fnulim(kt)
      if (fnu.le.ta) go to 480
      fn = fnu
      is = kt
c
c     uniform asymptotic expansion for nu to infinity
c
  130 continue
      i1 = abs(3-is)
      i1 = max(i1,1)
      flgjy = 1.0e0
      call asyjy(jairy,x,fn,flgjy,i1,temp(is),wk,iflw)
      if(iflw.ne.0) go to 380
      go to (320, 450, 620), is
  310 temp(1) = temp(3)
      kt = 1
  320 is = 2
      fni = fni - 1.0e0
      dfn = fni + fnf
      fn = dfn
      if(i1.eq.2) go to 450
      go to 130
c
c     series for (x/2)**2.le.nu+1
c
  330 continue
      gln = alngam(fnp1)
      arg = fn*xo2l - gln
      if (arg.lt.(-elim1)) go to 400
      earg = exp(arg)
  340 continue
      s = 1.0e0
      if (x.lt.tol) go to 360
      ak = 3.0e0
      t2 = 1.0e0
      t = 1.0e0
      s1 = fn
      do 350 k=1,17
        s2 = t2 + s1
        t = -t*sxo2/s2
        s = s + t
        if (abs(t).lt.tol) go to 360
        t2 = t2 + ak
        ak = ak + 2.0e0
        s1 = s1 + fn
  350 continue
  360 continue
      temp(is) = s*earg
      go to (370, 450, 610), is
  370 earg = earg*fn/xo2
      fni = fni - 1.0e0
      dfn = fni + fnf
      fn = dfn
      is = 2
      go to 340
c
c     set underflow value and update parameters
c     underflow can only occur for ns=0 since the order must be
c     larger than 36. therefore, ns need not be considered.
c
  380 y(nn) = 0.0e0
      nn = nn - 1
      fni = fni - 1.0e0
      dfn = fni + fnf
      fn = dfn
      if (nn-1) 440, 390, 130
  390 kt = 2
      is = 2
      go to 130
  400 y(nn) = 0.0e0
      nn = nn - 1
      fnp1 = fn
      fni = fni - 1.0e0
      dfn = fni + fnf
      fn = dfn
      if (nn-1) 440, 410, 420
  410 kt = 2
      is = 2
  420 if (sxo2.le.fnp1) go to 430
      go to 130
  430 arg = arg - xo2l + log(fnp1)
      if (arg.lt.(-elim1)) go to 400
      go to 330
  440 nz = n - nn
      return
c
c     backward recursion section
c
  450 continue
      if(ns.ne.0) go to 451
      nz = n - nn
      if (kt.eq.2) go to 470
c     backward recur from index alpha+nn-1 to alpha
      y(nn) = temp(1)
      y(nn-1) = temp(2)
      if (nn.eq.2) return
  451 continue
      trx = 2.0e0/x
      dtm = fni
      tm = (dtm+fnf)*trx
      ak=1.0e0
      ta=temp(1)
      tb=temp(2)
      if(abs(ta).gt.slim) go to 455
      ta=ta*rtol
      tb=tb*rtol
      ak=tol
  455 continue
      kk=2
      in=ns-1
      if(in.eq.0) go to 690
      if(ns.ne.0) go to 670
      k=nn-2
      do 460 i=3,nn
        s=tb
        tb=tm*tb-ta
        ta=s
        y(k)=tb*ak
        k=k-1
        dtm = dtm - 1.0e0
        tm = (dtm+fnf)*trx
  460 continue
      return
  470 y(1) = temp(2)
      return
c
c     asymptotic expansion for x to infinity with forward recursion in
c     oscillatory region x.gt.max(20, nu), provided the last member
c     of the sequence is also in the region.
c
  480 continue
      in = int(alpha-tau+2.0e0)
      if (in.le.0) go to 490
      idalp = ialp - in - 1
      kt = 1
      go to 500
  490 continue
      idalp = ialp
      in = 0
  500 is = kt
      fidal = idalp
      dalpha = fidal + fnf
      arg = x - pidt*dalpha - pdf
      sa = sin(arg)
      sb = cos(arg)
      coef = rttp/rtx
      etx = 8.0e0*x
  510 continue
      dtm = fidal + fidal
      dtm = dtm*dtm
      tm = 0.0e0
      if (fidal.eq.0.0e0 .and. abs(fnf).lt.tol) go to 520
      tm = 4.0e0*fnf*(fidal+fidal+fnf)
  520 continue
      trx = dtm - 1.0e0
      t2 = (trx+tm)/etx
      s2 = t2
      relb = tol*abs(t2)
      t1 = etx
      s1 = 1.0e0
      fn = 1.0e0
      ak = 8.0e0
      do 530 k=1,13
        t1 = t1 + etx
        fn = fn + ak
        trx = dtm - fn
        ap = trx + tm
        t2 = -t2*ap/t1
        s1 = s1 + t2
        t1 = t1 + etx
        ak = ak + 8.0e0
        fn = fn + ak
        trx = dtm - fn
        ap = trx + tm
        t2 = t2*ap/t1
        s2 = s2 + t2
        if (abs(t2).le.relb) go to 540
        ak = ak + 8.0e0
  530 continue
  540 temp(is) = coef*(s1*sb-s2*sa)
      if(is.eq.2) go to 560
      fidal = fidal + 1.0e0
      dalpha = fidal + fnf
      is = 2
      tb = sa
      sa = -sb
      sb = tb
      go to 510
c
c     forward recursion section
c
  560 if (kt.eq.2) go to 470
      s1 = temp(1)
      s2 = temp(2)
      tx = 2.0e0/x
      tm = dalpha*tx
      if (in.eq.0) go to 580
c
c     forward recur to index alpha
c
      do 570 i=1,in
        s = s2
        s2 = tm*s2 - s1
        tm = tm + tx
        s1 = s
  570 continue
      if (nn.eq.1) go to 600
      s = s2
      s2 = tm*s2 - s1
      tm = tm + tx
      s1 = s
  580 continue
c
c     forward recur from index alpha to alpha+n-1
c
      y(1) = s1
      y(2) = s2
      if (nn.eq.2) return
      do 590 i=3,nn
        y(i) = tm*y(i-1) - y(i-2)
        tm = tm + tx
  590 continue
      return
  600 y(1) = s2
      return
c
c     backward recursion with normalization by
c     asymptotic expansion for nu to infinity or power series.
c
  610 continue
c     computation of last order for series normalization
      akm = max(3.0e0-fn,0.0e0)
      km = int(akm)
      tfn = fn + km
      ta = (gln+tfn-0.9189385332e0-0.0833333333e0/tfn)/(tfn+0.5e0)
      ta = xo2l - ta
      tb = -(1.0e0-1.5e0/tfn)/tfn
      akm = tolln/(-ta+sqrt(ta*ta-tolln*tb)) + 1.5e0
      in = km + int(akm)
      go to 660
  620 continue
c     computation of last order for asymptotic expansion normalization
      gln = wk(3) + wk(2)
      if (wk(6).gt.30.0e0) go to 640
      rden = (pp(4)*wk(6)+pp(3))*wk(6) + 1.0e0
      rzden = pp(1) + pp(2)*wk(6)
      ta = rzden/rden
      if (wk(1).lt.0.10e0) go to 630
      tb = gln/wk(5)
      go to 650
  630 tb=(1.259921049e0+(0.1679894730e0+0.0887944358e0*wk(1))*wk(1))
     1 /wk(7)
      go to 650
  640 continue
      ta = 0.5e0*tolln/wk(4)
      ta=((0.0493827160e0*ta-0.1111111111e0)*ta+0.6666666667e0)*ta*wk(6)
      if (wk(1).lt.0.10e0) go to 630
      tb = gln/wk(5)
  650 in = int(ta/tb+1.5e0)
      if (in.gt.inlim) go to 310
  660 continue
      dtm = fni + in
      trx = 2.0e0/x
      tm = (dtm+fnf)*trx
      ta = 0.0e0
      tb = tol
      kk = 1
      ak=1.0e0
  670 continue
c
c     backward recur unindexed and scale when magnitudes are close to
c     underflow limits (less than slim=r1mach(1)*1.0e+3/tol)
c
      do 680 i=1,in
        s = tb
        tb = tm*tb - ta
        ta = s
        dtm = dtm - 1.0e0
        tm = (dtm+fnf)*trx
  680 continue
c     normalization
      if (kk.ne.1) go to 690
      s=temp(3)
      sa=ta/tb
      ta=s
      tb=s
      if(abs(s).gt.slim) go to 685
      ta=ta*rtol
      tb=tb*rtol
      ak=tol
  685 continue
      ta=ta*sa
      kk = 2
      in = ns
      if (ns.ne.0) go to 670
  690 y(nn) = tb*ak
      nz = n - nn
      if (nn.eq.1) return
      k = nn - 1
      s=tb
      tb = tm*tb - ta
      ta=s
      y(k)=tb*ak
      if (nn.eq.2) return
      dtm = dtm - 1.0e0
      tm = (dtm+fnf)*trx
      k=nn-2
c
c     backward recur indexed
c
      do 700 i=3,nn
        s=tb
        tb = tm*tb - ta
        ta=s
        y(k)=tb*ak
        dtm = dtm - 1.0e0
        tm = (dtm+fnf)*trx
        k = k - 1
  700 continue
      return
c
c
c
  710 continue
      call xermsg ('slatec', 'besj', 'order, alpha, less than zero.',
     +   2, 1)
      return
  720 continue
      call xermsg ('slatec', 'besj', 'n less than one.', 2, 1)
      return
  730 continue
      call xermsg ('slatec', 'besj', 'x less than zero.', 2, 1)
      return
      end
