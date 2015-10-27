*deck dpsifn
      subroutine dpsifn (x, n, kode, m, ans, nz, ierr)
c***begin prologue  dpsifn
c***purpose  compute derivatives of the psi function.
c***library   slatec
c***category  c7c
c***type      double precision (psifn-s, dpsifn-d)
c***keywords  derivatives of the gamma function, polygamma function,
c             psi function
c***author  amos, d. e., (snla)
c***description
c
c         the following definitions are used in dpsifn:
c
c      definition 1
c         psi(x) = d/dx (ln(gamma(x)), the first derivative of
c                  the log gamma function.
c      definition 2
c                     k   k
c         psi(k,x) = d /dx (psi(x)), the k-th derivative of psi(x).
c   ___________________________________________________________________
c      dpsifn computes a sequence of scaled derivatives of
c      the psi function; i.e. for fixed x and m it computes
c      the m-member sequence
c
c                    ((-1)**(k+1)/gamma(k+1))*psi(k,x)
c                       for k = n,...,n+m-1
c
c      where psi(k,x) is as defined above.   for kode=1, dpsifn returns
c      the scaled derivatives as described.  kode=2 is operative only
c      when k=0 and in that case dpsifn returns -psi(x) + ln(x).  that
c      is, the logarithmic behavior for large x is removed when kode=2
c      and k=0.  when sums or differences of psi functions are computed
c      the logarithmic terms can be combined analytically and computed
c      separately to help retain significant digits.
c
c         note that call dpsifn(x,0,1,1,ans) results in
c                   ans = -psi(x)
c
c     input      x is double precision
c           x      - argument, x .gt. 0.0d0
c           n      - first member of the sequence, 0 .le. n .le. 100
c                    n=0 gives ans(1) = -psi(x)       for kode=1
c                                       -psi(x)+ln(x) for kode=2
c           kode   - selection parameter
c                    kode=1 returns scaled derivatives of the psi
c                    function.
c                    kode=2 returns scaled derivatives of the psi
c                    function except when n=0. in this case,
c                    ans(1) = -psi(x) + ln(x) is returned.
c           m      - number of members of the sequence, m.ge.1
c
c    output     ans is double precision
c           ans    - a vector of length at least m whose first m
c                    components contain the sequence of derivatives
c                    scaled according to kode.
c           nz     - underflow flag
c                    nz.eq.0, a normal return
c                    nz.ne.0, underflow, last nz components of ans are
c                             set to zero, ans(m-k+1)=0.0, k=1,...,nz
c           ierr   - error flag
c                    ierr=0, a normal return, computation completed
c                    ierr=1, input error,     no computation
c                    ierr=2, overflow,        x too small or n+m-1 too
c                            large or both
c                    ierr=3, error,           n too large. dimensioned
c                            array trmr(nmax) is not large enough for n
c
c         the nominal computational accuracy is the maximum of unit
c         roundoff (=d1mach(4)) and 1.0d-18 since critical constants
c         are given to only 18 digits.
c
c         psifn is the single precision version of dpsifn.
c
c *long description:
c
c         the basic method of evaluation is the asymptotic expansion
c         for large x.ge.xmin followed by backward recursion on a two
c         term recursion relation
c
c                  w(x+1) + x**(-n-1) = w(x).
c
c         this is supplemented by a series
c
c                  sum( (x+k)**(-n-1) , k=0,1,2,... )
c
c         which converges rapidly for large n. both xmin and the
c         number of terms of the series are calculated from the unit
c         roundoff of the machine environment.
c
c***references  handbook of mathematical functions, national bureau
c                 of standards applied mathematics series 55, edited
c                 by m. abramowitz and i. a. stegun, equations 6.3.5,
c                 6.3.18, 6.4.6, 6.4.9 and 6.4.10, pp.258-260, 1964.
c               d. e. amos, a portable fortran subroutine for
c                 derivatives of the psi function, algorithm 610, acm
c                 transactions on mathematical software 9, 4 (1983),
c                 pp. 494-502.
c***routines called  d1mach, i1mach
c***revision history  (yymmdd)
c   820601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   891006  cosmetic changes to prologue.  (wrb)
c   891006  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dpsifn
      integer i, ierr, j, k, kode, m, mm, mx, n, nmax, nn, np, nx, nz,
     *  fn
      integer i1mach
      double precision ans, arg, b, den, elim, eps, fln,
     * fx, rln, rxsq, r1m4, r1m5, s, slope, t, ta, tk, tol, tols, trm,
     * trmr, tss, tst, tt, t1, t2, wdtol, x, xdmln, xdmy, xinc, xln,
     * xm, xmin, xq, yint
      double precision d1mach
      dimension b(22), trm(22), trmr(100), ans(*)
      save nmax, b
      data nmax /100/
c-----------------------------------------------------------------------
c             bernoulli numbers
c-----------------------------------------------------------------------
      data b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8), b(9), b(10),
     * b(11), b(12), b(13), b(14), b(15), b(16), b(17), b(18), b(19),
     * b(20), b(21), b(22) /1.00000000000000000d+00,
     * -5.00000000000000000d-01,1.66666666666666667d-01,
     * -3.33333333333333333d-02,2.38095238095238095d-02,
     * -3.33333333333333333d-02,7.57575757575757576d-02,
     * -2.53113553113553114d-01,1.16666666666666667d+00,
     * -7.09215686274509804d+00,5.49711779448621554d+01,
     * -5.29124242424242424d+02,6.19212318840579710d+03,
     * -8.65802531135531136d+04,1.42551716666666667d+06,
     * -2.72982310678160920d+07,6.01580873900642368d+08,
     * -1.51163157670921569d+10,4.29614643061166667d+11,
     * -1.37116552050883328d+13,4.88332318973593167d+14,
     * -1.92965793419400681d+16/
c
c***first executable statement  dpsifn
      ierr = 0
      nz=0
      if (x.le.0.0d0) ierr=1
      if (n.lt.0) ierr=1
      if (kode.lt.1 .or. kode.gt.2) ierr=1
      if (m.lt.1) ierr=1
      if (ierr.ne.0) return
      mm=m
      nx = min(-i1mach(15),i1mach(16))
      r1m5 = d1mach(5)
      r1m4 = d1mach(4)*0.5d0
      wdtol = max(r1m4,0.5d-18)
c-----------------------------------------------------------------------
c     elim = approximate exponential over and underflow limit
c-----------------------------------------------------------------------
      elim = 2.302d0*(nx*r1m5-3.0d0)
      xln = log(x)
   41 continue
      nn = n + mm - 1
      fn = nn
      t = (fn+1)*xln
c-----------------------------------------------------------------------
c     overflow and underflow test for small and large x
c-----------------------------------------------------------------------
      if (abs(t).gt.elim) go to 290
      if (x.lt.wdtol) go to 260
c-----------------------------------------------------------------------
c     compute xmin and the number of terms of the series, fln+1
c-----------------------------------------------------------------------
      rln = r1m5*i1mach(14)
      rln = min(rln,18.06d0)
      fln = max(rln,3.0d0) - 3.0d0
      yint = 3.50d0 + 0.40d0*fln
      slope = 0.21d0 + fln*(0.0006038d0*fln+0.008677d0)
      xm = yint + slope*fn
      mx = int(xm) + 1
      xmin = mx
      if (n.eq.0) go to 50
      xm = -2.302d0*rln - min(0.0d0,xln)
      arg = xm/n
      arg = min(0.0d0,arg)
      eps = exp(arg)
      xm = 1.0d0 - eps
      if (abs(arg).lt.1.0d-3) xm = -arg
      fln = x*xm/eps
      xm = xmin - x
      if (xm.gt.7.0d0 .and. fln.lt.15.0d0) go to 200
   50 continue
      xdmy = x
      xdmln = xln
      xinc = 0.0d0
      if (x.ge.xmin) go to 60
      nx = int(x)
      xinc = xmin - nx
      xdmy = x + xinc
      xdmln = log(xdmy)
   60 continue
c-----------------------------------------------------------------------
c     generate w(n+mm-1,x) by the asymptotic expansion
c-----------------------------------------------------------------------
      t = fn*xdmln
      t1 = xdmln + xdmln
      t2 = t + xdmln
      tk = max(abs(t),abs(t1),abs(t2))
      if (tk.gt.elim) go to 380
      tss = exp(-t)
      tt = 0.5d0/xdmy
      t1 = tt
      tst = wdtol*tt
      if (nn.ne.0) t1 = tt + 1.0d0/fn
      rxsq = 1.0d0/(xdmy*xdmy)
      ta = 0.5d0*rxsq
      t = (fn+1)*ta
      s = t*b(3)
      if (abs(s).lt.tst) go to 80
      tk = 2.0d0
      do 70 k=4,22
        t = t*((tk+fn+1)/(tk+1.0d0))*((tk+fn)/(tk+2.0d0))*rxsq
        trm(k) = t*b(k)
        if (abs(trm(k)).lt.tst) go to 80
        s = s + trm(k)
        tk = tk + 2.0d0
   70 continue
   80 continue
      s = (s+t1)*tss
      if (xinc.eq.0.0d0) go to 100
c-----------------------------------------------------------------------
c     backward recur from xdmy to x
c-----------------------------------------------------------------------
      nx = int(xinc)
      np = nn + 1
      if (nx.gt.nmax) go to 390
      if (nn.eq.0) go to 160
      xm = xinc - 1.0d0
      fx = x + xm
c-----------------------------------------------------------------------
c     this loop should not be changed. fx is accurate when x is small
c-----------------------------------------------------------------------
      do 90 i=1,nx
        trmr(i) = fx**(-np)
        s = s + trmr(i)
        xm = xm - 1.0d0
        fx = x + xm
   90 continue
  100 continue
      ans(mm) = s
      if (fn.eq.0) go to 180
c-----------------------------------------------------------------------
c     generate lower derivatives, j.lt.n+mm-1
c-----------------------------------------------------------------------
      if (mm.eq.1) return
      do 150 j=2,mm
        fn = fn - 1
        tss = tss*xdmy
        t1 = tt
        if (fn.ne.0) t1 = tt + 1.0d0/fn
        t = (fn+1)*ta
        s = t*b(3)
        if (abs(s).lt.tst) go to 120
        tk = 4 + fn
        do 110 k=4,22
          trm(k) = trm(k)*(fn+1)/tk
          if (abs(trm(k)).lt.tst) go to 120
          s = s + trm(k)
          tk = tk + 2.0d0
  110   continue
  120   continue
        s = (s+t1)*tss
        if (xinc.eq.0.0d0) go to 140
        if (fn.eq.0) go to 160
        xm = xinc - 1.0d0
        fx = x + xm
        do 130 i=1,nx
          trmr(i) = trmr(i)*fx
          s = s + trmr(i)
          xm = xm - 1.0d0
          fx = x + xm
  130   continue
  140   continue
        mx = mm - j + 1
        ans(mx) = s
        if (fn.eq.0) go to 180
  150 continue
      return
c-----------------------------------------------------------------------
c     recursion for n = 0
c-----------------------------------------------------------------------
  160 continue
      do 170 i=1,nx
        s = s + 1.0d0/(x+nx-i)
  170 continue
  180 continue
      if (kode.eq.2) go to 190
      ans(1) = s - xdmln
      return
  190 continue
      if (xdmy.eq.x) return
      xq = xdmy/x
      ans(1) = s - log(xq)
      return
c-----------------------------------------------------------------------
c     compute by series (x+k)**(-(n+1)) , k=0,1,2,...
c-----------------------------------------------------------------------
  200 continue
      nn = int(fln) + 1
      np = n + 1
      t1 = (n+1)*xln
      t = exp(-t1)
      s = t
      den = x
      do 210 i=1,nn
        den = den + 1.0d0
        trm(i) = den**(-np)
        s = s + trm(i)
  210 continue
      ans(1) = s
      if (n.ne.0) go to 220
      if (kode.eq.2) ans(1) = s + xln
  220 continue
      if (mm.eq.1) return
c-----------------------------------------------------------------------
c     generate higher derivatives, j.gt.n
c-----------------------------------------------------------------------
      tol = wdtol/5.0d0
      do 250 j=2,mm
        t = t/x
        s = t
        tols = t*tol
        den = x
        do 230 i=1,nn
          den = den + 1.0d0
          trm(i) = trm(i)/den
          s = s + trm(i)
          if (trm(i).lt.tols) go to 240
  230   continue
  240   continue
        ans(j) = s
  250 continue
      return
c-----------------------------------------------------------------------
c     small x.lt.unit round off
c-----------------------------------------------------------------------
  260 continue
      ans(1) = x**(-n-1)
      if (mm.eq.1) go to 280
      k = 1
      do 270 i=2,mm
        ans(k+1) = ans(k)/x
        k = k + 1
  270 continue
  280 continue
      if (n.ne.0) return
      if (kode.eq.2) ans(1) = ans(1) + xln
      return
  290 continue
      if (t.gt.0.0d0) go to 380
      nz=0
      ierr=2
      return
  380 continue
      nz=nz+1
      ans(mm)=0.0d0
      mm=mm-1
      if (mm.eq.0) return
      go to 41
  390 continue
      nz=0
      ierr=3
      return
      end
