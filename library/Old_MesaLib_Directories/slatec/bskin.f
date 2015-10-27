*deck bskin
      subroutine bskin (x, n, kode, m, y, nz, ierr)
c***begin prologue  bskin
c***purpose  compute repeated integrals of the k-zero bessel function.
c***library   slatec
c***category  c10f
c***type      single precision (bskin-s, dbskin-d)
c***keywords  bickley functions, exponential integral,
c             integrals of bessel functions, k-zero bessel function
c***author  amos, d. e., (snla)
c***description
c
c         the following definitions are used in bskin:
c
c   definition 1
c         ki(0,x) = k-zero bessel function.
c
c   definition 2
c         ki(n,x) = bickley function
c                 =  integral from x to infinity of ki(n-1,t)dt
c                     for x .ge. 0 and n = 1,2,...
c   ____________________________________________________________________
c      bskin computes sequences of bickley functions (repeated integrals
c      of the k0 bessel function); i.e. for fixed x and n and k=1,...,
c      bskin computes the m-member sequence
c
c                     y(k) =        ki(n+k-1,x) for kode=1
c      or
c                     y(k) = exp(x)*ki(n+k-1,x) for kode=2,
c
c      for n.ge.0 and x.ge.0 (n and x cannot be zero simultaneously).
c
c      input
c        x      - argument, x .ge. 0.0e0
c        n      - order of first member of the sequence n .ge. 0
c        kode   - selection parameter
c                 kode = 1 returns y(k)=       ki(n+k-1,x), k=1,m
c                      = 2 returns y(k)=exp(x)*ki(n+k-1,x), k=1,m
c        m      - number of members in the sequence, m.ge.1
c
c      output
c        y      - a vector of dimension at least m containing the
c                 sequence selected by kode.
c        nz     - underflow flag
c                 nz = 0 means computation completed
c                    = m means an exponential underflow occurred on
c                        kode=1.  y(k)=0.0e0, k=1,...,m is returned
c        ierr   - error flag
c                 ierr = 0, normal return, computation completed.
c                      = 1, input error,   no computation.
c                      = 2, error,         no computation.  the
c                           termination condition was not met.
c
c      the nominal computational accuracy is the maximum of unit
c      roundoff (=r1mach(4)) and 1.0e-18 since critical constants
c      are given to only 18 digits.
c
c      dbskin is the double precision version of bskin.
c
c *long description:
c
c         numerical recurrence on
c
c      (l-1)*ki(l,x) = x(ki(l-3,x) - ki(l-1,x)) + (l-2)*ki(l-2,x)
c
c         is stable where recurrence is carried forward or backward
c         away from int(x+0.5).  the power series for indices 0,1 and 2
c         on 0.le.x.le. 2 starts a stable recurrence for indices
c         greater than 2.  if n is sufficiently large (n.gt.nlim), the
c         uniform asymptotic expansion for n to infinity is more
c         economical.  on x.gt.2 the recursion is started by evaluating
c         the uniform expansion for the three members whose indices are
c         closest to int(x+0.5) within the set n,...,n+m-1.  forward
c         recurrence, backward recurrence or both, complete the
c         sequence depending on the relation of int(x+0.5) to the
c         indices n,...,n+m-1.
c
c***references  d. e. amos, uniform asymptotic expansions for
c                 exponential integrals e(n,x) and bickley functions
c                 ki(n,x), acm transactions on mathematical software,
c                 1983.
c               d. e. amos, a portable fortran subroutine for the
c                 bickley functions ki(n,x), algorithm 609, acm
c                 transactions on mathematical software, 1983.
c***routines called  bkias, bkisr, exint, gamrn, i1mach, r1mach
c***revision history  (yymmdd)
c   820601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891009  removed unreferenced statement label.  (wrb)
c   891009  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  bskin
      integer i, icase, ierr, il, i1m, k, kk, kode, ktrms, m,
     * m3, n, ne, nflg, nl, nlim, nn, np, ns, nt, nz
      integer i1mach
      real a, enlim, exi, fn, gr, h, hn, hrtpi, ss, tol, t1, t2, w, x,
     * xlim, xnlim, xp, y, ys, yss
      real gamrn, r1mach
      dimension exi(102), a(50), ys(3), yss(3), h(31), y(*)
      save a, hrtpi
c-----------------------------------------------------------------------
c             coefficients in series of exponential integrals
c-----------------------------------------------------------------------
      data a(1), a(2), a(3), a(4), a(5), a(6), a(7), a(8), a(9), a(10),
     * a(11), a(12), a(13), a(14), a(15), a(16), a(17), a(18), a(19),
     * a(20), a(21), a(22), a(23), a(24) /1.00000000000000000e+00,
     * 5.00000000000000000e-01,3.75000000000000000e-01,
     * 3.12500000000000000e-01,2.73437500000000000e-01,
     * 2.46093750000000000e-01,2.25585937500000000e-01,
     * 2.09472656250000000e-01,1.96380615234375000e-01,
     * 1.85470581054687500e-01,1.76197052001953125e-01,
     * 1.68188095092773438e-01,1.61180257797241211e-01,
     * 1.54981017112731934e-01,1.49445980787277222e-01,
     * 1.44464448094367981e-01,1.39949934091418982e-01,
     * 1.35833759559318423e-01,1.32060599571559578e-01,
     * 1.28585320635465905e-01,1.25370687619579257e-01,
     * 1.22385671247684513e-01,1.19604178719328047e-01,
     * 1.17004087877603524e-01/
      data a(25), a(26), a(27), a(28), a(29), a(30), a(31), a(32),
     * a(33), a(34), a(35), a(36), a(37), a(38), a(39), a(40), a(41),
     * a(42), a(43), a(44), a(45), a(46), a(47), a(48)
     * /1.14566502713486784e-01,1.12275172659217048e-01,
     * 1.10116034723462874e-01,1.08076848895250599e-01,
     * 1.06146905164978267e-01,1.04316786110409676e-01,
     * 1.02578173008569515e-01,1.00923686347140974e-01,
     * 9.93467537479668965e-02,9.78414999033007314e-02,
     * 9.64026543164874854e-02,9.50254735405376642e-02,
     * 9.37056752969190855e-02,9.24393823875012600e-02,
     * 9.12230747245078224e-02,9.00535481254756708e-02,
     * 8.89278787739072249e-02,8.78433924473961612e-02,
     * 8.67976377754033498e-02,8.57883629175498224e-02,
     * 8.48134951571231199e-02,8.38711229887106408e-02,
     * 8.29594803475290034e-02,8.20769326842574183e-02/
      data a(49), a(50) /8.12219646354630702e-02,8.03931690779583449e-02
     * /
c-----------------------------------------------------------------------
c             sqrt(pi)/2
c-----------------------------------------------------------------------
      data hrtpi /8.86226925452758014e-01/
c
c***first executable statement  bskin
      ierr = 0
      nz=0
      if (x.lt.0.0e0) ierr=1
      if (n.lt.0) ierr=1
      if (kode.lt.1 .or. kode.gt.2) ierr=1
      if (m.lt.1) ierr=1
      if (x.eq.0.0e0 .and. n.eq.0) ierr=1
      if (ierr.ne.0) return
      if (x.eq.0.0e0) go to 300
      i1m = -i1mach(12)
      t1 = 2.3026e0*r1mach(5)*i1m
      xlim = t1 - 3.228086e0
      t2 = t1 + n + m - 1
      if (t2.gt.1000.0e0) xlim = t1 - 0.5e0*(log(t2)-0.451583e0)
      if (x.gt.xlim .and. kode.eq.1) go to 320
      tol = max(r1mach(4),1.0e-18)
      i1m = i1mach(11)
c-----------------------------------------------------------------------
c     ln(nlim) = 0.125*ln(eps),   nlim = 2*ktrms+n
c-----------------------------------------------------------------------
      xnlim = 0.287823e0*(i1m-1)*r1mach(5)
      enlim = exp(xnlim)
      nlim = int(enlim) + 2
      nlim = min(100,nlim)
      nlim = max(20,nlim)
      m3 = min(m,3)
      nl = n + m - 1
      if (x.gt.2.0e0) go to 130
      if (n.gt.nlim) go to 280
c-----------------------------------------------------------------------
c     computation by series for 0.le.x.le.2
c-----------------------------------------------------------------------
      nflg = 0
      nn = n
      if (nl.le.2) go to 60
      m3 = 3
      nn = 0
      nflg = 1
   60 continue
      xp = 1.0e0
      if (kode.eq.2) xp = exp(x)
      do 80 i=1,m3
        call bkisr(x, nn, w, ierr)
      if(ierr.ne.0) return
        w = w*xp
        if (nn.lt.n) go to 70
        kk = nn - n + 1
        y(kk) = w
   70   continue
        ys(i) = w
        nn = nn + 1
   80 continue
      if (nflg.eq.0) return
      ns = nn
      xp = 1.0e0
   90 continue
c-----------------------------------------------------------------------
c     forward recursion scaled by exp(x) on icase=0,1,2
c-----------------------------------------------------------------------
      fn = ns - 1
      il = nl - ns + 1
      if (il.le.0) return
      do 110 i=1,il
        t1 = ys(2)
        t2 = ys(3)
        ys(3) = (x*(ys(1)-ys(3))+(fn-1.0e0)*ys(2))/fn
        ys(2) = t2
        ys(1) = t1
        fn = fn + 1.0e0
        if (ns.lt.n) go to 100
        kk = ns - n + 1
        y(kk) = ys(3)*xp
  100   continue
        ns = ns + 1
  110 continue
      return
c-----------------------------------------------------------------------
c     computation by asymptotic expansion for x.gt.2
c-----------------------------------------------------------------------
  130 continue
      w = x + 0.5e0
      nt = int(w)
      if (nl.gt.nt) go to 270
c-----------------------------------------------------------------------
c     case nl.le.nt, icase=0
c-----------------------------------------------------------------------
      icase = 0
      nn = nl
      nflg = min(m-m3,1)
  140 continue
      kk = (nlim-nn)/2
      ktrms = max(0,kk)
      ns = nn + 1
      np = nn - m3 + 1
      xp = 1.0e0
      if (kode.eq.1) xp = exp(-x)
      do 150 i=1,m3
        kk = i
        call bkias(x, np, ktrms, a, w, kk, ne, gr, h, ierr)
      if(ierr.ne.0) return
        ys(i) = w
        np = np + 1
  150 continue
c-----------------------------------------------------------------------
c     sum series of exponential integrals backward
c-----------------------------------------------------------------------
      if (ktrms.eq.0) go to 160
      ne = ktrms + ktrms + 1
      np = nn - m3 + 2
      call exint(x, np, 2, ne, tol, exi, nz, ierr)
      if(nz.ne.0) go to 320
      if(ierr.eq.2) return
  160 continue
      do 190 i=1,m3
        ss = 0.0e0
        if (ktrms.eq.0) go to 180
        kk = i + ktrms + ktrms - 2
        il = ktrms
        do 170 k=1,ktrms
          ss = ss + a(il)*exi(kk)
          kk = kk - 2
          il = il - 1
  170   continue
  180   continue
        ys(i) = ys(i) + ss
  190 continue
      if (icase.eq.1) go to 200
      if (nflg.ne.0) go to 220
  200 continue
      do 210 i=1,m3
        y(i) = ys(i)*xp
  210 continue
      if (icase.eq.1 .and. nflg.eq.1) go to 90
      return
  220 continue
c-----------------------------------------------------------------------
c     backward recursion scaled by exp(x) icase=0,2
c-----------------------------------------------------------------------
      kk = nn - n + 1
      k = m3
      do 230 i=1,m3
        y(kk) = ys(k)*xp
        yss(i) = ys(i)
        kk = kk - 1
        k = k - 1
  230 continue
      il = kk
      if (il.le.0) go to 250
      fn = nn - 3
      do 240 i=1,il
        t1 = ys(2)
        t2 = ys(1)
        ys(1) = ys(2) + ((fn+2.0e0)*ys(3)-(fn+1.0e0)*ys(1))/x
        ys(2) = t2
        ys(3) = t1
        y(kk) = ys(1)*xp
        kk = kk - 1
        fn = fn - 1.0e0
  240 continue
  250 continue
      if (icase.ne.2) return
      do 260 i=1,m3
        ys(i) = yss(i)
  260 continue
      go to 90
  270 continue
      if (n.lt.nt) go to 290
c-----------------------------------------------------------------------
c     icase=1, nt.le.n.le.nl with forward recursion
c-----------------------------------------------------------------------
  280 continue
      nn = n + m3 - 1
      nflg = min(m-m3,1)
      icase = 1
      go to 140
c-----------------------------------------------------------------------
c     icase=2, n.lt.nt.lt.nl with both forward and backward recursion
c-----------------------------------------------------------------------
  290 continue
      nn = nt + 1
      nflg = min(m-m3,1)
      icase = 2
      go to 140
c-----------------------------------------------------------------------
c     x=0 case
c-----------------------------------------------------------------------
  300 continue
      fn = n
      hn = 0.5e0*fn
      gr = gamrn(hn)
      y(1) = hrtpi*gr
      if (m.eq.1) return
      y(2) = hrtpi/(hn*gr)
      if (m.eq.2) return
      do 310 k=3,m
        y(k) = fn*y(k-2)/(fn+1.0e0)
        fn = fn + 1.0e0
  310 continue
      return
c-----------------------------------------------------------------------
c     underflow on kode=1, x.gt.xlim
c-----------------------------------------------------------------------
  320 continue
      nz=m
      do 330 i=1,m
        y(i) = 0.0e0
  330 continue
      return
      end
