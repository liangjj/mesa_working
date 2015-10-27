*deck dhkseq
      subroutine dhkseq (x, m, h, ierr)
c***begin prologue  dhkseq
c***subsidiary
c***purpose  subsidiary to dbskin
c***library   slatec
c***type      double precision (hkseq-s, dhkseq-d)
c***author  amos, d. e., (snla)
c***description
c
c   dhkseq is an adaptation of subroutine dpsifn described in the
c   reference below.  dhkseq generates the sequence
c   h(k,x) = (-x)**(k+1)*(psi(k,x) psi(k,x+0.5))/gamma(k+1), for
c            k=0,...,m.
c
c***see also  dbskin
c***references  d. e. amos, a portable fortran subroutine for
c                 derivatives of the psi function, algorithm 610, acm
c                 transactions on mathematical software 9, 4 (1983),
c                 pp. 494-502.
c***routines called  d1mach, i1mach
c***revision history  (yymmdd)
c   820601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c   920528  description and references sections revised.  (wrb)
c***end prologue  dhkseq
      integer i, ierr, j, k, m, mx, nx
      integer i1mach
      double precision b, fk, fln, fn, fnp, h, hrx, rln, rxsq, r1m5, s,
     * slope, t, tk, trm, trmh, trmr, tst, u, v, wdtol, x, xdmy, xh,
     * xinc, xm, xmin, yint
      double precision d1mach
      dimension b(22), trm(22), trmr(25), trmh(25), u(25), v(25), h(*)
      save b
c-----------------------------------------------------------------------
c             scaled bernoulli numbers 2.0*b(2k)*(1-2**(-2k))
c-----------------------------------------------------------------------
      data b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8), b(9), b(10),
     * b(11), b(12), b(13), b(14), b(15), b(16), b(17), b(18), b(19),
     * b(20), b(21), b(22) /1.00000000000000000d+00,
     * -5.00000000000000000d-01,2.50000000000000000d-01,
     * -6.25000000000000000d-02,4.68750000000000000d-02,
     * -6.64062500000000000d-02,1.51367187500000000d-01,
     * -5.06103515625000000d-01,2.33319091796875000d+00,
     * -1.41840972900390625d+01,1.09941936492919922d+02,
     * -1.05824747562408447d+03,1.23842434241771698d+04,
     * -1.73160495905935764d+05,2.85103429084961116d+06,
     * -5.45964619322445132d+07,1.20316174668075304d+09,
     * -3.02326315271452307d+10,8.59229286072319606d+11,
     * -2.74233104097776039d+13,9.76664637943633248d+14,
     * -3.85931586838450360d+16/
c
c***first executable statement  dhkseq
      ierr=0
      wdtol = max(d1mach(4),1.0d-18)
      fn = m - 1
      fnp = fn + 1.0d0
c-----------------------------------------------------------------------
c     compute xmin
c-----------------------------------------------------------------------
      r1m5 = d1mach(5)
      rln = r1m5*i1mach(14)
      rln = min(rln,18.06d0)
      fln = max(rln,3.0d0) - 3.0d0
      yint = 3.50d0 + 0.40d0*fln
      slope = 0.21d0 + fln*(0.0006038d0*fln+0.008677d0)
      xm = yint + slope*fn
      mx = int(xm) + 1
      xmin = mx
c-----------------------------------------------------------------------
c     generate h(m-1,xdmy)*xdmy**(m) by the asymptotic expansion
c-----------------------------------------------------------------------
      xdmy = x
      xinc = 0.0d0
      if (x.ge.xmin) go to 10
      nx = int(x)
      xinc = xmin - nx
      xdmy = x + xinc
   10 continue
      rxsq = 1.0d0/(xdmy*xdmy)
      hrx = 0.5d0/xdmy
      tst = 0.5d0*wdtol
      t = fnp*hrx
c-----------------------------------------------------------------------
c     initialize coefficient array
c-----------------------------------------------------------------------
      s = t*b(3)
      if (abs(s).lt.tst) go to 30
      tk = 2.0d0
      do 20 k=4,22
        t = t*((tk+fn+1.0d0)/(tk+1.0d0))*((tk+fn)/(tk+2.0d0))*rxsq
        trm(k) = t*b(k)
        if (abs(trm(k)).lt.tst) go to 30
        s = s + trm(k)
        tk = tk + 2.0d0
   20 continue
      go to 110
   30 continue
      h(m) = s + 0.5d0
      if (m.eq.1) go to 70
c-----------------------------------------------------------------------
c     generate lower derivatives, i.lt.m-1
c-----------------------------------------------------------------------
      do 60 i=2,m
        fnp = fn
        fn = fn - 1.0d0
        s = fnp*hrx*b(3)
        if (abs(s).lt.tst) go to 50
        fk = fnp + 3.0d0
        do 40 k=4,22
          trm(k) = trm(k)*fnp/fk
          if (abs(trm(k)).lt.tst) go to 50
          s = s + trm(k)
          fk = fk + 2.0d0
   40   continue
        go to 110
   50   continue
        mx = m - i + 1
        h(mx) = s + 0.5d0
   60 continue
   70 continue
      if (xinc.eq.0.0d0) return
c-----------------------------------------------------------------------
c     recur backward from xdmy to x
c-----------------------------------------------------------------------
      xh = x + 0.5d0
      s = 0.0d0
      nx = int(xinc)
      do 80 i=1,nx
        trmr(i) = x/(x+nx-i)
        u(i) = trmr(i)
        trmh(i) = x/(xh+nx-i)
        v(i) = trmh(i)
        s = s + u(i) - v(i)
   80 continue
      mx = nx + 1
      trmr(mx) = x/xdmy
      u(mx) = trmr(mx)
      h(1) = h(1)*trmr(mx) + s
      if (m.eq.1) return
      do 100 j=2,m
        s = 0.0d0
        do 90 i=1,nx
          trmr(i) = trmr(i)*u(i)
          trmh(i) = trmh(i)*v(i)
          s = s + trmr(i) - trmh(i)
   90   continue
        trmr(mx) = trmr(mx)*u(mx)
        h(j) = h(j)*trmr(mx) + s
  100 continue
      return
  110 continue
      ierr=2
      return
      end
