*deck dbsynu
      subroutine dbsynu (x, fnu, n, y)
c***begin prologue  dbsynu
c***subsidiary
c***purpose  subsidiary to dbesy
c***library   slatec
c***type      double precision (besynu-s, dbsynu-d)
c***author  amos, d. e., (snla)
c***description
c
c     abstract  **** a double precision routine ****
c         dbsynu computes n member sequences of y bessel functions
c         y/sub(fnu+i-1)/(x), i=1,n for non-negative orders fnu and
c         positive x. equations of the references are implemented on
c         small orders dnu for y/sub(dnu)/(x) and y/sub(dnu+1)/(x).
c         forward recursion with the three term recursion relation
c         generates higher orders fnu+i-1, i=1,...,n.
c
c         to start the recursion fnu is normalized to the interval
c         -0.5.le.dnu.lt.0.5. a special form of the power series is
c         implemented on 0.lt.x.le.x1 while the miller algorithm for the
c         k bessel function in terms of the confluent hypergeometric
c         function u(fnu+0.5,2*fnu+1,i*x) is implemented on x1.lt.x.le.x
c         here i is the complex number sqrt(-1.).
c         for x.gt.x2, the asymptotic expansion for large x is used.
c         when fnu is a half odd integer, a special formula for
c         dnu=-0.5 and dnu+1.0=0.5 is used to start the recursion.
c
c         the maximum number of significant digits obtainable
c         is the smaller of 14 and the number of digits carried in
c         double precision arithmetic.
c
c         dbsynu assumes that a significant digit sinh function is
c         available.
c
c     description of arguments
c
c         input
c           x      - x.gt.0.0d0
c           fnu    - order of initial y function, fnu.ge.0.0d0
c           n      - number of members of the sequence, n.ge.1
c
c         output
c           y      - a vector whose first n components contain values
c                    for the sequence y(i)=y/sub(fnu+i-1), i=1,n.
c
c     error conditions
c         improper input arguments - a fatal error
c         overflow - a fatal error
c
c***see also  dbesy
c***references  n. m. temme, on the numerical evaluation of the ordinary
c                 bessel function of the second kind, journal of
c                 computational physics 21, (1976), pp. 343-350.
c               n. m. temme, on the numerical evaluation of the modified
c                 bessel function of the third kind, journal of
c                 computational physics 19, (1975), pp. 324-337.
c***routines called  d1mach, dgamma, xermsg
c***revision history  (yymmdd)
c   800501  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   900328  added type section.  (wrb)
c   900727  added external statement.  (wrb)
c   910408  updated the author and references sections.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dbsynu
c
      integer i, inu, j, k, kk, n, nn
      double precision a,ak,arg,a1,a2,bk,cb,cbk,cc,cck,ck,coef,cpt,
     1 cp1, cp2, cs, cs1, cs2, cx, dnu, dnu2, etest, etx, f, fc, fhs,
     2 fk, fks, flrx, fmu, fn, fnu, fx, g, g1, g2, hpi, p, pi, pt, q,
     3 rb, rbk, rck, relb, rpt, rp1, rp2, rs, rs1, rs2, rthpi, rx, s,
     4 sa, sb, smu, ss, st, s1, s2, tb, tm, tol, t1, t2, x, x1, x2, y
      dimension a(120), rb(120), cb(120), y(*), cc(8)
      double precision dgamma, d1mach
      external dgamma
      save x1, x2,pi, rthpi, hpi, cc
      data x1, x2 / 3.0d0, 20.0d0 /
      data pi,rthpi        / 3.14159265358979d+00, 7.97884560802865d-01/
      data hpi             / 1.57079632679490d+00/
      data cc(1), cc(2), cc(3), cc(4), cc(5), cc(6), cc(7), cc(8)
     1                     / 5.77215664901533d-01,-4.20026350340952d-02,
     2-4.21977345555443d-02, 7.21894324666300d-03,-2.15241674114900d-04,
     3-2.01348547807000d-05, 1.13302723200000d-06, 6.11609500000000d-09/
c***first executable statement  dbsynu
      ak = d1mach(3)
      tol = max(ak,1.0d-15)
      if (x.le.0.0d0) go to 270
      if (fnu.lt.0.0d0) go to 280
      if (n.lt.1) go to 290
      rx = 2.0d0/x
      inu = int(fnu+0.5d0)
      dnu = fnu - inu
      if (abs(dnu).eq.0.5d0) go to 260
      dnu2 = 0.0d0
      if (abs(dnu).lt.tol) go to 10
      dnu2 = dnu*dnu
   10 continue
      if (x.gt.x1) go to 120
c
c     series for x.le.x1
c
      a1 = 1.0d0 - dnu
      a2 = 1.0d0 + dnu
      t1 = 1.0d0/dgamma(a1)
      t2 = 1.0d0/dgamma(a2)
      if (abs(dnu).gt.0.1d0) go to 40
c     series for f0 to resolve indeterminacy for small abs(dnu)
      s = cc(1)
      ak = 1.0d0
      do 20 k=2,8
        ak = ak*dnu2
        tm = cc(k)*ak
        s = s + tm
        if (abs(tm).lt.tol) go to 30
   20 continue
   30 g1 = -(s+s)
      go to 50
   40 continue
      g1 = (t1-t2)/dnu
   50 continue
      g2 = t1 + t2
      smu = 1.0d0
      fc = 1.0d0/pi
      flrx = log(rx)
      fmu = dnu*flrx
      tm = 0.0d0
      if (dnu.eq.0.0d0) go to 60
      tm = sin(dnu*hpi)/dnu
      tm = (dnu+dnu)*tm*tm
      fc = dnu/sin(dnu*pi)
      if (fmu.ne.0.0d0) smu = sinh(fmu)/fmu
   60 continue
      f = fc*(g1*cosh(fmu)+g2*flrx*smu)
      fx = exp(fmu)
      p = fc*t1*fx
      q = fc*t2/fx
      g = f + tm*q
      ak = 1.0d0
      ck = 1.0d0
      bk = 1.0d0
      s1 = g
      s2 = p
      if (inu.gt.0 .or. n.gt.1) go to 90
      if (x.lt.tol) go to 80
      cx = x*x*0.25d0
   70 continue
      f = (ak*f+p+q)/(bk-dnu2)
      p = p/(ak-dnu)
      q = q/(ak+dnu)
      g = f + tm*q
      ck = -ck*cx/ak
      t1 = ck*g
      s1 = s1 + t1
      bk = bk + ak + ak + 1.0d0
      ak = ak + 1.0d0
      s = abs(t1)/(1.0d0+abs(s1))
      if (s.gt.tol) go to 70
   80 continue
      y(1) = -s1
      return
   90 continue
      if (x.lt.tol) go to 110
      cx = x*x*0.25d0
  100 continue
      f = (ak*f+p+q)/(bk-dnu2)
      p = p/(ak-dnu)
      q = q/(ak+dnu)
      g = f + tm*q
      ck = -ck*cx/ak
      t1 = ck*g
      s1 = s1 + t1
      t2 = ck*(p-ak*g)
      s2 = s2 + t2
      bk = bk + ak + ak + 1.0d0
      ak = ak + 1.0d0
      s = abs(t1)/(1.0d0+abs(s1)) + abs(t2)/(1.0d0+abs(s2))
      if (s.gt.tol) go to 100
  110 continue
      s2 = -s2*rx
      s1 = -s1
      go to 160
  120 continue
      coef = rthpi/sqrt(x)
      if (x.gt.x2) go to 210
c
c     miller algorithm for x1.lt.x.le.x2
c
      etest = cos(pi*dnu)/(pi*x*tol)
      fks = 1.0d0
      fhs = 0.25d0
      fk = 0.0d0
      rck = 2.0d0
      cck = x + x
      rp1 = 0.0d0
      cp1 = 0.0d0
      rp2 = 1.0d0
      cp2 = 0.0d0
      k = 0
  130 continue
      k = k + 1
      fk = fk + 1.0d0
      ak = (fhs-dnu2)/(fks+fk)
      pt = fk + 1.0d0
      rbk = rck/pt
      cbk = cck/pt
      rpt = rp2
      cpt = cp2
      rp2 = rbk*rpt - cbk*cpt - ak*rp1
      cp2 = cbk*rpt + rbk*cpt - ak*cp1
      rp1 = rpt
      cp1 = cpt
      rb(k) = rbk
      cb(k) = cbk
      a(k) = ak
      rck = rck + 2.0d0
      fks = fks + fk + fk + 1.0d0
      fhs = fhs + fk + fk
      pt = max(abs(rp1),abs(cp1))
      fc = (rp1/pt)**2 + (cp1/pt)**2
      pt = pt*sqrt(fc)*fk
      if (etest.gt.pt) go to 130
      kk = k
      rs = 1.0d0
      cs = 0.0d0
      rp1 = 0.0d0
      cp1 = 0.0d0
      rp2 = 1.0d0
      cp2 = 0.0d0
      do 140 i=1,k
        rpt = rp2
        cpt = cp2
        rp2 = (rb(kk)*rpt-cb(kk)*cpt-rp1)/a(kk)
        cp2 = (cb(kk)*rpt+rb(kk)*cpt-cp1)/a(kk)
        rp1 = rpt
        cp1 = cpt
        rs = rs + rp2
        cs = cs + cp2
        kk = kk - 1
  140 continue
      pt = max(abs(rs),abs(cs))
      fc = (rs/pt)**2 + (cs/pt)**2
      pt = pt*sqrt(fc)
      rs1 = (rp2*(rs/pt)+cp2*(cs/pt))/pt
      cs1 = (cp2*(rs/pt)-rp2*(cs/pt))/pt
      fc = hpi*(dnu-0.5d0) - x
      p = cos(fc)
      q = sin(fc)
      s1 = (cs1*q-rs1*p)*coef
      if (inu.gt.0 .or. n.gt.1) go to 150
      y(1) = s1
      return
  150 continue
      pt = max(abs(rp2),abs(cp2))
      fc = (rp2/pt)**2 + (cp2/pt)**2
      pt = pt*sqrt(fc)
      rpt = dnu + 0.5d0 - (rp1*(rp2/pt)+cp1*(cp2/pt))/pt
      cpt = x - (cp1*(rp2/pt)-rp1*(cp2/pt))/pt
      cs2 = cs1*cpt - rs1*rpt
      rs2 = rpt*cs1 + rs1*cpt
      s2 = (rs2*q+cs2*p)*coef/x
c
c     forward recursion on the three term recursion relation
c
  160 continue
      ck = (dnu+dnu+2.0d0)/x
      if (n.eq.1) inu = inu - 1
      if (inu.gt.0) go to 170
      if (n.gt.1) go to 190
      s1 = s2
      go to 190
  170 continue
      do 180 i=1,inu
        st = s2
        s2 = ck*s2 - s1
        s1 = st
        ck = ck + rx
  180 continue
      if (n.eq.1) s1 = s2
  190 continue
      y(1) = s1
      if (n.eq.1) return
      y(2) = s2
      if (n.eq.2) return
      do 200 i=3,n
        y(i) = ck*y(i-1) - y(i-2)
        ck = ck + rx
  200 continue
      return
c
c     asymptotic expansion for large x, x.gt.x2
c
  210 continue
      nn = 2
      if (inu.eq.0 .and. n.eq.1) nn = 1
      dnu2 = dnu + dnu
      fmu = 0.0d0
      if (abs(dnu2).lt.tol) go to 220
      fmu = dnu2*dnu2
  220 continue
      arg = x - hpi*(dnu+0.5d0)
      sa = sin(arg)
      sb = cos(arg)
      etx = 8.0d0*x
      do 250 k=1,nn
        s1 = s2
        t2 = (fmu-1.0d0)/etx
        ss = t2
        relb = tol*abs(t2)
        t1 = etx
        s = 1.0d0
        fn = 1.0d0
        ak = 0.0d0
        do 230 j=1,13
          t1 = t1 + etx
          ak = ak + 8.0d0
          fn = fn + ak
          t2 = -t2*(fmu-fn)/t1
          s = s + t2
          t1 = t1 + etx
          ak = ak + 8.0d0
          fn = fn + ak
          t2 = t2*(fmu-fn)/t1
          ss = ss + t2
          if (abs(t2).le.relb) go to 240
  230   continue
  240   s2 = coef*(s*sa+ss*sb)
        fmu = fmu + 8.0d0*dnu + 4.0d0
        tb = sa
        sa = -sb
        sb = tb
  250 continue
      if (nn.gt.1) go to 160
      s1 = s2
      go to 190
c
c     fnu=half odd integer case
c
  260 continue
      coef = rthpi/sqrt(x)
      s1 = coef*sin(x)
      s2 = -coef*cos(x)
      go to 160
c
c
  270 call xermsg ('slatec', 'dbsynu', 'x not greater than zero', 2, 1)
      return
  280 call xermsg ('slatec', 'dbsynu', 'fnu not zero or positive', 2,
     +   1)
      return
  290 call xermsg ('slatec', 'dbsynu', 'n not greater than 0', 2, 1)
      return
      end
