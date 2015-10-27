*deck dbsknu
      subroutine dbsknu (x, fnu, kode, n, y, nz)
c***begin prologue  dbsknu
c***subsidiary
c***purpose  subsidiary to dbesk
c***library   slatec
c***type      double precision (besknu-s, dbsknu-d)
c***author  amos, d. e., (snla)
c***description
c
c     abstract  **** a double precision routine ****
c         dbsknu computes n member sequences of k bessel functions
c         k/sub(fnu+i-1)/(x), i=1,n for non-negative orders fnu and
c         positive x. equations of the references are implemented on
c         small orders dnu for k/sub(dnu)/(x) and k/sub(dnu+1)/(x).
c         forward recursion with the three term recursion relation
c         generates higher orders fnu+i-1, i=1,...,n. the parameter
c         kode permits k/sub(fnu+i-1)/(x) values or scaled values
c         exp(x)*k/sub(fnu+i-1)/(x), i=1,n to be returned.
c
c         to start the recursion fnu is normalized to the interval
c         -0.5.le.dnu.lt.0.5. a special form of the power series is
c         implemented on 0.lt.x.le.x1 while the miller algorithm for the
c         k bessel function in terms of the confluent hypergeometric
c         function u(fnu+0.5,2*fnu+1,x) is implemented on x1.lt.x.le.x2.
c         for x.gt.x2, the asymptotic expansion for large x is used.
c         when fnu is a half odd integer, a special formula for
c         dnu=-0.5 and dnu+1.0=0.5 is used to start the recursion.
c
c         the maximum number of significant digits obtainable
c         is the smaller of 14 and the number of digits carried in
c         double precision arithmetic.
c
c         dbsknu assumes that a significant digit sinh function is
c         available.
c
c     description of arguments
c
c         input      x,fnu are double precision
c           x      - x.gt.0.0d0
c           fnu    - order of initial k function, fnu.ge.0.0d0
c           n      - number of members of the sequence, n.ge.1
c           kode   - a parameter to indicate the scaling option
c                    kode= 1  returns
c                             y(i)=       k/sub(fnu+i-1)/(x)
c                                  i=1,...,n
c                        = 2  returns
c                             y(i)=exp(x)*k/sub(fnu+i-1)/(x)
c                                  i=1,...,n
c
c         output     y is double precision
c           y      - a vector whose first n components contain values
c                    for the sequence
c                    y(i)=       k/sub(fnu+i-1)/(x), i=1,...,n or
c                    y(i)=exp(x)*k/sub(fnu+i-1)/(x), i=1,...,n
c                    depending on kode
c           nz     - number of components set to zero due to
c                    underflow,
c                    nz= 0   , normal return
c                    nz.ne.0 , first nz components of y set to zero
c                              due to underflow, y(i)=0.0d0,i=1,...,nz
c
c     error conditions
c         improper input arguments - a fatal error
c         overflow - a fatal error
c         underflow with kode=1 - a non-fatal error (nz.ne.0)
c
c***see also  dbesk
c***references  n. m. temme, on the numerical evaluation of the modified
c                 bessel function of the third kind, journal of
c                 computational physics 19, (1975), pp. 324-337.
c***routines called  d1mach, dgamma, i1mach, xermsg
c***revision history  (yymmdd)
c   790201  date written
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
c***end prologue  dbsknu
c
      integer i, iflag, inu, j, k, kk, kode, koded, n, nn, nz
      integer i1mach
      double precision a,ak,a1,a2,b,bk,cc,ck,coef,cx,dk,dnu,dnu2,elim,
     1 etest, ex, f, fc, fhs, fk, fks, flrx, fmu, fnu, g1, g2, p, pi,
     2 pt, p1, p2, q, rthpi, rx, s, smu, sqk, st, s1, s2, tm, tol, t1,
     3 t2, x, x1, x2, y
      dimension a(160), b(160), y(*), cc(8)
      double precision dgamma, d1mach
      external dgamma
      save x1, x2, pi, rthpi, cc
      data x1, x2 / 2.0d0, 17.0d0 /
      data pi,rthpi        / 3.14159265358979d+00, 1.25331413731550d+00/
      data cc(1), cc(2), cc(3), cc(4), cc(5), cc(6), cc(7), cc(8)
     1                     / 5.77215664901533d-01,-4.20026350340952d-02,
     2-4.21977345555443d-02, 7.21894324666300d-03,-2.15241674114900d-04,
     3-2.01348547807000d-05, 1.13302723200000d-06, 6.11609500000000d-09/
c***first executable statement  dbsknu
      kk = -i1mach(15)
      elim = 2.303d0*(kk*d1mach(5)-3.0d0)
      ak = d1mach(3)
      tol = max(ak,1.0d-15)
      if (x.le.0.0d0) go to 350
      if (fnu.lt.0.0d0) go to 360
      if (kode.lt.1 .or. kode.gt.2) go to 370
      if (n.lt.1) go to 380
      nz = 0
      iflag = 0
      koded = kode
      rx = 2.0d0/x
      inu = int(fnu+0.5d0)
      dnu = fnu - inu
      if (abs(dnu).eq.0.5d0) go to 120
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
   30 g1 = -s
      go to 50
   40 continue
      g1 = (t1-t2)/(dnu+dnu)
   50 continue
      g2 = (t1+t2)*0.5d0
      smu = 1.0d0
      fc = 1.0d0
      flrx = log(rx)
      fmu = dnu*flrx
      if (dnu.eq.0.0d0) go to 60
      fc = dnu*pi
      fc = fc/sin(fc)
      if (fmu.ne.0.0d0) smu = sinh(fmu)/fmu
   60 continue
      f = fc*(g1*cosh(fmu)+g2*flrx*smu)
      fc = exp(fmu)
      p = 0.5d0*fc/t2
      q = 0.5d0/(fc*t1)
      ak = 1.0d0
      ck = 1.0d0
      bk = 1.0d0
      s1 = f
      s2 = p
      if (inu.gt.0 .or. n.gt.1) go to 90
      if (x.lt.tol) go to 80
      cx = x*x*0.25d0
   70 continue
      f = (ak*f+p+q)/(bk-dnu2)
      p = p/(ak-dnu)
      q = q/(ak+dnu)
      ck = ck*cx/ak
      t1 = ck*f
      s1 = s1 + t1
      bk = bk + ak + ak + 1.0d0
      ak = ak + 1.0d0
      s = abs(t1)/(1.0d0+abs(s1))
      if (s.gt.tol) go to 70
   80 continue
      y(1) = s1
      if (koded.eq.1) return
      y(1) = s1*exp(x)
      return
   90 continue
      if (x.lt.tol) go to 110
      cx = x*x*0.25d0
  100 continue
      f = (ak*f+p+q)/(bk-dnu2)
      p = p/(ak-dnu)
      q = q/(ak+dnu)
      ck = ck*cx/ak
      t1 = ck*f
      s1 = s1 + t1
      t2 = ck*(p-ak*f)
      s2 = s2 + t2
      bk = bk + ak + ak + 1.0d0
      ak = ak + 1.0d0
      s = abs(t1)/(1.0d0+abs(s1)) + abs(t2)/(1.0d0+abs(s2))
      if (s.gt.tol) go to 100
  110 continue
      s2 = s2*rx
      if (koded.eq.1) go to 170
      f = exp(x)
      s1 = s1*f
      s2 = s2*f
      go to 170
  120 continue
      coef = rthpi/sqrt(x)
      if (koded.eq.2) go to 130
      if (x.gt.elim) go to 330
      coef = coef*exp(-x)
  130 continue
      if (abs(dnu).eq.0.5d0) go to 340
      if (x.gt.x2) go to 280
c
c     miller algorithm for x1.lt.x.le.x2
c
      etest = cos(pi*dnu)/(pi*x*tol)
      fks = 1.0d0
      fhs = 0.25d0
      fk = 0.0d0
      ck = x + x + 2.0d0
      p1 = 0.0d0
      p2 = 1.0d0
      k = 0
  140 continue
      k = k + 1
      fk = fk + 1.0d0
      ak = (fhs-dnu2)/(fks+fk)
      bk = ck/(fk+1.0d0)
      pt = p2
      p2 = bk*p2 - ak*p1
      p1 = pt
      a(k) = ak
      b(k) = bk
      ck = ck + 2.0d0
      fks = fks + fk + fk + 1.0d0
      fhs = fhs + fk + fk
      if (etest.gt.fk*p1) go to 140
      kk = k
      s = 1.0d0
      p1 = 0.0d0
      p2 = 1.0d0
      do 150 i=1,k
        pt = p2
        p2 = (b(kk)*p2-p1)/a(kk)
        p1 = pt
        s = s + p2
        kk = kk - 1
  150 continue
      s1 = coef*(p2/s)
      if (inu.gt.0 .or. n.gt.1) go to 160
      go to 200
  160 continue
      s2 = s1*(x+dnu+0.5d0-p1/p2)/x
c
c     forward recursion on the three term recursion relation
c
  170 continue
      ck = (dnu+dnu+2.0d0)/x
      if (n.eq.1) inu = inu - 1
      if (inu.gt.0) go to 180
      if (n.gt.1) go to 200
      s1 = s2
      go to 200
  180 continue
      do 190 i=1,inu
        st = s2
        s2 = ck*s2 + s1
        s1 = st
        ck = ck + rx
  190 continue
      if (n.eq.1) s1 = s2
  200 continue
      if (iflag.eq.1) go to 220
      y(1) = s1
      if (n.eq.1) return
      y(2) = s2
      if (n.eq.2) return
      do 210 i=3,n
        y(i) = ck*y(i-1) + y(i-2)
        ck = ck + rx
  210 continue
      return
c     iflag=1 cases
  220 continue
      s = -x + log(s1)
      y(1) = 0.0d0
      nz = 1
      if (s.lt.-elim) go to 230
      y(1) = exp(s)
      nz = 0
  230 continue
      if (n.eq.1) return
      s = -x + log(s2)
      y(2) = 0.0d0
      nz = nz + 1
      if (s.lt.-elim) go to 240
      nz = nz - 1
      y(2) = exp(s)
  240 continue
      if (n.eq.2) return
      kk = 2
      if (nz.lt.2) go to 260
      do 250 i=3,n
        kk = i
        st = s2
        s2 = ck*s2 + s1
        s1 = st
        ck = ck + rx
        s = -x + log(s2)
        nz = nz + 1
        y(i) = 0.0d0
        if (s.lt.-elim) go to 250
        y(i) = exp(s)
        nz = nz - 1
        go to 260
  250 continue
      return
  260 continue
      if (kk.eq.n) return
      s2 = s2*ck + s1
      ck = ck + rx
      kk = kk + 1
      y(kk) = exp(-x+log(s2))
      if (kk.eq.n) return
      kk = kk + 1
      do 270 i=kk,n
        y(i) = ck*y(i-1) + y(i-2)
        ck = ck + rx
  270 continue
      return
c
c     asymptotic expansion for large x, x.gt.x2
c
c     iflag=0 means no underflow occurred
c     iflag=1 means an underflow occurred- computation proceeds with
c     koded=2 and a test for on scale values is made during forward
c     recursion
  280 continue
      nn = 2
      if (inu.eq.0 .and. n.eq.1) nn = 1
      dnu2 = dnu + dnu
      fmu = 0.0d0
      if (abs(dnu2).lt.tol) go to 290
      fmu = dnu2*dnu2
  290 continue
      ex = x*8.0d0
      s2 = 0.0d0
      do 320 k=1,nn
        s1 = s2
        s = 1.0d0
        ak = 0.0d0
        ck = 1.0d0
        sqk = 1.0d0
        dk = ex
        do 300 j=1,30
          ck = ck*(fmu-sqk)/dk
          s = s + ck
          dk = dk + ex
          ak = ak + 8.0d0
          sqk = sqk + ak
          if (abs(ck).lt.tol) go to 310
  300   continue
  310   s2 = s*coef
        fmu = fmu + 8.0d0*dnu + 4.0d0
  320 continue
      if (nn.gt.1) go to 170
      s1 = s2
      go to 200
  330 continue
      koded = 2
      iflag = 1
      go to 120
c
c     fnu=half odd integer case
c
  340 continue
      s1 = coef
      s2 = coef
      go to 170
c
c
  350 call xermsg ('slatec', 'dbsknu', 'x not greater than zero', 2, 1)
      return
  360 call xermsg ('slatec', 'dbsknu', 'fnu not zero or positive', 2,
     +   1)
      return
  370 call xermsg ('slatec', 'dbsknu', 'kode not 1 or 2', 2, 1)
      return
  380 call xermsg ('slatec', 'dbsknu', 'n not greater than 0', 2, 1)
      return
      end
