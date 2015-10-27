*deck besk
      subroutine besk (x, fnu, kode, n, y, nz)
c***begin prologue  besk
c***purpose  implement forward recursion on the three term recursion
c            relation for a sequence of non-negative order bessel
c            functions k/sub(fnu+i-1)/(x), or scaled bessel functions
c            exp(x)*k/sub(fnu+i-1)/(x), i=1,...,n for real, positive
c            x and non-negative orders fnu.
c***library   slatec
c***category  c10b3
c***type      single precision (besk-s, dbesk-d)
c***keywords  k bessel function, special functions
c***author  amos, d. e., (snla)
c***description
c
c     abstract
c         besk implements forward recursion on the three term
c         recursion relation for a sequence of non-negative order bessel
c         functions k/sub(fnu+i-1)/(x), or scaled bessel functions
c         exp(x)*k/sub(fnu+i-1)/(x), i=1,...,n for real x .gt. 0.0e0 and
c         non-negative orders fnu.  if fnu .lt. nulim, orders fnu and
c         fnu+1 are obtained from besknu to start the recursion.  if
c         fnu .ge. nulim, the uniform asymptotic expansion is used for
c         orders fnu and fnu+1 to start the recursion.  nulim is 35 or
c         70 depending on whether n=1 or n .ge. 2.  under and overflow
c         tests are made on the leading term of the asymptotic expansion
c         before any extensive computation is done.
c
c     description of arguments
c
c         input
c           x      - x .gt. 0.0e0
c           fnu    - order of the initial k function, fnu .ge. 0.0e0
c           kode   - a parameter to indicate the scaling option
c                    kode=1 returns y(i)=       k/sub(fnu+i-1)/(x),
c                                        i=1,...,n
c                    kode=2 returns y(i)=exp(x)*k/sub(fnu+i-1)/(x),
c                                        i=1,...,n
c           n      - number of members in the sequence, n .ge. 1
c
c         output
c           y      - a vector whose first n components contain values
c                    for the sequence
c                    y(i)=       k/sub(fnu+i-1)/(x), i=1,...,n  or
c                    y(i)=exp(x)*k/sub(fnu+i-1)/(x), i=1,...,n
c                    depending on kode
c           nz     - number of components of y set to zero due to
c                    underflow with kode=1,
c                    nz=0   , normal return, computation completed
c                    nz .ne. 0, first nz components of y set to zero
c                             due to underflow, y(i)=0.0e0, i=1,...,nz
c
c     error conditions
c         improper input arguments - a fatal error
c         overflow - a fatal error
c         underflow with kode=1 -  a non-fatal error (nz .ne. 0)
c
c***references  f. w. j. olver, tables of bessel functions of moderate
c                 or large orders, npl mathematical tables 6, her
c                 majesty's stationery office, london, 1962.
c               n. m. temme, on the numerical evaluation of the modified
c                 bessel function of the third kind, journal of
c                 computational physics 19, (1975), pp. 324-337.
c***routines called  asyik, besk0, besk0e, besk1, besk1e, besknu,
c                    i1mach, r1mach, xermsg
c***revision history  (yymmdd)
c   790201  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  besk
c
      integer i, j, k, kode, mz, n, nb, nd, nn, nud, nulim, nz
      integer i1mach
      real cn, dnu, elim, etx, flgik,fn, fnn, fnu,gln,gnu,rtz,s,s1,s2,
     1 t, tm, trx, w, x, xlim, y, zn
      real besk0, besk1, besk1e, besk0e, r1mach
      dimension w(2), nulim(2), y(*)
      save nulim
      data nulim(1),nulim(2) / 35 , 70 /
c***first executable statement  besk
      nn = -i1mach(12)
      elim = 2.303e0*(nn*r1mach(5)-3.0e0)
      xlim = r1mach(1)*1.0e+3
      if (kode.lt.1 .or. kode.gt.2) go to 280
      if (fnu.lt.0.0e0) go to 290
      if (x.le.0.0e0) go to 300
      if (x.lt.xlim) go to 320
      if (n.lt.1) go to 310
      etx = kode - 1
c
c     nd is a dummy variable for n
c     gnu is a dummy variable for fnu
c     nz = number of underflows on kode=1
c
      nd = n
      nz = 0
      nud = int(fnu)
      dnu = fnu - nud
      gnu = fnu
      nn = min(2,nd)
      fn = fnu + n - 1
      fnn = fn
      if (fn.lt.2.0e0) go to 150
c
c     overflow test  (leading exponential of asymptotic expansion)
c     for the last order, fnu+n-1.ge.nulim
c
      zn = x/fn
      if (zn.eq.0.0e0) go to 320
      rtz = sqrt(1.0e0+zn*zn)
      gln = log((1.0e0+rtz)/zn)
      t = rtz*(1.0e0-etx) + etx/(zn+rtz)
      cn = -fn*(t-gln)
      if (cn.gt.elim) go to 320
      if (nud.lt.nulim(nn)) go to 30
      if (nn.eq.1) go to 20
   10 continue
c
c     underflow test (leading exponential of asymptotic expansion)
c     for the first order, fnu.ge.nulim
c
      fn = gnu
      zn = x/fn
      rtz = sqrt(1.0e0+zn*zn)
      gln = log((1.0e0+rtz)/zn)
      t = rtz*(1.0e0-etx) + etx/(zn+rtz)
      cn = -fn*(t-gln)
   20 continue
      if (cn.lt.-elim) go to 230
c
c     asymptotic expansion for orders fnu and fnu+1.ge.nulim
c
      flgik = -1.0e0
      call asyik(x,gnu,kode,flgik,rtz,cn,nn,y)
      if (nn.eq.1) go to 240
      trx = 2.0e0/x
      tm = (gnu+gnu+2.0e0)/x
      go to 130
c
   30 continue
      if (kode.eq.2) go to 40
c
c     underflow test (leading exponential of asymptotic expansion in x)
c     for order dnu
c
      if (x.gt.elim) go to 230
   40 continue
      if (dnu.ne.0.0e0) go to 80
      if (kode.eq.2) go to 50
      s1 = besk0(x)
      go to 60
   50 s1 = besk0e(x)
   60 continue
      if (nud.eq.0 .and. nd.eq.1) go to 120
      if (kode.eq.2) go to 70
      s2 = besk1(x)
      go to 90
   70 s2 = besk1e(x)
      go to 90
   80 continue
      nb = 2
      if (nud.eq.0 .and. nd.eq.1) nb = 1
      call besknu(x, dnu, kode, nb, w, nz)
      s1 = w(1)
      if (nb.eq.1) go to 120
      s2 = w(2)
   90 continue
      trx = 2.0e0/x
      tm = (dnu+dnu+2.0e0)/x
c     forward recur from dnu to fnu+1 to get y(1) and y(2)
      if (nd.eq.1) nud = nud - 1
      if (nud.gt.0) go to 100
      if (nd.gt.1) go to 120
      s1 = s2
      go to 120
  100 continue
      do 110 i=1,nud
        s = s2
        s2 = tm*s2 + s1
        s1 = s
        tm = tm + trx
  110 continue
      if (nd.eq.1) s1 = s2
  120 continue
      y(1) = s1
      if (nd.eq.1) go to 240
      y(2) = s2
  130 continue
      if (nd.eq.2) go to 240
c     forward recur from fnu+2 to fnu+n-1
      do 140 i=3,nd
        y(i) = tm*y(i-1) + y(i-2)
        tm = tm + trx
  140 continue
      go to 240
c
  150 continue
c     underflow test for kode=1
      if (kode.eq.2) go to 160
      if (x.gt.elim) go to 230
  160 continue
c     overflow test
      if (fn.le.1.0e0) go to 170
      if (-fn*(log(x)-0.693e0).gt.elim) go to 320
  170 continue
      if (dnu.eq.0.0e0) go to 180
      call besknu(x, fnu, kode, nd, y, mz)
      go to 240
  180 continue
      j = nud
      if (j.eq.1) go to 210
      j = j + 1
      if (kode.eq.2) go to 190
      y(j) = besk0(x)
      go to 200
  190 y(j) = besk0e(x)
  200 if (nd.eq.1) go to 240
      j = j + 1
  210 if (kode.eq.2) go to 220
      y(j) = besk1(x)
      go to 240
  220 y(j) = besk1e(x)
      go to 240
c
c     update parameters on underflow
c
  230 continue
      nud = nud + 1
      nd = nd - 1
      if (nd.eq.0) go to 240
      nn = min(2,nd)
      gnu = gnu + 1.0e0
      if (fnn.lt.2.0e0) go to 230
      if (nud.lt.nulim(nn)) go to 230
      go to 10
  240 continue
      nz = n - nd
      if (nz.eq.0) return
      if (nd.eq.0) go to 260
      do 250 i=1,nd
        j = n - i + 1
        k = nd - i + 1
        y(j) = y(k)
  250 continue
  260 continue
      do 270 i=1,nz
        y(i) = 0.0e0
  270 continue
      return
c
c
c
  280 continue
      call xermsg ('slatec', 'besk', 'scaling option, kode, not 1 or 2'
     +   , 2, 1)
      return
  290 continue
      call xermsg ('slatec', 'besk', 'order, fnu, less than zero', 2,
     +   1)
      return
  300 continue
      call xermsg ('slatec', 'besk', 'x less than or equal to zero', 2,
     +   1)
      return
  310 continue
      call xermsg ('slatec', 'besk', 'n less than one', 2, 1)
      return
  320 continue
      call xermsg ('slatec', 'besk',
     +   'overflow, fnu or n too large or x too small', 6, 1)
      return
      end
