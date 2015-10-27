*deck dexint
      subroutine dexint (x, n, kode, m, tol, en, nz, ierr)
c***begin prologue  dexint
c***purpose  compute an m member sequence of exponential integrals
c            e(n+k,x), k=0,1,...,m-1 for n .ge. 1 and x .ge. 0.
c***library   slatec
c***category  c5
c***type      double precision (exint-s, dexint-d)
c***keywords  exponential integral, special functions
c***author  amos, d. e., (snla)
c***description
c
c         dexint computes m member sequences of exponential integrals
c         e(n+k,x), k=0,1,...,m-1 for n .ge. 1 and x .ge. 0.  the
c         exponential integral is defined by
c
c         e(n,x)=integral on (1,infinity) of exp(-xt)/t**n
c
c         where x=0.0 and n=1 cannot occur simultaneously.  formulas
c         and notation are found in the nbs handbook of mathematical
c         functions (ref. 1).
c
c         the power series is implemented for x .le. xcut and the
c         confluent hypergeometric representation
c
c                     e(a,x) = exp(-x)*(x**(a-1))*u(a,a,x)
c
c         is computed for x .gt. xcut.  since sequences are computed in
c         a stable fashion by recurring away from x, a is selected as
c         the integer closest to x within the constraint n .le. a .le.
c         n+m-1.  for the u computation, a is further modified to be the
c         nearest even integer.  indices are carried forward or
c         backward by the two term recursion relation
c
c                     k*e(k+1,x) + x*e(k,x) = exp(-x)
c
c         once e(a,x) is computed.  the u function is computed by means
c         of the backward recursive miller algorithm applied to the
c         three term contiguous relation for u(a+k,a,x), k=0,1,...
c         this produces accurate ratios and determines u(a+k,a,x), and
c         hence e(a,x), to within a multiplicative constant c.
c         another contiguous relation applied to c*u(a,a,x) and
c         c*u(a+1,a,x) gets c*u(a+1,a+1,x), a quantity proportional to
c         e(a+1,x).  the normalizing constant c is obtained from the
c         two term recursion relation above with k=a.
c
c         the maximum number of significant digits obtainable
c         is the smaller of 14 and the number of digits carried in
c         double precision arithmetic.
c
c     description of arguments
c
c         input     * x and tol are double precision *
c           x       x .gt. 0.0 for n=1 and  x .ge. 0.0 for n .ge. 2
c           n       order of the first member of the sequence, n .ge. 1
c                   (x=0.0 and n=1 is an error)
c           kode    a selection parameter for scaled values
c                   kode=1   returns        e(n+k,x), k=0,1,...,m-1.
c                       =2   returns exp(x)*e(n+k,x), k=0,1,...,m-1.
c           m       number of exponential integrals in the sequence,
c                   m .ge. 1
c           tol     relative accuracy wanted, etol .le. tol .le. 0.1
c                   etol is the larger of double precision unit
c                   roundoff = d1mach(4) and 1.0d-18
c
c         output    * en is a double precision vector *
c           en      a vector of dimension at least m containing values
c                   en(k) = e(n+k-1,x) or exp(x)*e(n+k-1,x), k=1,m
c                   depending on kode
c           nz      underflow indicator
c                   nz=0   a normal return
c                   nz=m   x exceeds xlim and an underflow occurs.
c                          en(k)=0.0d0 , k=1,m returned on kode=1
c           ierr    error flag
c                   ierr=0, normal return, computation completed
c                   ierr=1, input error,   no computation
c                   ierr=2, error,         no computation
c                           algorithm termination condition not met
c
c***references  m. abramowitz and i. a. stegun, handbook of
c                 mathematical functions, nbs ams series 55, u.s. dept.
c                 of commerce, 1955.
c               d. e. amos, computation of exponential integrals, acm
c                 transactions on mathematical software 6, (1980),
c                 pp. 365-377 and pp. 420-428.
c***routines called  d1mach, dpsixn, i1mach
c***revision history  (yymmdd)
c   800501  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900326  removed duplicate information from description section.
c           (wrb)
c   910408  updated the references section.  (wrb)
c   920207  updated with code with a revision date of 880811 from
c           d. amos.  included correction of argument list.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dexint
      double precision a,aa,aams,ah,ak,at,b,bk,bt,cc,cnorm,ct,em,emx,en,
     1                 etol,fnm,fx,pt,p1,p2,s,tol,tx,x,xcut,xlim,xtol,y,
     2                 yt,y1,y2
      double precision d1mach,dpsixn
      integer i,ic,icase,ict,ierr,ik,ind,ix,i1m,jset,k,kk,kn,kode,ks,m,
     1        ml,mu,n,nd,nm,nz
      integer i1mach
      dimension en(*), a(99), b(99), y(2)
      save xcut
      data xcut / 2.0d0 /
c***first executable statement  dexint
      ierr = 0
      nz = 0
      etol = max(d1mach(4),0.5d-18)
      if (x.lt.0.0d0) ierr = 1
      if (n.lt.1) ierr = 1
      if (kode.lt.1 .or. kode.gt.2) ierr = 1
      if (m.lt.1) ierr = 1
      if (tol.lt.etol .or. tol.gt.0.1d0) ierr = 1
      if (x.eq.0.0d0 .and. n.eq.1) ierr = 1
      if(ierr.ne.0) return
      i1m = -i1mach(15)
      pt = 2.3026d0*i1m*d1mach(5)
      xlim = pt - 6.907755d0
      bt = pt + (n+m-1)
      if (bt.gt.1000.0d0) xlim = pt - log(bt)
c
      if (x.gt.xcut) go to 100
      if (x.eq.0.0d0 .and. n.gt.1) go to 80
c-----------------------------------------------------------------------
c     series for e(n,x) for x.le.xcut
c-----------------------------------------------------------------------
      tx = x + 0.5d0
      ix = tx
c-----------------------------------------------------------------------
c     icase=1 means integer closest to x is 2 and n=1
c     icase=2 means integer closest to x is 0,1, or 2 and n.ge.2
c-----------------------------------------------------------------------
      icase = 2
      if (ix.gt.n) icase = 1
      nm = n - icase + 1
      nd = nm + 1
      ind = 3 - icase
      mu = m - ind
      ml = 1
      ks = nd
      fnm = nm
      s = 0.0d0
      xtol = 3.0d0*tol
      if (nd.eq.1) go to 10
      xtol = 0.3333d0*tol
      s = 1.0d0/fnm
   10 continue
      aa = 1.0d0
      ak = 1.0d0
      ic = 35
      if (x.lt.etol) ic = 1
      do 50 i=1,ic
        aa = -aa*x/ak
        if (i.eq.nm) go to 30
        s = s - aa/(ak-fnm)
        if (abs(aa).le.xtol*abs(s)) go to 20
        ak = ak + 1.0d0
        go to 50
   20   continue
        if (i.lt.2) go to 40
        if (nd-2.gt.i .or. i.gt.nd-1) go to 60
        ak = ak + 1.0d0
        go to 50
   30   s = s + aa*(-log(x)+dpsixn(nd))
        xtol = 3.0d0*tol
   40   ak = ak + 1.0d0
   50 continue
      if (ic.ne.1) go to 340
   60 if (nd.eq.1) s = s + (-log(x)+dpsixn(1))
      if (kode.eq.2) s = s*exp(x)
      en(1) = s
      emx = 1.0d0
      if (m.eq.1) go to 70
      en(ind) = s
      aa = ks
      if (kode.eq.1) emx = exp(-x)
      go to (220, 240), icase
   70 if (icase.eq.2) return
      if (kode.eq.1) emx = exp(-x)
      en(1) = (emx-s)/x
      return
   80 continue
      do 90 i=1,m
        en(i) = 1.0d0/(n+i-2)
   90 continue
      return
c-----------------------------------------------------------------------
c     backward recursive miller algorithm for
c              e(n,x)=exp(-x)*(x**(n-1))*u(n,n,x)
c     with recursion away from n=integer closest to x.
c     u(a,b,x) is the second confluent hypergeometric function
c-----------------------------------------------------------------------
  100 continue
      emx = 1.0d0
      if (kode.eq.2) go to 130
      if (x.le.xlim) go to 120
      nz = m
      do 110 i=1,m
        en(i) = 0.0d0
  110 continue
      return
  120 emx = exp(-x)
  130 continue
      tx = x + 0.5d0
      ix = tx
      kn = n + m - 1
      if (kn.le.ix) go to 140
      if (n.lt.ix .and. ix.lt.kn) go to 170
      if (n.ge.ix) go to 160
      go to 340
  140 icase = 1
      ks = kn
      ml = m - 1
      mu = -1
      ind = m
      if (kn.gt.1) go to 180
  150 ks = 2
      icase = 3
      go to 180
  160 icase = 2
      ind = 1
      ks = n
      mu = m - 1
      if (n.gt.1) go to 180
      if (kn.eq.1) go to 150
      ix = 2
  170 icase = 1
      ks = ix
      ml = ix - n
      ind = ml + 1
      mu = kn - ix
  180 continue
      ik = ks/2
      ah = ik
      jset = 1 + ks - (ik+ik)
c-----------------------------------------------------------------------
c     start computation for
c              en(ind) = c*u( a , a ,x)    jset=1
c              en(ind) = c*u(a+1,a+1,x)    jset=2
c     for an even integer a.
c-----------------------------------------------------------------------
      ic = 0
      aa = ah + ah
      aams = aa - 1.0d0
      aams = aams*aams
      tx = x + x
      fx = tx + tx
      ak = ah
      xtol = tol
      if (tol.le.1.0d-3) xtol = 20.0d0*tol
      ct = aams + fx*ah
      em = (ah+1.0d0)/((x+aa)*xtol*sqrt(ct))
      bk = aa
      cc = ah*ah
c-----------------------------------------------------------------------
c     forward recursion for p(ic),p(ic+1) and index ic for backward
c     recursion
c-----------------------------------------------------------------------
      p1 = 0.0d0
      p2 = 1.0d0
  190 continue
      if (ic.eq.99) go to 340
      ic = ic + 1
      ak = ak + 1.0d0
      at = bk/(bk+ak+cc+ic)
      bk = bk + ak + ak
      a(ic) = at
      bt = (ak+ak+x)/(ak+1.0d0)
      b(ic) = bt
      pt = p2
      p2 = bt*p2 - at*p1
      p1 = pt
      ct = ct + fx
      em = em*at*(1.0d0-tx/ct)
      if (em*(ak+1.0d0).gt.p1*p1) go to 190
      ict = ic
      kk = ic + 1
      bt = tx/(ct+fx)
      y2 = (bk/(bk+cc+kk))*(p1/p2)*(1.0d0-bt+0.375d0*bt*bt)
      y1 = 1.0d0
c-----------------------------------------------------------------------
c     backward recurrence for
c              y1=             c*u( a ,a,x)
c              y2= c*(a/(1+a/2))*u(a+1,a,x)
c-----------------------------------------------------------------------
      do 200 k=1,ict
        kk = kk - 1
        yt = y1
        y1 = (b(kk)*y1-y2)/a(kk)
        y2 = yt
  200 continue
c-----------------------------------------------------------------------
c     the contiguous relation
c              x*u(b,c+1,x)=(c-b)*u(b,c,x)+u(b-1,c,x)
c     with  b=a+1 , c=a is used for
c              y(2) = c * u(a+1,a+1,x)
c     x is incorporated into the normalizing relation
c-----------------------------------------------------------------------
      pt = y2/y1
      cnorm = 1.0e0 - pt*(ah+1.0e0)/aa
      y(1) = 1.0e0/(cnorm*aa+x)
      y(2) = cnorm*y(1)
      if (icase.eq.3) go to 210
      en(ind) = emx*y(jset)
      if (m.eq.1) return
      aa = ks
      go to (220, 240), icase
c-----------------------------------------------------------------------
c     recursion section  n*e(n+1,x) + x*e(n,x)=emx
c-----------------------------------------------------------------------
  210 en(1) = emx*(1.0e0-y(1))/x
      return
  220 k = ind - 1
      do 230 i=1,ml
        aa = aa - 1.0d0
        en(k) = (emx-aa*en(k+1))/x
        k = k - 1
  230 continue
      if (mu.le.0) return
      aa = ks
  240 k = ind
      do 250 i=1,mu
        en(k+1) = (emx-x*en(k))/aa
        aa = aa + 1.0d0
        k = k + 1
  250 continue
      return
  340 continue
      ierr = 2
      return
      end
